# ===============================================================================
# REACH Secondary Analysis: Baseline Mortality Estimation
# ===============================================================================
# Purpose: Estimate community-level baseline U5MR and IMR using spatial hierarchical models
# Author: William Msemburi
# Date: August 2025
# ===============================================================================

rm(list = ls())

# Working directory setup
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

# Load required libraries
required_packages <- c(
  "dplyr", "tidyr", "tibble", "purrr", "sf", "ggokabeito",
  "cmdstanr", "posterior", "bayesplot", "ggplot2", "showtext", "INLA"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed.", pkg))
  }
}

# Load utility functions
source("../utils/reach_analysis_utils.R")

# ===============================================================================
# DATA PREPARATION
# ===============================================================================

cat("Loading and preparing data...\n")

# Load inputs
# 1) Trial microdata (cleaned with counts/exposure)
if (!file.exists("../../data/Clean/data_crt_geospatial.rda")) {
  stop("Clean trial data not found at ../../data/Clean/data_crt_geospatial.rda")
}
load("../../data/Clean/data_crt_geospatial.rda")

# 2) DHS prior: trial-country hazard priors
if (!file.exists("../../results/trial_hazard_priors.rds")) {
  stop("Trial hazard priors not found. Run derive_trial_priors.R first.")
}
final_priors <- readRDS("../../results/trial_hazard_priors.rds")

# 3) DHS pooled neonatal regression
if (!file.exists("../../results/dhs_h0_pooled_model.rds")) {
  stop("DHS H0 pooled model not found. Run derive_trial_priors.R first.")
}
pooled_model <- readRDS("../../results/dhs_h0_pooled_model.rds")
stopifnot(all(c("coef","vcov","terms","sigma_hat") %in% names(pooled_model)))

# ===============================================================================
# STEP 1: Aggregate trial data to 1–59m counts/exposure per cluster
# ===============================================================================

cat("Aggregating trial data...\n")

agg_159 <- data_crt_geospatial %>%
  dplyr::filter(!is.na(trt), !is.na(person_time),
                !is.na(latitude), !is.na(longitude),
                age_start >= 1, age_start < 60) %>%
  dplyr::mutate(
    # Recode before composing trial_country; keeps Niger I/II together
    trial_phase = dplyr::case_when(
      country == "Niger" & trial_phase != "AVENIR" ~ "MORDOR I/II",
      TRUE ~ trial_phase
    ),
    E_years = person_time / 365.25
  ) %>%
  dplyr::group_by(community_id, trial_phase, country, trt, latitude, longitude) %>%
  dplyr::summarise(y_159 = sum(n), E_yrs = sum(E_years), .groups = "drop")

# Create stable indices
group_index   <- create_group_index(agg_159)
country_index <- create_country_index(agg_159)

df.mod <- agg_159 %>%
  dplyr::mutate(group = paste(trial_phase, country, sep = " ")) %>%
  dplyr::inner_join(group_index,  by = c("group","trial_phase","country")) %>%
  dplyr::inner_join(country_index, by = "country") %>%
  dplyr::arrange(group_id, community_id) %>%
  dplyr::mutate(
    row_id = row_number(),
    is_obs = as.integer(trt == 0L),                       # placebo-only contributes to likelihood
    logE   = log(pmax(E_yrs, 1e-8))
  )

cat("Data aggregated: ", nrow(df.mod), " clusters across ", 
    dplyr::n_distinct(df.mod$group_id), " trial-countries\n")

# ===============================================================================
# STEP 2: Group-level prior for total 1–59m log-rate from DHS H1_59 priors
# ===============================================================================

cat("Setting up group-level priors...\n")

gp <- final_priors %>%
  dplyr::mutate(group = trial_country) %>%
  dplyr::inner_join(group_index, by = c("group" = "group")) %>%
  dplyr::arrange(group_id) %>%
  dplyr::transmute(
    group_id,
    mu_prior = log_h1_59_prior,                    # log cumulative hazard H1_59
    sd_prior = pmax(log_h1_59_prior_sd, 1e-6)
  )

stopifnot(nrow(gp) == dplyr::n_distinct(df.mod$group_id))

# Convert cumulative hazard over 1–59m to per-year log-rate:
#   H1_59 = exp(log_lambda_total) * (59/12)  ⇒ log_lambda_total = log(H1_59) - log(59/12)
log_59_12 <- log(59/12)
gp <- gp %>%
  dplyr::mutate(
    mu_logLambda_prior = mu_prior - log_59_12,
    sd_logLambda_prior = sd_prior
  )

# ===============================================================================
# STEP 3: Age-split prior p_g ~ Beta(a0,b0) from DHS cumulative hazards
# ===============================================================================

cat("Deriving age-split priors...\n")

set.seed(1)
S <- 50000L

beta_prior_by_group <- final_priors %>%
  dplyr::transmute(
    trial_country,
    p_draws = purrr::pmap(
      list(log_h1_11_prior, log_h1_11_prior_sd,
           log_h12_59_prior, log_h12_59_prior_sd),
      function(m1, s1, m2, s2) {
        H111  <- rlnorm(S, meanlog = m1, sdlog = s1)      # cumulative hazards
        H1259 <- rlnorm(S, meanlog = m2, sdlog = s2)
        H111 / (H111 + H1259)                             # share of H1–59 in 1–11
      }
    )
  ) %>%
  dplyr::mutate(
    m = purrr::map_dbl(p_draws, mean),
    v = purrr::map_dbl(p_draws, var),
    v = pmin(v, 0.999 * m * (1 - m)),                    # Beta-compatible variance
    k = m * (1 - m) / v - 1,
    a0 = m * k,
    b0 = (1 - m) * k
  ) %>%
  dplyr::select(trial_country, a0, b0)

# Map to group_id order
beta_ab_by_gid <- beta_prior_by_group %>%
  dplyr::mutate(group = trial_country) %>%
  dplyr::inner_join(group_index, by = c("group" = "group")) %>%
  dplyr::arrange(group_id) %>%
  dplyr::transmute(a0 = a0, b0 = b0)

# ===============================================================================
# STEP 4: DHS pooled neonatal regression -> country-specific parameter draws
# ===============================================================================

cat("Setting up DHS neonatal bridge...\n")

nm     <- pooled_model$terms
coefs  <- pooled_model$coef
V      <- pooled_model$vcov

i_b0    <- match("(Intercept)", nm)
i_b111  <- match("log(H_111)",  nm)
i_b1259 <- match("log(H_1259)", nm)
stopifnot(all(!is.na(c(i_b0, i_b111, i_b1259))))

fe_name <- function(cn) paste0("factor(country)", cn)
trial_countries <- as.character(country_index$country)  # in country_id order
Kc <- length(trial_countries)

alpha_row <- function(cn){
  A <- rep(0, length(nm)); A[i_b0] <- 1
  j <- match(fe_name(cn), nm, nomatch = 0)
  if (j > 0) A[j] <- 1
  alpha_c <- sum(A * coefs)
  list(alpha = alpha_c, A = A)
}

ab_mean <- matrix(NA_real_, nrow = Kc, ncol = 3)
L_ab    <- array(NA_real_, dim = c(Kc, 3, 3))

for (k in seq_len(Kc)) {
  cn <- trial_countries[k]
  ax <- alpha_row(cn); A <- ax$A
  
  # mean vector: [alpha_c, beta_logH111, beta_logH1259]
  ab_mean[k, ] <- c(ax$alpha, coefs[i_b111], coefs[i_b1259])
  
  # covariance via linear combos from full VCOV
  var_a     <- as.numeric(t(A) %*% V %*% A)
  cov_a_b1  <- as.numeric(t(A) %*% V[, i_b111])
  cov_a_b2  <- as.numeric(t(A) %*% V[, i_b1259])
  var_b1    <- V[i_b111,  i_b111]
  var_b2    <- V[i_b1259, i_b1259]
  cov_b1_b2 <- V[i_b111,  i_b1259]
  
  S <- matrix(c(
    var_a,     cov_a_b1,  cov_a_b2,
    cov_a_b1,  var_b1,    cov_b1_b2,
    cov_a_b2,  cov_b1_b2, var_b2
  ), nrow = 3, byrow = TRUE)
  
  # numerical safety
  S <- S + 1e-8 * diag(3)
  L_ab[k, , ] <- t(chol(S))   # lower Cholesky for Stan
}
sigma0_hat <- rep(pooled_model$sigma_hat, Kc)  # fixed residual scale per country

# ===============================================================================
# STEP 5: HSGP basis per trial-country group
# ===============================================================================

cat("Computing HSGP spatial basis...\n")

# HSGP computation function
compute_group_hsgp <- function(df_g, m1 = 25, m2 = 25, c_factor = 1.1) {
  n_g <- nrow(df_g)
  if (n_g <= 1) return(list(
    PHI = matrix(0, nrow = n_g, ncol = 1),
    lambda_hsgp = array(1), m1 = 1L, m2 = 1L, L1 = 1, L2 = 1
  ))
  
  tryCatch({
    epsg <- pick_utm_epsg(mean(df_g$longitude), mean(df_g$latitude))
    sf_points <- sf::st_as_sf(df_g, coords = c("longitude","latitude"), crs = 4326)
    xy_m <- sf::st_coordinates(sf::st_transform(sf_points, epsg))
    xr <- range(xy_m[,1]); yr <- range(xy_m[,2])
    x1_span <- max(diff(xr), 1000); x2_span <- max(diff(yr), 1000)
    x1s <- (xy_m[,1] - xr[1]) / x1_span; x2s <- (xy_m[,2] - yr[1]) / x2_span
    L1 <- max(c_factor * max(x1s), 1.0); L2 <- max(c_factor * max(x2s), 1.0)
    m1 <- max(5L, as.integer(m1)); m2 <- max(5L, as.integer(m2))
    
    lambda_hsgp <- numeric(m1*m2)
    PHI <- matrix(0, nrow = n_g, ncol = m1*m2)
    idx <- 1L
    for (i in 1:m1) for (j in 1:m2) {
      lambda_hsgp[idx] <- (i*pi/(2*L1))^2 + (j*pi/(2*L2))^2  # ω^2
      PHI[, idx] <- sin(i*pi*x1s/(2*L1)) * sin(j*pi*x2s/(2*L2))
      idx <- idx + 1L
    }
    PHI <- PHI * sqrt(4/(L1*L2))  # orthonormal scaling
    list(PHI = PHI, lambda_hsgp = lambda_hsgp,
         m1 = as.integer(m1), m2 = as.integer(m2), L1 = L1, L2 = L2)
  }, error = function(e){
    warning("HSGP basis failed, returning small random basis: ", e$message)
    list(PHI = matrix(rnorm(n_g*4, 0, 0.1), nrow = n_g, ncol = 4),
         lambda_hsgp = rep(1, 4), m1 = 2L, m2 = 2L, L1 = 1, L2 = 1)
  })
}

split_idx <- split(seq_len(nrow(df.mod)), df.mod$group_id)
groups <- sort(unique(df.mod$group_id))

hsgp_results <- lapply(seq_along(groups), function(h) {
  compute_group_hsgp(df.mod[split_idx[[h]], ], m1 = 25, m2 = 25, c_factor = 1.1)
})

# Combine HSGP results
total_basis <- sum(sapply(hsgp_results, \(x) ncol(x$PHI)))
C <- nrow(df.mod)
PHI_combined <- matrix(0, nrow = C, ncol = total_basis)
lambda_combined <- numeric(total_basis)
group_basis_start <- integer(length(groups))
group_basis_size  <- integer(length(groups))

col_start <- 1L
for (h in seq_along(groups)) {
  rows <- split_idx[[h]]
  nb <- ncol(hsgp_results[[h]]$PHI)
  group_basis_start[h] <- col_start
  group_basis_size[h]  <- nb
  PHI_combined[rows, col_start:(col_start+nb-1)] <- hsgp_results[[h]]$PHI
  lambda_combined[col_start:(col_start+nb-1)]     <- hsgp_results[[h]]$lambda_hsgp
  col_start <- col_start + nb
}

cat("HSGP basis constructed: ", total_basis, " basis functions\n")

# ===============================================================================
# STEP 6: Stan data preparation
# ===============================================================================

idx_obs <- which(df.mod$is_obs == 1L); stopifnot(length(idx_obs) > 0)

stan_data <- list(
  C = nrow(df.mod),
  K = dplyr::n_distinct(df.mod$group_id),
  Kc = Kc,
  group_id = as.integer(df.mod$group_id),
  country_id = as.integer(df.mod$country_id),
  y = as.integer(df.mod$y_159),
  logE = as.vector(df.mod$logE),
  is_obs = as.integer(df.mod$is_obs),
  C_obs = length(idx_obs),
  idx_obs = as.array(idx_obs),
  
  # HSGP
  M_hsgp = as.integer(total_basis),
  PHI = PHI_combined,
  lambda_hsgp = lambda_combined,
  n_groups = length(groups),
  group_basis_start = as.array(group_basis_start),
  group_basis_size  = as.array(group_basis_size),
  
  # DHS neonatal bridge
  ab_mean = ab_mean,
  L_ab = L_ab,
  sigma0_hat = sigma0_hat,   # fixed, no extra noise
  
  # Group mean priors for log per-year rate
  mu_logLambda_prior = as.array(gp$mu_logLambda_prior),
  sd_logLambda_prior = as.array(gp$sd_logLambda_prior),
  
  # Age-split Beta prior (no hand-tuned kappa; strength from DHS)
  a0 = as.array(beta_ab_by_gid$a0),
  b0 = as.array(beta_ab_by_gid$b0),
  
  # Regularized horseshoe slab scale (computational stability)
  slab_scale = 1.0
)

# Consistency checks
stopifnot(nrow(stan_data$PHI) == stan_data$C)
stopifnot(ncol(stan_data$PHI) == stan_data$M_hsgp)
stopifnot(length(stan_data$lambda_hsgp) == stan_data$M_hsgp)
stopifnot(length(stan_data$a0) == stan_data$K && length(stan_data$b0) == stan_data$K)

cat("Stan data prepared successfully\n")

# ===============================================================================
# STEP 7: Stan model definition
# ===============================================================================

stancode <- "
data {
  int<lower=1> C;
  int<lower=1> K;
  int<lower=1> Kc;
  array[C] int<lower=1,upper=K>  group_id;
  array[C] int<lower=1,upper=Kc> country_id;
  int<lower=0> C_obs;
  array[C_obs] int<lower=1,upper=C> idx_obs;

  array[C] int<lower=0> y;
  vector[C] logE;
  array[C] int<lower=0,upper=1> is_obs;

  // HSGP
  int<lower=1> M_hsgp;
  matrix[C, M_hsgp] PHI;
  vector<lower=0>[M_hsgp] lambda_hsgp; // omega^2

  int<lower=1> n_groups;
  array[n_groups] int<lower=1> group_basis_start;
  array[n_groups] int<lower=1> group_basis_size;

  real<lower=0> slab_scale;

  // DHS neonatal bridge
  matrix[Kc,3] ab_mean;
  array[Kc] matrix[3,3] L_ab;
  vector<lower=0>[Kc] sigma0_hat;       // fixed residual SD

  // Group mean priors (per-year log-rate)
  array[K] real mu_logLambda_prior;
  array[K] real<lower=0> sd_logLambda_prior;

  // Age-split prior (DHS-only)
  array[K] real<lower=0> a0;
  array[K] real<lower=0> b0;
}
parameters {
  vector[K] mu_g;
  real<lower=0.1, upper=5.0> sigma;          // GP marginal SD
  real<lower=0.1, upper=10.0> rho;           // GP range
  vector[M_hsgp] beta_hsgp;                  // HSGP weights

  vector[C] z_v;
  vector<lower=0>[C] lambda;
  vector<lower=0>[K] tau_g;

  vector<lower=0,upper=1>[K] p_g;            // age split (sampled)
}
transformed parameters {
  vector[C] u;
  vector[C] v;
  vector[C] log_lambda_total;
  vector[C] log_lambda_1_11;
  vector[C] log_lambda_12_59;

  // Matérn 3/2 in 2D: spectral density S(w) ∝ kappa^3 (kappa^2 + w^2)^(-2.5)
  real<lower=0> kappa = sqrt(3.0) / rho;
  {
    vector[M_hsgp] sqrt_eigenvalues;
    vector[M_hsgp] beta_scaled;
    for (m in 1:M_hsgp) {
      real sqrt_spd = sigma
        * sqrt( 3.0 * pi()
                * pow(kappa, 3)
                * pow(kappa^2 + lambda_hsgp[m], -2.5) );
      sqrt_eigenvalues[m] = sqrt_spd;
      beta_scaled[m] = sqrt_eigenvalues[m] * beta_hsgp[m];
    }
    u = PHI * beta_scaled;
  }
  
  for (i in 1:C) {
    real tau_lam = tau_g[group_id[i]] * lambda[i];
    real c2 = square(slab_scale);
    real lt2 = square(slab_scale * lambda[i]) / (c2 + square(tau_lam));
    real lambda_tilde = sqrt(fmax(lt2, 1e-10));
    v[i] = tau_g[group_id[i]] * lambda_tilde * z_v[i];
    log_lambda_total[i] = mu_g[group_id[i]] + u[i] + v[i];
  }
  
  // Age-specific log rates via sampled p_g
  for (i in 1:C) {
    int g = group_id[i];
    real p = p_g[g];
    log_lambda_1_11[i]  = log_lambda_total[i] + log(p)   + log(59.0/11.0);
    log_lambda_12_59[i] = log_lambda_total[i] + log1m(p) + log(59.0/48.0);
  }
}
model {
  // group means (from DHS H1_59 priors)
  for (g in 1:K)
    mu_g[g] ~ normal(mu_logLambda_prior[g], sd_logLambda_prior[g]);

  // age split: *only* DHS-implied Beta prior
  p_g ~ beta(a0, b0);

  // GP + horseshoe
  sigma    ~ normal(0, 1);
  rho      ~ lognormal(0, 1);
  beta_hsgp ~ normal(0, 1);
  z_v      ~ normal(0, 1);
  lambda   ~ student_t(3, 0, 1);
  tau_g    ~ student_t(3, 0, 1);

  // likelihood (placebo clusters only)
  y[idx_obs] ~ poisson_log(logE[idx_obs] + log_lambda_total[idx_obs]);
}
generated quantities {
  array[Kc] vector[3] ab_rng;             // country-specific coef draws
  array[C] real H1_11_true;
  array[C] real H12_59_true;
  array[C] real H1_total_true;
  array[C] real H1_11_pred;
  array[C] real H12_59_pred;
  array[C] real H1_total_pred;
  array[C] real H0_pred;
  array[C] real imr_0_11;
  array[C] real u5mr_0_59;

  // Draw country-level neonatal coefficients; residual SD fixed to sigma0_hat
  for (k in 1:Kc) {
    vector[3] mu_k = to_vector(ab_mean[k]');
    ab_rng[k] = multi_normal_cholesky_rng(mu_k, L_ab[k]);
  }

  // Predictive cluster hazards and mortality
  real c2 = square(slab_scale);
  for (i in 1:C) {
    int g  = group_id[i];
    int kc = country_id[i];

    // latent hazards (fitted)
    real logH1_total_true = log_lambda_total[i] + log(59.0/12.0);
    real logH1_11_true    = log_lambda_1_11[i] + log(11.0/12.0);
    real logH12_59_true   = log_lambda_12_59[i] + log(48.0/12.0);
    H1_total_true[i] = exp(logH1_total_true);
    H1_11_true[i]    = exp(logH1_11_true);
    H12_59_true[i]   = exp(logH12_59_true);

    // predictive new horseshoe innovation
    real lambda_new = sqrt(square(student_t_rng(3, 0, 1)));
    real z_new = normal_rng(0, 1);
    real tau_lam_new = tau_g[g] * lambda_new;
    real lt2_new = square(slab_scale * lambda_new) / (c2 + square(tau_lam_new));
    real lambda_tilde_new = sqrt(lt2_new);
    real v_new = tau_g[g] * lambda_tilde_new * z_new;

    real log_lambda_total_pred = mu_g[g] + u[i] + v_new;
    real p = p_g[g];
    real log_lambda_1_11_pred  = log_lambda_total_pred + log(p)   + log(59.0/11.0);
    real log_lambda_12_59_pred = log_lambda_total_pred + log1m(p) + log(59.0/48.0);

    real logH1_total_pred = log_lambda_total_pred + log(59.0/12.0);
    real logH1_11_pred    = log_lambda_1_11_pred + log(11.0/12.0);
    real logH12_59_pred   = log_lambda_12_59_pred + log(48.0/12.0);
    H1_total_pred[i] = exp(logH1_total_pred);
    H1_11_pred[i]    = exp(logH1_11_pred);
    H12_59_pred[i]   = exp(logH12_59_pred);

    // Neonatal hazard via DHS bridge
    {
      vector[3] ab = ab_rng[kc];
      real mu_logH0 = ab[1]
                    + ab[2] * log(H1_11_pred[i])
                    + ab[3] * log(H12_59_pred[i]);
      real logH0    = normal_rng(mu_logH0, sigma0_hat[kc]);  // fixed residual SD
      H0_pred[i]    = exp(logH0);
    }

    // Mortality rates
    imr_0_11[i]  = 1 - exp(-(H0_pred[i] + H1_11_pred[i]));
    u5mr_0_59[i] = 1 - exp(-(H0_pred[i] + H1_11_pred[i] + H12_59_pred[i]));
  }
}
"

# ===============================================================================
# STEP 8: Model compilation and fitting
# ===============================================================================

cat("Compiling and fitting Stan model...\n")

# Initialization function
set.seed(42)
init_fun <- function() {
  # Initialize p_g at prior mean for stability
  p0 <- with(beta_ab_by_gid, a0 / (a0 + b0))
  list(
    mu_g    = rnorm(stan_data$K, stan_data$mu_logLambda_prior,
                    pmax(0.05 * stan_data$sd_logLambda_prior, 0.01)),
    sigma   = runif(1, 0.5, 1.5),
    rho     = runif(1, 1.0, 3.0),
    beta_hsgp = rnorm(stan_data$M_hsgp, 0, 0.1),
    z_v     = rnorm(stan_data$C, 0, 0.1),
    lambda  = pmax(abs(rt(stan_data$C, df = 3) * 0.1), 0.01),
    tau_g   = pmax(abs(rt(stan_data$K, df = 3) * 0.1), 0.01),
    p_g     = pmin(pmax(p0, 1e-3), 1-1e-3)
  )
}

# CmdStan setup (adjust path as needed)
if (Sys.info()["sysname"] == "Windows") {
  Sys.setenv(CMDSTANR_USE_WSL = "TRUE")
}

model_dir <- "../../models/stan"
stan_file <- cmdstanr::write_stan_file(stancode, dir = model_dir, basename = "hsgp_age_model")
mod <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)

# Fit the model
cat("Starting MCMC sampling...\n")
fit <- mod$sample(
  data = stan_data, init = init_fun,
  chains = 4, parallel_chains = 4, threads_per_chain = 2,
  iter_warmup = 1000, iter_sampling = 1500,
  adapt_delta = 0.90, max_treedepth = 12,
  refresh = 500, show_messages = TRUE, show_exceptions = TRUE
)

cat("Sampling completed!\n")

# ===============================================================================
# STEP 9: Posterior processing and diagnostics
# ===============================================================================

cat("Processing posterior draws...\n")

# Extract key quantities
dr <- fit$draws(variables = c(
  "H1_total_true","H1_11_true","H12_59_true",
  "H1_total_pred","H1_11_pred","H12_59_pred",
  "imr_0_11","u5mr_0_59"
))
rv <- posterior::as_draws_rvars(dr)

# Extract draws matrices
H1_total_pred_draws <- posterior::draws_of(rv$H1_total_pred)
H1_11_pred_draws    <- posterior::draws_of(rv$H1_11_pred)
H12_59_pred_draws   <- posterior::draws_of(rv$H12_59_pred)
H1_total_true_draws <- posterior::draws_of(rv$H1_total_true)
H1_11_true_draws    <- posterior::draws_of(rv$H1_11_true)
H12_59_true_draws   <- posterior::draws_of(rv$H12_59_true)
imr_draws           <- posterior::draws_of(rv$imr_0_11)
u5mr_draws          <- posterior::draws_of(rv$u5mr_0_59)

# Create posterior summary
post_summary <- tibble::tibble(
  row_id       = df.mod$row_id,
  community_id = df.mod$community_id,
  trial_phase  = df.mod$trial_phase,
  country      = df.mod$country,
  group_id     = df.mod$group_id,
  trt          = df.mod$trt,
  latitude     = df.mod$latitude,
  longitude    = df.mod$longitude,
  
  H1_total_pred_mean = apply(H1_total_pred_draws, 2, mean),
  H1_total_true_mean = apply(H1_total_true_draws,  2, mean),
  H1_total_pred_sd   = apply(H1_total_pred_draws,  2, sd),
  H1_total_true_sd   = apply(H1_total_true_draws,  2, sd),
  
  H1_11_pred_mean    = apply(H1_11_pred_draws, 2, mean),
  H1_11_true_mean    = apply(H1_11_true_draws, 2, mean),
  H1_11_pred_sd      = apply(H1_11_pred_draws, 2, sd),
  H1_11_true_sd      = apply(H1_11_true_draws, 2, sd),
  
  H12_59_pred_mean   = apply(H12_59_pred_draws, 2, mean),
  H12_59_true_mean   = apply(H12_59_true_draws, 2, mean),
  H12_59_pred_sd     = apply(H12_59_pred_draws, 2, sd),
  H12_59_true_sd     = apply(H12_59_true_draws, 2, sd),
  
  post_imr           = apply(imr_draws,  2, mean),
  post_imr_sd        = apply(imr_draws,  2, sd),
  post_u5mr          = apply(u5mr_draws, 2, mean),
  post_u5mr_sd       = apply(u5mr_draws, 2, sd),
  
  post_u5mr_log      = apply(log(pmax(u5mr_draws, 1e-12)), 2, mean),
  post_u5mr_log_sd   = apply(log(pmax(u5mr_draws, 1e-12)), 2, sd),
  post_imr_log       = apply(log(pmax(imr_draws,  1e-12)), 2, mean),
  post_imr_log_sd    = apply(log(pmax(imr_draws,  1e-12)), 2, sd)
)

# Merge with original data, handling duplicates
df.mod <- df.mod %>%
  dplyr::left_join(post_summary,
            by = dplyr::join_by(row_id, community_id, trial_phase, country, trt, group_id)) %>%
  dplyr::group_by(community_id) %>% 
  dplyr::arrange(desc(trial_phase), .by_group = TRUE) %>%
  dplyr::slice(1) %>% 
  dplyr::ungroup()

cat("Posterior processing completed\n")

# ===============================================================================
# STEP 10: Save results
# ===============================================================================

cat("Saving results...\n")

# Save key outputs
save(fit, file = "../../results/posterior_hsgp_fit.rda")
save(df.mod, gp, file = "../../results/baseline_mortality_estimates.rda")
save(beta_prior_by_group, file = "../../results/age_beta_priors_from_dhs.rda")

# Summary statistics
age_decomp_summary <- df.mod %>%
  dplyr::group_by(trial_phase, country, trt) %>%
  dplyr::summarise(
    n_clusters = n(),
    imr_mean   = mean(post_imr  * 1000),
    imr_sd     = sd(post_imr    * 1000),
    u5mr_mean  = mean(post_u5mr * 1000),
    u5mr_sd    = sd(post_u5mr   * 1000),
    ratio_imr_u5mr = imr_mean / u5mr_mean,
    .groups = "drop"
  ) %>%
  dplyr::mutate(trial_country = paste(trial_phase, country))

# Create results directories
results_dir <- "../../results"
# Export CSV results
readr::write_csv(age_decomp_summary, file.path(results_dir, "age_decomp_summary.csv"))

cat("\n=== BASELINE MORTALITY ESTIMATION SUMMARY ===\n")
print(age_decomp_summary)

cat("\nFiles saved:\n")
cat("  - ../../results/posterior_hsgp_fit.rda\n")
cat("  - ../../results/baseline_mortality_estimates.rda\n")
cat("  - ../../results/age_beta_priors_from_dhs.rda\n")

cat("\nBaseline mortality estimation completed successfully!\n")

#########################################################################################
### Some maps of observed hazards, predicte hazard, predicted imr and predicted u5mr

map.df <- df.mod %>%
  mutate(y_tilde = ifelse(trt == 0, y_159 + 0.5, NA),
         log_h1_obs = ifelse(trt == 0 & y_tilde > 0 & E_yrs > 0,
                             log((59/12) * (y_tilde / E_yrs)), NA),
         log_h1_post = log(H1_total_pred_mean),
         post_imr = 1e3*post_imr, post_u5mr = 1e3*post_u5mr,
         latitude = latitude.x, longitude = longitude.x) %>%
  select(trial_phase, country, community_id, latitude, longitude, 
         log_h1_obs, log_h1_post, post_imr, post_u5mr) |>
  left_join(data_crt_geospatial |> 
              select(longitude, latitude, ihme_u5m_2015) |>
              mutate(`IHME 2015` = 1e3*ihme_u5m_2015))

get.surfaces.plot <- function(trial, cnt){
  
  require(sp)         # Always first: S4 classes used by most spatial packages
  require(raster)     # Next: depends on sp
  require(sf)         # Simple features (can go after raster/sp)
  require(gstat)      # If you use it, after sp/raster
  
  select <- dplyr::select
  
  iso = countrycode::countrycode(cnt, "country.name", "iso3c")
  
  # Load administrative shapefile
  admin <- sf::st_read(paste0("../../data/shps/gadm41_", iso, "_shp/gadm41_", iso, "_2.shp"), quiet = TRUE) |> 
    st_transform(4326)
  
  data_all <- map.df |> filter(trial_phase == trial, country == cnt) 
  
  data_hazard_plot <- data_all |>
    rename("Observed logH1" = log_h1_obs, 
           "Smoothed logH1" = log_h1_post, 
           "Post IMR" = post_imr, 
           "Post U5MR" = post_u5mr) |>
    gather(method, est, -c(trial_phase, country, community_id, latitude, longitude))
  
  data_hazard_sf <- st_as_sf(data_hazard_plot, coords = c("longitude", "latitude"), crs = 4326) |>
    left_join(data_hazard_plot |> select(trial_phase, country, community_id, longitude, latitude) |> distinct())
  
  palette <- c("#0571b0", "#92c5de", "#ffffbf", "#f4a582", "#ca0020")
  
  plot_df <- data_hazard_sf %>% filter(method == "Post U5MR")
  coords <- st_coordinates(plot_df)
  
  # Build spatial points object
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = plot_df %>% st_drop_geometry(),
                                 proj4string = CRS("+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs"))
  
  # Define padded bounding box for interpolation
  bbox <- st_bbox(plot_df)
  xpad <- diff(c(bbox["xmin"], bbox["xmax"])) * 0.05
  ypad <- diff(c(bbox["ymin"], bbox["ymax"])) * 0.05
  xlim <- c(bbox["xmin"] - xpad, bbox["xmax"] + xpad)
  ylim <- c(bbox["ymin"] - ypad, bbox["ymax"] + ypad)
  
  plot.map <- function(meth){
    df <- data_hazard_sf |> filter(method == meth, !is.na(est))
    
    ggplot() +
      geom_sf(data = admin, fill = NA, color = "gray30", size = 0.15) +
      geom_point(data = df, aes(x = longitude, y = latitude, fill = est),
                 size = 2.2, shape = 21, color = "black", stroke = 0.2) +
      scale_fill_gradientn(colors = palette, #limits = limits_shared, 
                           name = meth) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme_void() + 
      theme(
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
      ) 
  }
  
  gg.obs <- plot.map("Observed logH1")
  gg.smooth <- plot.map("Smoothed logH1")
  gg.ihme <- plot.map("IHME 2015")
  gg.imr <- plot.map("Post IMR")
  gg.u5mr <- plot.map("Post U5MR") 
  
  thetrial = ifelse(trial %in% c("MORDOR I/II", "MORDOR I"), "MORDOR", trial)
  titl <- paste(thetrial, cnt)
  
  require(patchwork)
  all.gg <- patchwork::wrap_plots(gg.obs, gg.smooth, 
                                  gg.ihme, gg.u5mr, ncol = 2) +
    plot_annotation(
      title = titl,
      subtitle = "Observed & Smoothed log hazards, IHME & posterior U5MR per 1,000"
    )
  
  message("5. Printing the combined figure")
  
  figs_dir <- "../../results/figs"
  
  pdf(file.path(figs_dir, paste0("Map ", titl, ".pdf")), 
      width = 7, height = 5)
  print(all.gg)
  dev.off()
}

get.surfaces.plot("MORDOR I/II", "Niger")
get.surfaces.plot("AVENIR", "Niger")
get.surfaces.plot("MORDOR I", "Tanzania")
get.surfaces.plot("MORDOR I", "Malawi")
get.surfaces.plot("CHAT", "Burkina Faso")

#############################################################################
# Lets create the posterior density plots but lets also include the 
# INLA SPDE to compare with the HSGP + HS
#############################################################################

spde.surfaces <- list()

get.spde.surface <- function(trial, cnt) {
  message(paste("Getting INLA surfaces for", trial, cnt, "..."))
  
  # Observed: all clusters for prediction
  data_obs <- data_analysis %>%
    filter(trial_phase == trial, country == cnt)
  
  print("Creating INLA mesh for SPDE")
  
  coords_obs <- as.matrix(data_obs %>% select(longitude, latitude))
  
  # Fit: placebo-only clusters
  data_fit <- data_obs %>% filter(trt == 0) %>% mutate(log_offset = log(person_time))
  coords_fit <- as.matrix(data_fit %>% select(longitude, latitude))
  
  # Mesh and SPDE
  mesh <- inla.mesh.2d(loc = coords_fit, max.edge = c(0.1, 0.5), cutoff = 0.05)
  spde <- inla.spde2.pcmatern(mesh, prior.range = c(1.2, 0.5), prior.sigma = c(0.3, 0.5))
  A_fit <- inla.spde.make.A(mesh = mesh, loc = coords_fit)
  A_pred <- inla.spde.make.A(mesh = mesh, loc = coords_obs)
  s.index <- inla.spde.make.index("spatial", spde$n.spde)
  
  # Stack for fitting (placebo only)
  stack_fit <- inla.stack(
    data = list(y = data_fit$n),
    A = list(A_fit, 1),
    effects = list(spatial = s.index, intercept = rep(1, nrow(data_fit))),
    tag = "fit"
  )
  
  # Stack for prediction (all clusters)
  stack_pred <- inla.stack(
    data = list(y = NA),
    A = list(A_pred, 1),
    effects = list(spatial = s.index, intercept = rep(1, nrow(data_obs))),
    tag = "pred"
  )
  
  stack_full <- inla.stack(stack_fit, stack_pred)
  
  print("Fitting INLA model")
  
  # Fit model
  formula <- y ~ 0 + intercept + f(spatial, model = spde)
  
  fit <- inla(
    formula,
    family = "nbinomial",  # <-- change here
    data = inla.stack.data(stack_full),
    control.predictor = list(A = inla.stack.A(stack_full), compute = TRUE),
    E = c(data_fit$person_time, rep(1, nrow(data_obs))),
    control.family = list(  # <-- use this to estimate overdispersion
      list(hyper = list(theta = list(prior = "loggamma", param = c(1, 0.01))))
    )
  )
  
  # Extract predictions for all clusters
  idx_pred <- inla.stack.index(stack_full, "pred")$data
  fitted_vals <- fit$summary.fitted.values[idx_pred, ]
  
  print("Creating data frame")
  
 data_obs %>%
    mutate(
      log_spde = log(59/12) + fitted_vals$mean) %>%
    select(community_id, latitude, longitude, trial_phase, country, log_spde)
}

data_analysis <-   data_crt_geospatial |>
  filter(!is.na(trt), !is.na(latitude), !is.na(longitude), !is.na(person_time)) |> 
  mutate(person_time = person_time/365.25) |>
  group_by(community_id, latitude, longitude, trt, trial_phase, country) |> 
  summarize(n = sum(n), person_time = sum(person_time), .groups = "drop") |> 
  mutate(trial_phase = ifelse(country == "Niger" & trial_phase != "AVENIR", "MORDOR I/II", trial_phase))

# Get all the SPDE surfaces

spde.surfaces[[1]] <- get.spde.surface("MORDOR I/II", "Niger")
spde.surfaces[[2]] <- get.spde.surface("AVENIR", "Niger")
spde.surfaces[[3]] <- get.spde.surface("MORDOR I", "Tanzania")
spde.surfaces[[4]] <- get.spde.surface("MORDOR I", "Malawi")
spde.surfaces[[5]] <- get.spde.surface("CHAT", "Burkina Faso")

spde.surfaces.df <- data.table::rbindlist(spde.surfaces)

plot_df <- df.mod %>%
  left_join(spde.surfaces.df) |>
  mutate(trial_country = paste(trial_phase, country),
         y_tilde = ifelse(trt == 0, y_159 + 0.5, NA),
         log_h1_obs = ifelse(trt == 0 & y_tilde > 0 & E_yrs > 0,
                             log((59/12) * (y_tilde / E_yrs)), NA),
         log_h1_post = log(H1_total_pred_mean),
         log_h1_spde = log_spde, 
         u5_obs = ifelse(!is.na(log_h1_obs), 1e3 * (1 - exp(-exp(log_h1_obs))), NA),
         u5_post = 1e3 * (1 - exp(-exp(log_h1_post)))) %>%
  select(trial_country, group_id, log_h1_obs, log_h1_post, log_h1_spde)

dens_long <- plot_df %>%
  rename(Observed = log_h1_obs, Posterior = log_h1_post, SPDE = log_h1_spde) %>%
  tidyr::pivot_longer(c(Observed, Posterior, SPDE), names_to = "source", values_to = "log_h1") %>%
  filter(!is.na(log_h1))
dens_long_h <- dens_long %>% mutate(H1 = exp(log_h1))

prior_df <- gp %>%
  left_join(df.mod %>% distinct(group_id, trial_phase, country) %>%
              mutate(trial_country = paste(trial_phase, country)),
            by = "group_id") %>%
  distinct(trial_country, group_id, .keep_all = TRUE) %>%
  transmute(trial_country, group_id,
            mu_prior, mu_prior_l = mu_prior - 1.96*sd_prior, mu_prior_u = mu_prior + 1.96*sd_prior)
prior_df_h <- prior_df %>% mutate(mu_prior_H1 = exp(mu_prior),
                                  mu_prior_H1_l = exp(mu_prior_l),
                                  mu_prior_H1_u = exp(mu_prior_u))
prior_df_u5m <- prior_df_h %>%
  distinct(trial_country, mu_prior_H1, mu_prior_H1_l, mu_prior_H1_u) %>%
  mutate(x   = 1e3*(1 - exp(-mu_prior_H1)),
         x_l = 1e3*(1 - exp(-mu_prior_H1_l)),
         x_u = 1e3*(1 - exp(-mu_prior_H1_u)),
         label_ci = sprintf("%.1f (%.1f–%.1f)", x, x_l, x_u))

x_range <- dens_long_h %>%
  mutate(x = 1e3*(1 - exp(-H1))) %>%
  group_by(trial_country) %>%
  summarise(xmin = min(x, na.rm=TRUE), xmax = max(x, na.rm=TRUE), .groups="drop")

prior_df_u5m2 <- prior_df_u5m %>%
  left_join(x_range, by="trial_country") %>%
  mutate(dx = 0.06*(xmax - xmin),
         x_text = pmin(x + dx, xmax - 0.02*(xmax - xmin)),
         y_text = Inf)

library(ggplot2)

plot_theme <- function() {
  theme_bw() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.margin = margin(b = 5),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      plot.margin = margin(5, 10, 5, 10)
    )
}

# Plot 1: Observed vs Posterior vs Prior
plot1 <- tryCatch({
  # Filter data for first plot
  dens_plot1 <- dens_long_h %>% 
    filter(source %in% c("Observed", "Posterior"))
  
  ggplot(dens_plot1, aes(x = 1e3*(1 - exp(-H1)), fill = source)) +
    # Prior interval background
    geom_rect(data = prior_df_u5m2,
              aes(xmin=x_l, xmax=x_u, ymin=-Inf, ymax=Inf),
              inherit.aes = FALSE, fill = "grey70", alpha = 0.3) +
    
    # Density curves
    geom_density(alpha = 0.7, adjust = 1, size = 0.5) +
    
    # Prior mean line
    geom_vline(data = prior_df_u5m2, aes(xintercept = x),
               linetype = "dashed", linewidth = 0.8, color = "black") +
    
    # Prior interval marker at bottom
    geom_segment(data = prior_df_u5m2, aes(x=x_l, xend=x_u, y=0, yend=0),
                 inherit.aes = FALSE, linewidth = 2, lineend = "round", color = "black") +
    
    # Prior CI labels
    geom_text(data = prior_df_u5m2, 
              aes(x = 5, y = Inf, label = paste0("Prior: ", label_ci)),
              inherit.aes = FALSE, 
              vjust = 1.3, hjust = 0, 
              size = 3, color = "black") +
    
    # Colors for observed vs posterior
    scale_fill_manual(values = c("Observed" = "#E69F00", 
                                 "Posterior" = "#0072B2"),
                      labels = c("Observed data", 
                                 "Model posterior")) +
    
    facet_wrap(~ trial_country, scales = "free", ncol = 3) +
    
    labs(title = "Observed Data vs Model Posterior",
         subtitle = "Grey area and dashed line show prior distribution",
         x = NULL,  # Remove x-axis title from top plot
         y = "Density") +
    
    scale_x_continuous(limits = c(0, 250)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    
    plot_theme()
}, error = function(e) { 
  warning("Plot 1 failed: ", e$message); 
  NULL 
})

# Plot 2: Posterior vs SPDE vs Prior  
plot2 <- tryCatch({
  # Filter data for second plot
  dens_plot2 <- dens_long_h %>% 
    filter(source %in% c("Posterior", "SPDE"))
  
  ggplot(dens_plot2, aes(x = 1e3*(1 - exp(-H1)), fill = source)) +
    # Prior interval background
    geom_rect(data = prior_df_u5m2,
              aes(xmin=x_l, xmax=x_u, ymin=-Inf, ymax=Inf),
              inherit.aes = FALSE, fill = "grey70", alpha = 0.3) +
    
    # Density curves
    geom_density(alpha = 0.7, adjust = 1, size = 0.5) +
    
    # Prior mean line
    geom_vline(data = prior_df_u5m2, aes(xintercept = x),
               linetype = "dashed", linewidth = 0.8, color = "black") +
    
    # Prior interval marker at bottom
    geom_segment(data = prior_df_u5m2, aes(x=x_l, xend=x_u, y=0, yend=0),
                 inherit.aes = FALSE, linewidth = 2, lineend = "round", color = "black") +
    

    # Colors for posterior vs SPDE
    scale_fill_manual(values = c("Posterior" = "#0072B2",
                                 "SPDE" = "#009E73"),
                      labels = c("Model posterior", 
                                 "SPDE comparison")) +
    
    facet_wrap(~ trial_country, scales = "free", ncol = 3) +
    
    labs(title = "Model Posterior vs SPDE Method",
         subtitle = "Grey area and dashed line show prior distribution",
         x = "Deaths per 1,000 surviving neonates", 
         y = "Density") +
    
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    
    plot_theme()
}, error = function(e) { 
  warning("Plot 2 failed: ", e$message); 
  NULL 
})


figs_dir <- "../../results/figs"

pdf(file.path(figs_dir, "prior_post_u5m.pdf"), width = 7, height = 5)
if (!is.null(plot1)) print(plot1)
dev.off()

pdf(file.path(figs_dir, "spde_post_u5m.pdf"), width = 7, height = 5)
if (!is.null(plot2)) print(plot2)
dev.off()


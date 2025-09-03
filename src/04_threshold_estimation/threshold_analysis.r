# ===============================================================================
# REACH Secondary Analysis: Threshold Analysis
# ===============================================================================
# Purpose: Derive age- and context-specific cessation thresholds
# Author: William Msemburi
# Date: August 2025
# ===============================================================================

library(tidyverse)
library(posterior)
library(RColorBrewer)
library(patchwork)
library(data.table)
library(cmdstanr)
library(bayesplot)
library(glue)

rm(list = ls())

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

start_time <- Sys.time()

# Load utility functions
source("../utils/reach_analysis_utils.R")

# ===============================================================================
# DATA PREPARATION
# ===============================================================================

cat("Loading data for threshold analysis...\n")

# Load trial data and baseline mortality estimates
if (!file.exists("../../data/Clean/data_crt_geospatial.rda")) {
  stop("Clean trial data not found")
}
load("../../data/Clean/data_crt_geospatial.rda")

if (!file.exists("../../results/baseline_mortality_estimates.rda")) {
  stop("Baseline mortality estimates not found. Run bayesian_hsgp_model.R first.")
}
load("../../results/baseline_mortality_estimates.rda")

cat("Data loaded successfully\n")

# Filter to single age group (1-59 months) for threshold analysis
data_analysis <- data_crt_geospatial %>%
  dplyr::filter(!is.na(trt), !trial_phase %in% "MORDOR III",
                !is.na(latitude), !is.na(longitude), !is.na(person_time), person_time > 0,
                age_start >= 1 & age_start <= 59)

# Aggregate and join with mortality estimates
df <- data_analysis %>%
  sum_deaths_and_person_time(c("country", "trial_phase", "community_id", "trt")) %>%
  dplyr::left_join(df.mod %>% 
                     dplyr::select(country, community_id, 
                                   post_u5mr, post_u5mr_log, post_u5mr_log_sd,
                                   post_imr, post_imr_log, post_imr_log_sd) %>% 
                     dplyr::distinct(),
                   by = dplyr::join_by(country, community_id)) %>%
  dplyr::filter(!is.na(post_u5mr_log), !is.na(post_u5mr_log_sd), 
                !is.na(post_imr_log), !is.na(post_imr_log_sd),
                deaths >= 0, person_time > 0)

# Create group indices and prepare data
group_index <- create_group_index(df)
df <- df %>%
  dplyr::left_join(group_index, by = c("trial_phase", "country")) %>%
  dplyr::arrange(group_id, community_id) %>%
  dplyr::filter(!is.na(group_id)) 

cat("Data prepared: ", nrow(df), " clusters across ", 
    dplyr::n_distinct(df$group_id), " trial-countries\n")

# ===============================================================================
# PLOTTING DATA PREPARATION
# ===============================================================================

# Prepare data for plotting
df_plot <- df %>% 
  dplyr::filter(deaths > 0) %>%
  dplyr::mutate(
    post_u5mr = round(post_u5mr, digits = 3),
    post_imr = round(post_imr, digits = 3)
  ) %>%
  dplyr::group_by(country, trial_phase, trt, post_u5mr, post_imr) %>%
  dplyr::summarise(deaths = sum(deaths), person_time = sum(person_time), 
                   .groups = "drop") %>%
  dplyr::mutate(
    post_u5mr_1000 = 1e3 * post_u5mr,
    post_imr_1000 = 1e3 * post_imr,
    person_time_years = person_time / 365.25,
    mortality_rate = log(deaths / person_time_years),
    trial_country = paste(trial_phase, country),
    trt = factor(trt, levels = c(0, 1), labels = c("Placebo", "Azithromycin"))
  ) %>%
  dplyr::filter(!is.na(post_u5mr_1000), !is.na(post_imr_1000), !is.na(mortality_rate))

u5_min <- 0.035
u5_max <- 0.235
im_min <- 0.020
im_max <- 0.110

cat("Mortality ranges - U5MR: [", round(u5_min*1000, 1), ", ", round(u5_max*1000, 1), 
    "], IMR: [", round(im_min*1000, 1), ", ", round(im_max*1000, 1), "]\n")

# ===============================================================================
# STAN MODEL FOR THRESHOLD ANALYSIS
# ===============================================================================

stancode_variable_model <- "
functions {
  real crossing_from_preds(int Np, array[,] real pred, vector log_x) {
    real neg_inf = negative_infinity();
    for (u in 1:(Np-1)) {
      real d0 = log(pred[u,2])   - log(pred[u,1]);
      real d1 = log(pred[u+1,2]) - log(pred[u+1,1]);
      if (abs(d0) < 1e-12) return exp(log_x[u]);
      if ((d0 < 0 && d1 > 0) || (d0 > 0 && d1 < 0)) {
        real t  = d0 / (d0 - d1);
        real lx = log_x[u] + t * (log_x[u+1] - log_x[u]);
        return exp(lx);
      }
    }
    return neg_inf;
  }
}

data {
  int<lower=1> N, J, N_pred;
  array[N] int<lower=0> deaths;
  vector<lower=0>[N] person_time;
  array[N] int<lower=0,upper=1> trt;
  array[N] int<lower=1,upper=J> group;
  vector[N] log_mortality_mean;
  vector<lower=0>[N] log_mortality_sd;
  vector[N_pred] log_mortality_pred;
}

parameters {
  vector[N] log_mortality_raw;
  
  // Hierarchical group effects
  vector[J] alpha_group_raw;
  real mu_alpha_group;
  real<lower=0> sigma_group;
  
  // Model coefficients
  real beta_mortality;
  real beta_trt;
  real beta_trt_mortality;
}

transformed parameters {
  vector[N] log_mortality = log_mortality_mean + log_mortality_sd .* log_mortality_raw;
  vector[J] alpha_group = mu_alpha_group + sigma_group * alpha_group_raw;
}

model {
  // Priors
  log_mortality_raw ~ std_normal();
  alpha_group_raw ~ std_normal();
  mu_alpha_group ~ normal(0, 1);
  sigma_group ~ normal(0, 1);
  
  beta_mortality ~ normal(0, 1);
  beta_trt ~ normal(0, 0.5);
  beta_trt_mortality ~ normal(0, 0.5);
  
  // Likelihood
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = alpha_group[group[i]] + 
            beta_mortality * log_mortality[i] +
            trt[i] * (beta_trt + beta_trt_mortality * log_mortality[i]) +
            log(person_time[i]);
  }
  deaths ~ poisson_log(mu);
}

generated quantities {
  // Predictions
  array[N_pred, 2] real pred_lambda;
  
  // Crossing-based threshold
  real threshold;
  int has_cross;
  
  real mean_alpha = mean(alpha_group);
  
  // Generate predictions
  for (u in 1:N_pred) {
    real mu_base = mean_alpha +
                   beta_mortality * log_mortality_pred[u] +
                   log(365.25);
    real mu_trt = mu_base +
                  beta_trt +
                  beta_trt_mortality * log_mortality_pred[u];
    
    pred_lambda[u, 1] = exp(mu_base);  // placebo
    pred_lambda[u, 2] = exp(mu_trt);   // treatment
  }
  
  // Calculate crossing threshold
  threshold = crossing_from_preds(N_pred, pred_lambda, log_mortality_pred);
  has_cross = (threshold != negative_infinity());
}
"

stancode_fixed_model <- "
functions {
  real crossing_from_preds(int Np, array[,] real pred, vector log_x) {
    real neg_inf = negative_infinity();
    for (u in 1:(Np-1)) {
      real d0 = log(pred[u,2])   - log(pred[u,1]);
      real d1 = log(pred[u+1,2]) - log(pred[u+1,1]);
      if (abs(d0) < 1e-12) return exp(log_x[u]);
      if ((d0 < 0 && d1 > 0) || (d0 > 0 && d1 < 0)) {
        real t  = d0 / (d0 - d1);
        real lx = log_x[u] + t * (log_x[u+1] - log_x[u]);
        return exp(lx);
      }
    }
    return neg_inf;
  }
}
data {
  int<lower=1> N, J, N_pred;
  array[N] int<lower=0> deaths;
  vector<lower=0>[N] person_time;
  array[N] int<lower=0,upper=1> trt;
  array[N] int<lower=1,upper=J> group;
  vector[N] log_mortality_mean;
  vector<lower=0>[N] log_mortality_sd;
  vector[N_pred] log_mortality_pred;
}
parameters {
  // Hierarchical group effects
  vector[J] alpha_group_raw;
  real mu_alpha_group;
  real<lower=0> sigma_group;
  
  // Model coefficients
  real beta_mortality;
  real beta_trt;
  real beta_trt_mortality;
}
transformed parameters {
  vector[J] alpha_group = mu_alpha_group + sigma_group * alpha_group_raw;
}
model {
  // Priors
  alpha_group_raw ~ std_normal();
  mu_alpha_group ~ normal(0, 1);
  sigma_group ~ normal(0, 1);
  
  beta_mortality ~ normal(0, 1);
  beta_trt ~ normal(0, 0.5);
  beta_trt_mortality ~ normal(0, 0.5);
  
  // Likelihood
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = alpha_group[group[i]] + 
            beta_mortality * log_mortality_mean[i] +
            trt[i] * (beta_trt + beta_trt_mortality * log_mortality_mean[i]) +
            log(person_time[i]);
  }
  deaths ~ poisson_log(mu);
}
generated quantities {
  // Predictions
  array[N_pred, 2] real pred_lambda;
  
  // Crossing-based threshold
  real threshold;
  int has_cross;
  
  real mean_alpha = mean(alpha_group);
  
  // Generate predictions
  for (u in 1:N_pred) {
    real mu_base = mean_alpha +
                   beta_mortality * log_mortality_pred[u] +
                   log(365.25);
    real mu_trt = mu_base +
                  beta_trt +
                  beta_trt_mortality * log_mortality_pred[u];
    
    pred_lambda[u, 1] = exp(mu_base);  // placebo
    pred_lambda[u, 2] = exp(mu_trt);   // treatment
  }
  
  // Calculate crossing threshold
  threshold = crossing_from_preds(N_pred, pred_lambda, log_mortality_pred);
  has_cross = (threshold != negative_infinity());
}
"

# ===============================================================================
# MODEL FITTING
# ===============================================================================

stan_data_u5mr <- make_stan_data(
  df, "post_u5mr_log", "post_u5mr_log_sd",
  u5_min, u5_max, n_pred = 100L
)

stan_data_imr <- make_stan_data(
  df, "post_imr_log", "post_imr_log_sd",
  im_min, im_max, n_pred = 100L
)

# Fit U5MR model
init_values_u5mr <- compute_smart_inits_single(df, "post_u5mr_log", n_chains = 4)
init_values_imr <- compute_smart_inits_single(df, "post_imr_log", n_chains = 4)

# Create model directory if it doesn't exist
model_dir <- "../../models/stan"
if (!dir.exists(model_dir)) {
  dir.create(model_dir, recursive = TRUE)
}

cat("Compiling Variable Stan model...\n")

# Write and compile Stan model
stan_file <- cmdstanr::write_stan_file(stancode_variable_model, dir = model_dir, 
                                       basename = "reach_mortality_model")
mod <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

cat("Fitting U5MR threshold model (variable mortality)...\n")

fit_u5mr <- mod$sample(
  data = stan_data_u5mr,
  init = init_values_u5mr,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99,
  max_treedepth = 15,
  refresh = 1000
)

cat("Fitting IMR threshold model (variable mortality)...\n")

# Fit IMR model

fit_imr <- mod$sample(
  data = stan_data_imr,
  init = init_values_imr,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99,
  max_treedepth = 15,
  refresh = 1000
)

cat("Model fitting completed\n")

cat("Compiling Fixed Stan models...\n")

stan_file_fixed <- cmdstanr::write_stan_file(stancode_fixed_model, dir = model_dir, 
                                             basename = "reach_mortality_model_fixed")
mod_fixed <- cmdstanr::cmdstan_model(stan_file_fixed, cpp_options = list(stan_threads = TRUE))

cat("Fitting U5MR threshold model (fixed mortality)...\n")

fit_u5mr_fixed <- mod_fixed$sample(
  data = stan_data_u5mr,
  init = init_values_u5mr,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99,
  max_treedepth = 15,
  refresh = 1000
)

cat("Fitting IMR threshold model (fixed mortality)...\n")

fit_imr_fixed <- mod_fixed$sample(
  data = stan_data_imr,
  init = init_values_imr,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99,
  max_treedepth = 15,
  refresh = 1000
)

cat("Model fitting completed\n")

# force-read so the objects don’t depend on CSVs later
invisible(fit_u5mr$draws()); 
invisible(fit_imr$draws());
invisible(fit_u5mr_fixed$draws());
invisible(fit_imr_fixed$draws());

# ===============================================================================
# THRESHOLD EXTRACTION
# ===============================================================================

cat("Extracting thresholds from Stan models...\n")

# Function to extract finite thresholds
extract_finite_thresholds <- function(fit) {
  dd <- fit$draws(variables = c("threshold", "has_cross"), format = "draws_df") %>%
    as_draws_df()
  # Convert -Inf to NA to play nice with summaries
  dd$threshold[is.infinite(dd$threshold)] <- NA_real_
  dd %>% filter(has_cross == 1)
}

# Extract thresholds from both models
thresholds_u5mr <- extract_finite_thresholds(fit_u5mr)
thresholds_imr <- extract_finite_thresholds(fit_imr)

thresholds_u5mr_fixed <- extract_finite_thresholds(fit_u5mr_fixed)
thresholds_imr_fixed <- extract_finite_thresholds(fit_imr_fixed)

# Summary function for quantiles
qs95 <- function(df){
  df %>%
    summarise(
      n = n(),
      threshold_mean_1000 = 1e3*mean(threshold, na.rm = TRUE),
      threshold_med_1000 = 1e3*median(threshold, na.rm = TRUE),
      sd = sd(threshold, na.rm = TRUE),
      threshold_lo_1000 = 1e3*quantile(threshold, 0.025, na.rm = TRUE),
      threshold_hi_1000 = 1e3*quantile(threshold, 0.975, na.rm = TRUE)
    )
}

# Create threshold results summary
threshold_results <- rbind(data.frame(qs95(thresholds_u5mr)),
                           data.frame(qs95(thresholds_imr))) |>
  mutate(model = c("U5MR", "IMR"))

threshold_results_fixed <- rbind(data.frame(qs95(thresholds_u5mr_fixed)),
                                 data.frame(qs95(thresholds_imr_fixed))) |>
  mutate(model = c("U5MR", "IMR"))  
# ===============================================================================
# COEFFICIENT EXTRACTION
# ===============================================================================

# Extract coefficients from both models
coefficients_u5mr <- fit_u5mr$draws(variables = c("beta_mortality", "beta_trt", "beta_trt_mortality"), 
                                    format = "draws_df") %>%
  posterior::summarise_draws(mean, mad, ~quantile(.x, probs = c(0.025, 0.975))) %>%
  dplyr::mutate(model = "U5MR")

coefficients_imr <- fit_imr$draws(variables = c("beta_mortality", "beta_trt", "beta_trt_mortality"), 
                                  format = "draws_df") %>%
  posterior::summarise_draws(mean, mad, ~quantile(.x, probs = c(0.025, 0.975))) %>%
  dplyr::mutate(model = "IMR")

all_coefficients <- dplyr::bind_rows(coefficients_u5mr, coefficients_imr)

# ===============================================================================
# PREDICTION EXTRACTION
# ===============================================================================

# Function to extract predictions
extract_predictions_single <- function(fit, stan_data) {
  # Extract prediction draws
  pred_draws <- fit$draws(variables = "pred_lambda", format = "draws_matrix")
  
  # Get array indices from parameter names
  param_names <- colnames(pred_draws)
  index_pattern <- "pred_lambda\\[(\\d+),(\\d+)\\]"
  index_mat <- do.call(rbind, stringr::str_match_all(param_names, index_pattern))
  index_df <- index_mat[,2:3] |> apply(2, as.integer) |> as.data.frame()
  names(index_df) <- c("i", "t")
  
  mortality_pred_seq <- exp(stan_data$log_mortality_pred)
  
  # Create prediction grid
  pred_grid <- index_df |>
    dplyr::mutate(
      mortality = mortality_pred_seq[i],
      trt = factor(t - 1, levels = 0:1, labels = c("Placebo", "Azithromycin")),
      mortality_1000 = mortality * 1000
    )
  
  # Summarize posterior predictions
  lambda_summary <- apply(pred_draws, 2, function(x) {
    c(mean = mean(x), lwr = quantile(x, 0.025), upr = quantile(x, 0.975))
  }) |> t() |> as.data.frame()
  
  # Merge predictions with grid
  lambda_preds_df <- dplyr::bind_cols(pred_grid, lambda_summary) |>
    dplyr::mutate(
      pred_mean = log(mean),
      pred_lb   = log(`lwr.2.5%`),
      pred_ub   = log(`upr.97.5%`)
    )
  
  return(lambda_preds_df)
}

# Extract predictions for both models
lambda_preds_u5mr <- extract_predictions_single(fit_u5mr, stan_data_u5mr)
lambda_preds_imr <- extract_predictions_single(fit_imr, stan_data_imr)

# ===============================================================================
# VISUALIZATION
# ===============================================================================

cat("Creating visualizations...\n")

library(ggplot2)
library(ggokabeito)
library(showtext)

# One-time per session: fetch & register the font with showtext/sysfonts
font_add_google("News Cycle", "News Cycle")
showtext_auto()  # turn on showtext for all devices

theme_clean <- function() {
  theme_minimal(base_family = "News Cycle") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(face = "bold"))
}

# Function to create mortality plots
create_mortality_plot <- function(df_plot, lambda_preds_df, threshold_results_row, model_type) {
  
  # Extract threshold values
  threshold_mean <- threshold_results_row$threshold_mean_1000 / 1000
  threshold_lo <- threshold_results_row$threshold_lo_1000 / 1000
  threshold_hi <- threshold_results_row$threshold_hi_1000 / 1000
  
  # Prepare observed data
  if (model_type == "u5mr") {
    mortality_col <- "post_u5mr_1000"
    mortality_label <- "U5MR (per 1,000 live births)"
    plot_title <- "Observed and Predicted Mortality vs. U5MR"
  } else {
    mortality_col <- "post_imr_1000"  
    mortality_label <- "IMR (per 1,000 live births)"
    plot_title <- "Observed and Predicted Mortality vs. IMR"
  }
  
  # Get values for annotations
  ovals <- lambda_preds_df |> 
    dplyr::filter(mortality == max(mortality)) |> 
    dplyr::pull(pred_mean)
  
  # Trial labels
  label_df <- tibble::tibble(
    trial_country = c(
      "MORDOR I Tanzania", "MORDOR I Malawi", "CHAT Burkina Faso",
      "AVENIR Niger", "MORDOR I/II Niger", "Threshold", "treat", "control"
    ),
    label = c("TZA MORDOR", "MWI MORDOR", "BFA CHAT", "NER AVENIR", 
              "NER MORDOR", 
              paste("Baseline", toupper(model_type), "Threshold"), 
              "Control", "Treated"),
    x = c(
      if (model_type == "u5mr") c(22, 57, 63, 155, 130, 1e3*threshold_mean + 2, 237, 237)
      else c(12, 35, 42, 88, 78, 1e3*threshold_mean + 2, 112, 112)
    ),
    y = 1e3*exp(c(-4, -5.3, -5.7, -4.0, -3, -6.5, max(ovals), min(ovals))),
    label_color = c(
      RColorBrewer::brewer.pal(5, "Set1")[c(4, 3, 2, 1, 5)],
      "navy", "grey40", "firebrick"
    ),
    hjust = 0
  )
  
  # Create plot
  p <- ggplot2::ggplot(df_plot, 
                       ggplot2::aes(x = .data[[mortality_col]], y = 1e3*exp(mortality_rate))) +
    
    # Observed points
    ggplot2::geom_point(
      ggplot2::aes(
        fill = trial_country,
        alpha = trt,
        shape = trt,
        size = person_time_years
      ),
      stroke = 0.2,
      color = "black"
    ) +
    
    # Prediction ribbons
    ggplot2::geom_ribbon(
      data = lambda_preds_df |> dplyr::filter(trt == "Placebo"),
      ggplot2::aes(x = mortality_1000, ymin = 1e3*exp(pred_lb), ymax = 1e3*exp(pred_ub)),
      fill = "grey40", alpha = 0.2, inherit.aes = FALSE
    ) +
    ggplot2::geom_ribbon(
      data = lambda_preds_df |> dplyr::filter(trt == "Azithromycin"),
      ggplot2::aes(x = mortality_1000, ymin = 1e3*exp(pred_lb), ymax = 1e3*exp(pred_ub)),
      fill = "firebrick", alpha = 0.2, inherit.aes = FALSE
    ) +
    
    # Predicted lines
    ggplot2::geom_line(
      data = lambda_preds_df |> dplyr::filter(trt == "Placebo"),
      ggplot2::aes(x = mortality_1000, y = 1e3*exp(pred_mean)),
      color = "grey40", linewidth = 0.6
    ) +
    ggplot2::geom_line(
      data = lambda_preds_df |> dplyr::filter(trt == "Azithromycin"),
      ggplot2::aes(x = mortality_1000, y = 1e3*exp(pred_mean)),
      color = "firebrick", linewidth = 0.6
    ) +
    
    # Threshold lines
    ggplot2::geom_vline(xintercept = 1000 * threshold_mean, lty = 2, color = "navy", linewidth = 1.0) +
    ggplot2::geom_vline(xintercept = 1000 * threshold_lo, lty = 3, color = "navy", linewidth = 0.5, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 1000 * threshold_hi, lty = 3, color = "navy", linewidth = 0.5, alpha = 0.8) +
    
    # Trial and line annotations
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = x, y = y, label = label, color = label_color),
      inherit.aes = FALSE,
      size = 3,
      hjust = label_df$hjust,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    
    # Scales and aesthetics
    ggplot2::scale_alpha_manual(values = c("Placebo" = 0.2, "Azithromycin" = 0.6)) +
    ggplot2::scale_shape_manual(values = c("Placebo" = 22, "Azithromycin" = 21)) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    
    # Labels
    ggplot2::labs(
      title = plot_title,
      x = mortality_label,
      y = "Deaths per 1,000 Person-Yr"
    ) +
    
    # Theme
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_y_log10() +
    
    # Adjust x-axis limits based on model type
    ggplot2::xlim(if (model_type == "u5mr") c(10, 240) else c(10, 120))
  
  # Add embedded legend
  legend_x_min <- if (model_type == "u5mr") 160 else 85
  legend_x_max <- if (model_type == "u5mr") 213 else 118
  
  p <- p + 
    ggplot2::annotate("rect", xmin = legend_x_min, xmax = legend_x_max, 
                      ymin = 1e3*exp(-6), ymax = 1e3*exp(-5.15),
                      fill = "white", color = "grey60", linewidth = 0.3) +
    ggplot2::annotate("point", x = legend_x_min + 5, y = 1e3*exp(-5.4), 
                      shape = 22, size = 3.5, fill = "grey60", alpha = 0.4) +
    ggplot2::annotate("point", x = legend_x_min + 5, y = 1e3*exp(-5.65), 
                      shape = 21, size = 3.5, fill = "grey60", alpha = 1) +
    ggplot2::annotate("text", x = legend_x_min + 10, y = 1e3*exp(-5.4), 
                      label = "Observed placebo", hjust = 0, size = 3) +
    ggplot2::annotate("text", x = legend_x_min + 10, y = 1e3*exp(-5.65), 
                      label = "Observed azithromycin", hjust = 0, size = 3)
  
  return(p)
}

# Create plots
plot_u5mr <- create_mortality_plot(df_plot, lambda_preds_u5mr, 
                                   threshold_results[threshold_results$model == "U5MR", ], "u5mr")

plot_imr <- create_mortality_plot(df_plot, lambda_preds_imr, 
                                  threshold_results[threshold_results$model == "IMR", ], "imr")

# Summary plots
threshold_plot <- threshold_results %>%
  ggplot2::ggplot(ggplot2::aes(x = model, y = threshold_mean_1000, color = model)) +
  ggplot2::geom_point(size = 4) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = threshold_lo_1000, ymax = threshold_hi_1000), 
                         width = 0.2, linewidth = 1.2) +
  ggplot2::labs(
    title = "MDA Stopping Thresholds",
    subtitle = "U5MR and IMR models with crossing-based threshold calculation",
    x = "Baseline Mortality Measure",
    y = "Threshold (per 1,000 live births)",
    caption = "Points: mean thresholds; Error bars: 95% credible intervals", 
    color = "Model" 
  ) +
  ggplot2::scale_color_manual(values = c("IMR" = "darkred", "U5MR" = "darkblue")) +
  ggplot2::theme_minimal(base_size = 12) +   scale_fill_okabe_ito() +
  scale_color_okabe_ito() +  theme_clean() + 
  ggplot2::theme(legend.position = "none")  

coef_plot <- all_coefficients %>%
  dplyr::mutate(
    coefficient_name = dplyr::case_when(
      stringr::str_detect(variable, "beta_mortality") ~ "Mortality Effect",
      stringr::str_detect(variable, "beta_trt$") ~ "Treatment Effect",
      stringr::str_detect(variable, "beta_trt_mortality") ~ "Treatment × Mortality"
    )
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = coefficient_name, y = mean, color = model)) +
  ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(width = 0.3)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`), 
                         width = 0.2, position = ggplot2::position_dodge(width = 0.3)) +
  ggplot2::labs(
    title = "Coefficient Comparison",
    x = "Coefficient Type",
    y = "Coefficient Value",
    color = "Model"
  ) +
  ggplot2::scale_color_manual(values = c("IMR" = "darkred", "U5MR" = "darkblue")) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +   scale_fill_okabe_ito() +
  scale_color_okabe_ito() +  theme_clean() + theme(legend.position = "top") 

# Combined visualizations
analysis_plots <- threshold_plot | coef_plot 

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

cat("Saving results...\n")

# Create results directories
results_dir <- "../../results"
figs_dir <- "../../results/figs"

# Save model fits and data
save(fit_u5mr, fit_imr, 
     fit_u5mr_fixed, fit_imr_fixed,  
     stan_data_u5mr, stan_data_imr, df, df_plot, group_index,
     file = file.path(results_dir, "threshold_models.rda"))

# Save threshold results
save(threshold_results, all_coefficients, 
     file = file.path(results_dir, "threshold_analysis_results.rda"))

# Save plots

pdf(paste0(figs_dir, "/threshold_analysis_plots.pdf"), height = 5, width = 8)
print(analysis_plots)
dev.off()

pdf(paste0(figs_dir, "/threshold_u5mr.pdf"), height = 5, width = 8)
print(plot_u5mr)
dev.off()

pdf(paste0(figs_dir, "/threshold_imr.pdf"), height = 5, width = 8)
print(plot_imr)
dev.off()

# Export CSV results
readr::write_csv(threshold_results, file.path(results_dir, "threshold_results.csv"))
readr::write_csv(all_coefficients, file.path(results_dir, "threshold_coefficients.csv"))

# ===============================================================================
# MODEL DIAGNOSTICS AND SUMMARY
# ===============================================================================

cat("\n=== MODEL DIAGNOSTICS ===\n")
cat("U5MR Model:\n")
u5mr_diagnostics <- fit_u5mr$diagnostic_summary()
print(u5mr_diagnostics)

cat("\nIMR Model:\n")
imr_diagnostics <- fit_imr$diagnostic_summary()
print(imr_diagnostics)

if (any(c(u5mr_diagnostics$num_divergent, imr_diagnostics$num_divergent) > 0)) {
  cat("Warning: One or more models have divergent transitions. Consider increasing adapt_delta.\n")
}

# Print results summary
cat("\n=== THRESHOLD ANALYSIS RESULTS ===\n")
print(threshold_results)

cat("\n=== COEFFICIENT ESTIMATES ===\n")
coef_display <- all_coefficients %>% 
  dplyr::select(model, variable, mean, `2.5%`, `97.5%`) %>%
  dplyr::arrange(variable, model)
print(coef_display)

cat("\n=== SUMMARY ===\n")
cat("U5MR threshold:", round(threshold_results$threshold_mean_1000[threshold_results$model == "U5MR"], 1), 
    "per 1,000 live births (95% CI:", round(threshold_results$threshold_lo_1000[threshold_results$model == "U5MR"], 1),
    "-", round(threshold_results$threshold_hi_1000[threshold_results$model == "U5MR"], 1), ")\n")
cat("IMR threshold: ", round(threshold_results$threshold_mean_1000[threshold_results$model == "IMR"], 1), 
    "per 1,000 live births (95% CI:", round(threshold_results$threshold_lo_1000[threshold_results$model == "IMR"], 1),
    "-", round(threshold_results$threshold_hi_1000[threshold_results$model == "IMR"], 1), ")\n\n")

cat("\nFiles saved:\n")
cat("  - ../../results/threshold_models.rda\n")
cat("  - ../../results/threshold_analysis_results.rda\n")
cat("  - ../../results/threshold_results.csv\n")
cat("  - ../../results/threshold_coefficients.csv\n")
cat("  - ../../results/figs/threshold_analysis_plots.png\n")
cat("  - ../../results/figs/threshold_mortality_plots.png\n")

end_time <- Sys.time()
time_elapsed <- round(as.numeric(end_time - start_time, units = 'mins'), 1)

cat("\nThreshold analysis completed in", time_elapsed, "minutes\n")

message("Threshold analysis completed successfully!")

reg_dens <- function(){
  require(dplyr)
  require(tidyr)
  require(bayesplot)
  require(patchwork)
  
  df <- bind_rows(
    posterior::as_draws_df(fit_u5mr$draws("threshold")) |>
      transmute(value = 1000 * threshold, model = "U5MR"),
    posterior::as_draws_df(fit_imr$draws("threshold")) |>
      transmute(value = 1000 * threshold, model = "IMR")
  ) |>
    filter(is.finite(value))
  
  wide <- df %>%
    group_by(model) %>% mutate(.row = row_number()) %>% ungroup() %>%
    pivot_wider(id_cols = .row, names_from = model, values_from = value) %>%
    select(-.row) %>%
    drop_na()                                # ensure rectangular
  
  mat <- as.matrix(wide)                      # iterations × parameters
  colnames(mat) <- colnames(wide)             # e.g., "IMR", "U5MR"
  
  thresh <- bayesplot::mcmc_areas(mat, prob = 0.95) +
    labs(x = "Threshold (per 1,000 live births)", y = NULL, title = "Crossing thresholds")
  
  df <- bind_rows(
    posterior::as_draws_df(fit_u5mr$draws("beta_trt")) |>
      transmute(value = beta_trt, model = "U5MR", coeff = "Treatment"),
    posterior::as_draws_df(fit_imr$draws("beta_trt")) |>
      transmute(value = beta_trt, model = "IMR", coeff = "Treatment"),
    posterior::as_draws_df(fit_u5mr$draws("beta_trt_mortality")) |>
      transmute(value = beta_trt_mortality, model = "U5MR", coeff = "Treatment x Mortality"), 
    posterior::as_draws_df(fit_imr$draws("beta_trt_mortality")) |>
      transmute(value = beta_trt_mortality, model = "IMR", coeff = "Treatment x Mortality")
  ) |> filter(is.finite(value)) |>
    mutate(model = paste(model, coeff), coeff = NULL)
  
  wide <- df %>%
    group_by(model) %>% mutate(.row = row_number()) %>% ungroup() %>%
    pivot_wider(id_cols = .row, names_from = model, values_from = value) %>%
    select(-.row) %>%
    drop_na()  
  
  mat <- as.matrix(wide)                      # iterations × parameters
  colnames(mat) <- colnames(wide)             # e.g., "IMR", "U5MR"
  
  regs <- bayesplot::mcmc_areas(mat, prob = 0.95) +
    labs(x = "Effects", y = NULL, title = "Regression coefficients")
  
  out <- patchwork::wrap_plots(regs, thresh)
}

figs_dir <- "../../results/figs"
pdf(file.path(figs_dir, "regression_densities.pdf"), 
    width = 7, height = 5)
print(reg_dens())
dev.off()

#=======================================================================================
# Some diagnostics
#=======================================================================================

few.diagnostics <- function(fit){
  # Check standard Stan diagnostics
  mod_diagnostics <- fit$diagnostic_summary()
  cat("Divergent transitions:", mod_diagnostics$num_divergent, "\n")
  cat("Max treedepth hits:", mod_diagnostics$num_max_treedepth, "\n")
  
  # Detailed convergence diagnostics
  mod_summary <- fit$summary(variables = c("threshold", "beta_trt", "beta_trt_mortality"))
  print(mod_summary[c("variable", "rhat", "ess_bulk", "ess_tail")])
}

few.diagnostics(fit_u5mr)
few.diagnostics(fit_imr)
few.diagnostics(fit_u5mr_fixed)
few.diagnostics(fit_imr_fixed)

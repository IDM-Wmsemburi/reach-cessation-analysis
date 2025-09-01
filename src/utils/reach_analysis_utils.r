# ===============================================================================
# REACH Secondary Analysis: Core Utility Functions
# ===============================================================================
# Purpose: Shared utility functions used across multiple analysis scripts
# Author: William Msemburi
# Date: August 2025
# Note: Contains ONLY functions that are called in multiple pipeline stages
# ===============================================================================

# Required packages for the utility functions
required_packages <- c("data.table", "dplyr", "sf")

# Load required packages quietly
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not found: ", pkg)
  }
}

# ===============================================================================
# DATA AGGREGATION FUNCTIONS (used in contextual analysis & thresholds)
# ===============================================================================

#' Aggregate deaths and person time by grouping variables
#' Used in: contextual_heterogeneity.R, threshold_analysis.R
#' @param data Input dataset
#' @param grouping_vars Variables to group by
#' @param death_var Name of death count variable (default: "n")
#' @param person_time_var Name of person time variable (default: "person_time")
#' @return Aggregated dataset with deaths and person_time columns
sum_deaths_and_person_time <- function(data, grouping_vars, death_var = "n", person_time_var = "person_time") {
  dt <- data.table::as.data.table(data)
  result <- dt[, .(deaths = sum(get(death_var), na.rm = TRUE),
                   person_time = sum(get(person_time_var), na.rm = TRUE)), 
               by = grouping_vars]
  dplyr::as_tibble(result)
}

# ===============================================================================
# INDEX CREATION FUNCTIONS (used in baseline mortality & thresholds)
# ===============================================================================

#' Create stable group indices for trial-country combinations
#' Used in: bayesian_hsgp_model.R, threshold_analysis.R
#' @param df Dataset with trial_phase and country columns
#' @return Tibble with group, trial_phase, country, and group_id
create_group_index <- function(df) {
  df %>%
    dplyr::transmute(group = paste(trial_phase, country, sep = " "), trial_phase, country) %>%
    dplyr::distinct() %>% 
    dplyr::arrange(group) %>% 
    dplyr::mutate(group_id = row_number())
}

#' Create stable country indices
#' Used in: bayesian_hsgp_model.R
#' @param df Dataset with country column
#' @return Tibble with country and country_id
create_country_index <- function(df) {
  df %>% 
    dplyr::distinct(country) %>% 
    dplyr::arrange(country) %>% 
    dplyr::mutate(country_id = row_number())
}

# ===============================================================================
# SPATIAL UTILITY FUNCTIONS (used in baseline mortality)
# ===============================================================================

#' Pick appropriate UTM EPSG code based on coordinates
#' Used in: bayesian_hsgp_model.R (for HSGP spatial basis)
#' @param lon Longitude
#' @param lat Latitude  
#' @return EPSG code for UTM projection
pick_utm_epsg <- function(lon, lat) {
  zone <- floor((lon + 180)/6) + 1
  if (lat >= 0) 32600 + zone else 32700 + zone
}

# ===============================================================================
# STAN MODEL UTILITIES (used in threshold analysis)
# ===============================================================================

# Safe coefficient extraction for model initialization
get_or0 <- function(coef_vec, coef_name) {
  if (coef_name %in% names(coef_vec)) coef_vec[[coef_name]] else 0
}

compute_smart_inits_single <- function(data, mortality_col, n_chains = 4) {
  # GLM fit for initialization
  fml <- as.formula(paste("deaths ~ offset(log(person_time)) +", mortality_col, "+ trt + factor(group_id)"))
  
  tryCatch({
    glm_fit <- glm(fml, data = data, family = poisson())
    coefs <- coef(glm_fit)
  }, error = function(e) {
    warning("GLM initialization failed, using default values")
    coefs <- c("(Intercept)" = 0, trt = 0)
    names(coefs)[2] <- mortality_col
  })
  
  purrr::map(1:n_chains, function(i) {
    list(
      log_mortality_raw = rnorm(nrow(data), 0, 0.1),
      alpha_group_raw = rnorm(max(data$group_id), 0, 0.1),
      sigma_group = abs(rnorm(1, 0.5, 0.1)),
      mu_alpha_group = get_or0(coefs, "(Intercept)") + rnorm(1, 0, 0.1),
      beta_mortality = get_or0(coefs, mortality_col) + rnorm(1, 0, 0.1),
      beta_trt = get_or0(coefs, "trt") + rnorm(1, 0, 0.1),
      beta_trt_mortality = rnorm(1, 0, 0.1)
    )
  })
}


#' Generic Stan data construction for mortality models (with MC integration)
#' Used in: threshold_analysis.R
#' @param df Input dataset
#' @param log_mean_col Column name for log mortality mean
#' @param log_sd_col Column name for log mortality standard deviation
#' @param xmin Minimum value for prediction grid
#' @param xmax Maximum value for prediction grid
#' @param n_pred Number of prediction points (default: 41)
make_stan_data <- function(df, log_mean_col, log_sd_col, xmin, xmax, n_pred = 41L) {
  list(
    N = nrow(df),
    deaths = as.integer(df$deaths),
    person_time = df$person_time,
    trt = as.integer(df$trt),
    group = as.integer(df$group_id),
    J = as.integer(max(df$group_id)),
    log_mortality_mean = df[[log_mean_col]],
    log_mortality_sd = df[[log_sd_col]],
    N_pred = n_pred,
    log_mortality_pred = log(seq(xmin, xmax, length.out = n_pred))
  )
}

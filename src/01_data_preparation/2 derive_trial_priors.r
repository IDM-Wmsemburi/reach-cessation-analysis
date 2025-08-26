# ===============================================================================
# REACH Secondary Analysis: Trial Prior Derivation
# ===============================================================================
# Purpose: Generate prior distributions from REACH trial placebo data and DHS
# Author: William Msemburi
# Date: August 2025
# ===============================================================================

# Environment setup
library(tidyverse)
library(dplyr)
library(tidyr)
library(purrr)
library(sandwich)  # for vcovCL

rm(list = ls())  # Clear environment

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

# Load utility functions
source("../utils/reach_analysis_utils.R")

# ===============================================================================
# STEP 1: Link H1_11 & H12_59 to H0 by country using DHS data
# ===============================================================================

#' Fit pooled H0 model using DHS hazard relationships
#' @return Tibble with pooled model coefficients and metadata
get_dhs_h0_fits <- function() {
  # Load DHS hazard data
  if (!file.exists("../../results/dhs_hazard_data.rda")) {
    stop("DHS hazard data not found. Run prepare_dhs_hazards.R first.")
  }
  
  load("../../results/dhs_hazard_data.rda")
  rm(dhs.ests.df)  # Only need hazards.df
  
  # Process hazards summary
  hazards_summary <- hazards.df |>
    dplyr::mutate(
      hsi = hs * ns, 
      id = substr(isid, 1, 2),
      country = dplyr::case_when(
        id == "BF" ~ "Burkina Faso", 
        id == "MW" ~ "Malawi", 
        id == "NG" ~ "Nigeria", 
        id == "NI" ~ "Niger", 
        id == "TZ" ~ "Tanzania"
      )
    ) |> 
    dplyr::filter(region == "All") |> 
    dplyr::group_by(isid, country, region, years) |>
    dplyr::summarise(
      H_0    = sum(hsi[age_group == "0"]), 
      H_111  = sum(hsi[age_group == "1-11"]),
      H_1259 = sum(hsi[age_group %in% c("12-23", "24-35", "36-47", "48-59")]), 
      .groups = "drop"
    )
  
  # Fit pooled model with country fixed effects
  fit_pooled_h0_model <- function(df) {
    # Pooled model with country fixed effects
    pooled_fit <- lm(log(H_0) ~ log(H_111) + log(H_1259) + factor(country), data = df)
    
    # Extract coefficients
    coefs <- coef(pooled_fit)
    country_names <- c("Burkina Faso", "Malawi", "Niger", "Nigeria", "Tanzania")
    
    # Create country-specific intercepts
    base_intercept <- coefs["(Intercept)"]
    country_intercepts <- c(
      base_intercept,  # Burkina Faso (reference)
      base_intercept + ifelse("factor(country)Malawi" %in% names(coefs), coefs["factor(country)Malawi"], 0),
      base_intercept + ifelse("factor(country)Niger" %in% names(coefs), coefs["factor(country)Niger"], 0),
      base_intercept + ifelse("factor(country)Nigeria" %in% names(coefs), coefs["factor(country)Nigeria"], 0),
      base_intercept + ifelse("factor(country)Tanzania" %in% names(coefs), coefs["factor(country)Tanzania"], 0)
    )
    
    # Return in structured format for Stan compatibility
    tibble::tibble(
      model_type = "pooled_h0",
      beta_h111  = coefs["log(H_111)"],
      beta_h1259 = coefs["log(H_1259)"],
      sigma_hat  = sigma(pooled_fit),
      r2         = 100 * summary(pooled_fit)$r.squared,
      obs        = nrow(df),
      country    = list(country_names),
      alpha_hat  = list(country_intercepts),
      vcov       = list(vcov(pooled_fit)),
      model      = list(pooled_fit)
    )
  }
  
  fit_pooled_h0_model(hazards_summary)
}

# ===============================================================================
# STEP 2: Generate hazard priors from REACH trial placebo data
# ===============================================================================

#' Extract hazard priors from trial placebo data
#' @return Tibble with trial-country specific hazard priors
get_hazard_priors <- function() {
  # Load trial data
  if (!file.exists("../../data/Clean/data_crt_geospatial.rda")) {
    stop("Clean trial data not found. Please ensure data_crt_geospatial.rda exists.")
  }
  
  load("../../data/Clean/data_crt_geospatial.rda")
  
  # Prepare trial data - placebo only for priors
  trial_data <- data_crt_geospatial |>
    dplyr::filter(
      !is.na(trt), !is.na(latitude), !is.na(longitude), !is.na(person_time),
      trt == 0,  # Placebo only for priors
      age_start >= 1, age_start <= 59,
      !is.na(n), person_time > 0
    ) |> 
    dplyr::mutate(
      person_time_years = person_time/365.25,
      trial_country = paste(trial_phase, country)
    )
  
  message("Data prepared: ", nrow(trial_data), " observations")
  
  # Function to fit model for specific age range and trial-country
  fit_age_stratum_model <- function(data, age_min, age_max, trial_ctry, age_name) {
    # Filter for specific age range and trial-country
    subset_data <- data %>%
      dplyr::filter(
        trial_country == trial_ctry,
        age_start >= age_min, 
        age_start <= age_max
      )
    
    # Aggregate by cluster (community_id)
    cluster_data <- subset_data %>%
      dplyr::group_by(community_id) %>%
      dplyr::summarise(
        deaths = sum(n),
        person_years = sum(person_time_years),
        .groups = "drop"
      ) %>%
      dplyr::filter(person_years > 0) %>%
      dplyr::mutate(log_offset = log(person_years))
    
    # Fit quasipoisson model
    model <- glm(
      deaths ~ 1 + offset(log_offset), 
      data = cluster_data, 
      family = quasipoisson()
    )
    
    # Get cluster-robust standard errors (analytical)
    vcov_cl <- sandwich::vcovCL(model, cluster = cluster_data$community_id)
    se_cl <- sqrt(diag(vcov_cl))
    
    # Extract intercept (log hazard rate per person-year)
    intercept_coef <- coef(model)["(Intercept)"]
    intercept_se <- se_cl["(Intercept)"]
    
    # Convert to cumulative hazard over age period
    time_period <- (age_max - age_min + 1) / 12  # Convert months to years
    log_cum_hazard_mean <- intercept_coef + log(time_period)
    log_cum_hazard_sd <- intercept_se  # SE of log doesn't change with additive constant
    
    tibble::tibble(
      trial_country = trial_ctry,
      age_stratum = age_name,
      log_h_mean = log_cum_hazard_mean,
      log_h_sd = log_cum_hazard_sd,
      method = "quasipoisson_analytical",
      n_obs = nrow(subset_data),
      n_clusters = nrow(cluster_data)
    )
  }
  
  # Get all unique trial-countries
  trial_countries <- trial_data %>%
    dplyr::distinct(trial_country) %>%
    dplyr::pull(trial_country)
  
  # Age stratum definitions
  age_strata <- list(
    list(age_min = 1, age_max = 11, age_name = "1-11"),
    list(age_min = 12, age_max = 59, age_name = "12-59"),
    list(age_min = 1, age_max = 59, age_name = "1-59")
  )
  
  message("Fitting models for ", length(trial_countries), " trial-countries and ", length(age_strata), " age strata...")
  
  # Fit all models
  successful_results <- purrr::map_dfr(trial_countries, function(tc) {
    purrr::map_dfr(age_strata, function(stratum) {
      message("  Fitting: ", tc, " - ", stratum$age_name, " months")
      fit_age_stratum_model(
        trial_data, 
        stratum$age_min, 
        stratum$age_max, 
        tc, 
        stratum$age_name
      )
    })
  })
  
  # Reshape to wide format for Stan
  wide_results <- successful_results %>%
    dplyr::select(trial_country, age_stratum, log_h_mean, log_h_sd) %>%
    tidyr::pivot_wider(
      names_from = age_stratum,
      values_from = c(log_h_mean, log_h_sd),
      names_sep = "_"
    )
  
  # Standardize column names for consistency
  names(wide_results) <- gsub("1-11", "1_11", names(wide_results))
  names(wide_results) <- gsub("12-59", "12_59", names(wide_results))
  names(wide_results) <- gsub("1-59", "1_59", names(wide_results))
  
  # Final formatting - direct rename using expected column names
  final_results <- wide_results %>%
    dplyr::rename(
      log_h1_11_prior = log_h_mean_1_11,
      log_h1_11_prior_sd = log_h_sd_1_11,
      log_h12_59_prior = log_h_mean_12_59,
      log_h12_59_prior_sd = log_h_sd_12_59,
      log_h1_59_prior = log_h_mean_1_59,
      log_h1_59_prior_sd = log_h_sd_1_59
    ) 
  
  return(final_results)
}

# ===============================================================================
# MAIN EXECUTION
# ===============================================================================

cat("Starting trial prior derivation...\n")

# Step 1: Generate DHS-based H0 model fits
cat("Fitting DHS H0 pooled model...\n")
system.time({
  fits <- get_dhs_h0_fits()
})

cat("DHS model summary:\n")
cat("  R-squared: ", round(fits$r2, 1), "%\n")
cat("  Observations: ", fits$obs, "\n")
cat("  Residual SD: ", round(fits$sigma_hat, 3), "\n")

# Step 2: Generate hazard priors using analytical method
cat("Extracting hazard priors from trial placebo data...\n")
system.time({
  final_priors <- get_hazard_priors()
})

cat("Trial prior summary:\n")
print(final_priors)

# ===============================================================================
# SAVE PROCESSED MODELS
# ===============================================================================

# Build compact, explicit object for the pooled DHS model
pooled_model <- {
  stopifnot("model" %in% names(fits), length(fits$model) >= 1L)
  m <- fits$model[[1]]
  stopifnot(inherits(m, "lm"))
  co <- coef(m)
  vc <- vcov(m)
  mf <- model.frame(m)
  
  list(
    schema_version = "h0-pool-v1",
    call           = m$call,
    coef           = co,            # (Intercept), log(H_111), log(H_1259), FEs
    vcov           = vc,            # full VCOV over those terms
    terms          = names(co),     # to align indices robustly later
    sigma_hat      = stats::sigma(m),      # residual SD on log(H0)
    r2             = unname(summary(m)$r.squared),
    countries_seen = sort(unique(mf$country)),
    reference_level = levels(mf$country)[1],
    created_at     = format(Sys.time(), tz = "UTC")
  )
}

# Save results
saveRDS(pooled_model, file = "../../results/dhs_h0_pooled_model.rds")
saveRDS(final_priors, file = "../../results/trial_hazard_priors.rds")

# ===============================================================================
# VALIDATION AND SUMMARY
# ===============================================================================

cat("\n=== PRIOR DERIVATION SUMMARY ===\n")
cat("DHS H0 Pooled Model:\n")
cat("  Schema version:", pooled_model$schema_version, "\n")
cat("  Countries included:", paste(pooled_model$countries_seen, collapse = ", "), "\n")
cat("  R-squared:", round(pooled_model$r2 * 100, 1), "%\n")
cat("  Residual SD:", round(pooled_model$sigma_hat, 3), "\n")

cat("\nTrial Hazard Priors:\n")
cat("  Trial-countries:", nrow(final_priors), "\n")
cat("  Age strata: 1-11m, 12-59m, 1-59m\n")

# Display prior ranges
prior_summary <- final_priors %>%
  dplyr::select(trial_country, contains("_prior"), -contains("_sd")) %>%
  dplyr::mutate(
    h1_11_median = exp(log_h1_11_prior),
    h12_59_median = exp(log_h12_59_prior),
    h1_59_median = exp(log_h1_59_prior)
  ) %>%
  dplyr::select(trial_country, h1_11_median, h12_59_median, h1_59_median)

cat("\nHazard prior medians (cumulative hazards):\n")
print(prior_summary)

cat("\nFiles saved:\n")
cat("  - ../../results/dhs_h0_pooled_model.rds\n")
cat("  - ../../results/trial_hazard_priors.rds\n")

message("Prior derivation completed successfully!")
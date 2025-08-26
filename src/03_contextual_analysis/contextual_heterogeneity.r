# ===============================================================================
# REACH Secondary Analysis: Contextual Heterogeneity
# ===============================================================================
# Purpose: Explore vaccination coverage, malaria burden, and other contextual factors
# Author: William Msemburi
# Date: August 2025
# ===============================================================================

rm(list = ls())  # Clear environment

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

# Load required libraries
library(tidyverse)
library(lme4)
library(future)
library(furrr)
library(broom.mixed)
library(raster)
library(sf)
library(corrplot)
library(scales)
library(ggthemes)

# Load utility functions
source("../utils/reach_analysis_utils.R")

# ===============================================================================
# DATA LOADING AND PREPARATION
# ===============================================================================

cat("Loading data...\n")

# Load trial data and baseline mortality estimates
if (!file.exists("../../data/Clean/data_crt_geospatial.rda")) {
  stop("Clean trial data not found")
}
load("../../data/Clean/data_crt_geospatial.rda")

if (!file.exists("../../results/baseline_mortality_estimates.rda")) {
  stop("Baseline mortality estimates not found. Run bayesian_hsgp_model.R first.")
}
load("../../results/baseline_mortality_estimates.rda")

# Filter and prepare analysis dataset
data_analysis <- data_crt_geospatial |>
  dplyr::filter(!is.na(trt), 
                !is.na(latitude), !is.na(longitude),  
                !is.na(person_time))

cat("Data loaded: ", nrow(data_analysis), " records\n")

# ===============================================================================
# EXTRACT CONTEXTUAL COVARIATES FROM RASTERS
# ===============================================================================

cat("Extracting contextual covariates...\n")

# Define paths to raster files (adjust paths as needed)
raster_paths <- list(
  bcg1   = "../../data/IHME/bcg1_cov_mean_raked_2000_2024.TIF",
  mcv1   = "../../data/IHME/mcv1_cov_mean_raked_2000_2024.TIF", 
  dpt1   = "../../data/IHME/dpt1_cov_mean_raked_2000_2024.TIF",
  dpt3   = "../../data/IHME/dpt3_cov_mean_raked_2000_2024.TIF",
  polio3 = "../../data/IHME/polio3_cov_mean_raked_2000_2024.TIF",
  itn_use = "../../data/Malaria/2024_GBD2023_Africa_ITN_Use_Rate_2019.TIF",
  itn_acc = "../../data/Malaria/2024_GBD2023_Africa_ITN_Access_2019.TIF",
  mal_treat = "../../data/Malaria/2024_GBD2023_Global_Antimalarial_EFT_2019.TIF",
  mal_inc = "../../data/Malaria/2024_GBD2023_Global_Pf_Incidence_Rate_2019.TIF",
  mal_mort = "../../data/Malaria/2024_GBD2023_Global_Pf_Mortality_Rate_2019.TIF"
)

# Function to safely load and extract raster values
extract_raster_values <- function(data, raster_path, var_name, as_percentage = TRUE, invert = FALSE) {
  if (!file.exists(raster_path)) {
    warning("Raster file not found: ", raster_path)
    return(rep(NA, nrow(data)))
  }
  
  tryCatch({
    rast <- raster::raster(raster_path)
    target_sf <- sf::st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)
    values <- raster::extract(rast, target_sf)
    
    if (invert) {
      values <- 1 - values
    }
    
    if (as_percentage) {
      values <- 100 * values
    }
    
    return(values)
  }, error = function(e) {
    warning("Error extracting ", var_name, ": ", e$message)
    return(rep(NA, nrow(data)))
  })
}

# Extract vaccination coverage (as unmet need - inverted)
data_analysis$bcg1   <- extract_raster_values(data_analysis, raster_paths$bcg1, "BCG1", TRUE, TRUE)
data_analysis$polio3 <- extract_raster_values(data_analysis, raster_paths$polio3, "POLIO3", TRUE, TRUE)
data_analysis$mcv1   <- extract_raster_values(data_analysis, raster_paths$mcv1, "MCV1", TRUE, TRUE)
data_analysis$dpt1   <- extract_raster_values(data_analysis, raster_paths$dpt1, "DPT1", TRUE, TRUE)
data_analysis$dpt3   <- extract_raster_values(data_analysis, raster_paths$dpt3, "DPT3", TRUE, TRUE)

# Extract malaria indicators
data_analysis$itn.use    <- extract_raster_values(data_analysis, raster_paths$itn_use, "ITN use", TRUE, TRUE)
data_analysis$itn.acc    <- extract_raster_values(data_analysis, raster_paths$itn_acc, "ITN access", TRUE, TRUE)
data_analysis$mal.treat  <- extract_raster_values(data_analysis, raster_paths$mal_treat, "Malaria treatment", TRUE, TRUE)
data_analysis$mal.inc    <- extract_raster_values(data_analysis, raster_paths$mal_inc, "Malaria incidence", FALSE, FALSE)
data_analysis$mal.mort   <- extract_raster_values(data_analysis, raster_paths$mal_mort, "Malaria mortality", FALSE, FALSE)

# Calculate derived vaccination variables
data_analysis <- data_analysis |>
  dplyr::mutate(
    # DTP dropout measures
    dpt.rel.drop = 100 * ((100 - dpt1) - (100 - dpt3))/(100 - dpt1),
    dpt.abs.drop = (dpt3 - dpt1),
    
    # Geometric mean coverage (unmet)
    gm_cov = exp(rowMeans(log(data_analysis[, c("bcg1", "dpt1", "dpt3", "mcv1", "polio3")]), na.rm = TRUE)),
    
    # Log-transform baseline U5MR from IHME
    ihme_u5m_2015_log = log(ihme_u5m_2015)
  )

# Principal components analysis for vaccination
pc1 <- prcomp(data_analysis[, c("bcg1", "dpt1", "dpt3", "mcv1", "polio3")], scale. = TRUE)
data_analysis$pc1 <- pc1$x[,1]

# Join with baseline mortality estimates
data_analysis <- data_analysis |> 
  dplyr::left_join(df.mod |> 
                   dplyr::select(community_id, trial_phase, country, 
                                post_u5mr, post_u5mr_log, post_imr, post_imr_log),
                   by = c("community_id", "trial_phase", "country")) |>
  dplyr::rename(post_u5m_log = post_u5mr_log)

cat("Contextual covariates extracted\n")

# ===============================================================================
# DEFINE COVARIATE METADATA
# ===============================================================================

# Define covariate names and labels for analysis
covariate_names <- c(
  "ihme_u5m_2015_log", "post_u5m_log", "post_imr_log", 
  "dpt1", "dpt3", "dpt.rel.drop", "dpt.abs.drop", 
  "bcg1", "polio3", "mcv1", "gm_cov", "pc1", 
  "mal.mort", "mal.inc", "itn.acc", "itn.use"
)

covariate_titles <- c(
  "IHME U5M Values", "Posterior U5M Values","Posterior IMR Values", 
  "DPT1 unmet", "DPT3 unmet", "Rel DPT1-3 dropout", "Abs DPT1-3 dropout", 
  "BCG unmet", "POLIO3 unmet", "MCV1 unmet", "Unmet geo mean", "Unmet PC1", 
  "PF mortality", "PF incidence", "ITN access", "ITN use"
)

covariate_short_names <- c(
  "ihme_u5m_log", "post_u5m_log", "post_imr_log", 
  "dpt1", "dpt3", "dpt.rel.drop", "dpt.abs.drop", 
  "bcg1", "polio3", "mcv1", "gm_cov", "pc1", 
  "mal.mort", "mal.inc", "itn.acc", "itn.use"
)

covariate_short_titles <- c(
  "IHME U5M (log)", "Posterior U5M (log)", "Posterior IMR (log)", 
  "DPT1 unmet", "DPT3 unmet", "Rel DPT1-3 dropout", "Abs DPT1-3 dropout", 
  "BCG unmet", "POLIO3 unmet", "MCV1 unmet", "Unmet geo mean", "Unmet PC1", 
  "PF mortality", "PF incidence", "ITN access", "ITN use"
)

# Helper functions for covariate labels
make_short_name <- function(c_name) {
  covariate_short_names[match(c_name, covariate_names)]
}

make_short_title <- function(c_name) {
  covariate_short_titles[match(c_name, covariate_names)]
}

make_long_title <- function(c_name) {
  covariate_titles[match(c_name, covariate_names)]
}

# ===============================================================================
# UNIVARIATE SCATTER PLOTS AND CORRELATION ANALYSIS
# ===============================================================================

cat("Creating scatter plots and correlation analysis...\n")

# Create summary dataset
data_summary <- data_analysis |> 
  sum_deaths_and_person_time(c("country", "trial_phase", "community_id", "trt", covariate_names)) |>
  dplyr::ungroup()

# Function to create binned data for scatter plots
make_binned_data <- function(data, covariate_name, covariate_title, width_divide = 10) {
  data <- data |>
    dplyr::mutate(trt = factor(trt, levels = c(0, 1), labels = c("Control", "Azithromycin"))) |>
    dplyr::filter(!is.na(trt))
  
  width_use <- data |>
    dplyr::group_by(country, trial_phase) |>
    dplyr::summarize(covariate_width = diff(c(min(.data[[covariate_name]], na.rm = TRUE),
                                             quantile(.data[[covariate_name]], 0.95, na.rm = TRUE))),
                     covariate_width = covariate_width / width_divide,
                     .groups = "drop")
  
  data |> 
    dplyr::inner_join(width_use, by = c("country", "trial_phase")) |>
    dplyr::group_by(country, trial_phase) |>
    dplyr::mutate(
      var = covariate_name, 
      binned_covariate = ifelse(var != "mcv1_group", 
                               cut_width(.data[[covariate_name]], width = covariate_width[1]), 
                               as.numeric(.data[[covariate_name]]))
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(trt, country, trial_phase, binned_covariate) |>
    dplyr::summarise(deaths = sum(deaths),
                    person_time = sum(person_time),
                    covariate_mean = mean(.data[[covariate_name]]),
                    nclusters = n(),
                    .groups = "drop") |>
    dplyr::mutate(rate = deaths * 365 * 1000 / person_time, 
                  variable = covariate_name)
}

# Create binned scatter data
binned_scatter_data <- purrr::map2(covariate_names, covariate_titles,
                                   ~make_binned_data(data_summary, .x, .y))

# Function to create scatter plots
make_scatter_plots <- function(df.plot, covariate_name, covariate_title) {
  output_dir <- "../../results/figs"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (! covariate_name %in% c("ihme_u5m_2015_log", "post_u5m_log", "post_imr_log")) {
    covariate_title <- paste(covariate_title)
  }
  
  png(paste0(output_dir, "/scatter_", covariate_name, ".png"), 
      width = 2000, height = 1500, res = 300)
  
  p <- ggplot(df.plot, aes(x = covariate_mean, y = rate)) +
    geom_point(aes(col = trt, size = person_time / (365 * 1000)), alpha = 0.4) +
    labs(x = covariate_title,
         y = "Deaths / Person Time (1000 yr)") +
    facet_wrap(~paste0(country, "\n", trial_phase), scale = "free") +
    guides(size = guide_legend(title = "Person time"),
           col = guide_legend(title = "Treatment")) +
    theme_bw() +
    theme(legend.position = "right")
  
  print(p)
  dev.off()
}

# Create all scatter plots
for (j in 1:length(covariate_names)) {
  make_scatter_plots(binned_scatter_data[[j]], 
                     covariate_names[j],
                     covariate_titles[j])
}

# ===============================================================================
# CORRELATION ANALYSIS
# ===============================================================================

cat("Computing correlation matrix...\n")

# Prepare correlation data
cor.data <- list()
for (j in 1:length(covariate_names)) {
  cor.data[[covariate_names[j]]] <- binned_scatter_data[[j]] 
}
cor.data <- dplyr::bind_rows(cor.data) 

# Create correlation matrix
cor.data <- cor.data %>%
  dplyr::filter(variable %in% covariate_names) %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(obs_id = row_number()) %>%
  dplyr::ungroup()

# Pivot to wide format
covariate_wide <- cor.data %>%
  dplyr::select(obs_id, variable, covariate_mean) %>%
  tidyr::pivot_wider(names_from = variable, values_from = covariate_mean)

rate_wide <- cor.data %>%
  dplyr::select(obs_id, variable, rate) %>%
  tidyr::pivot_wider(names_from = variable, values_from = rate)

# Average rate across variables
avg_rate <- rate_wide %>%
  dplyr::mutate(rate = rowMeans(dplyr::across(dplyr::all_of(covariate_names)), na.rm = TRUE)) %>%
  dplyr::select(obs_id, rate)

# Create full matrix for correlation
full_matrix <- covariate_wide %>%
  dplyr::inner_join(avg_rate, by = "obs_id") %>%
  dplyr::select(rate, dplyr::all_of(covariate_names)) |>
  dplyr::rename(death_rate = rate, ihme_u5m_log = ihme_u5m_2015_log)

# Compute correlation matrix
cor_matrix <- cor(full_matrix, use = "pairwise.complete.obs", method = "spearman")

# Create correlation plot
output_dir <- "../../results/figs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

png(paste0(output_dir, "/correlation_matrix.png"), width = 1500, height = 1500, res = 300)
par(mar = c(0, 0, 0, 0))
corrplot::corrplot(cor_matrix, method = "pie", type = "upper",
                   order = "original",
                   tl.cex = 0.8, cl.cex = 0.8,
                   number.cex = 0.7,
                   col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(200))
dev.off()

# Variables to use in further regression

cov_subset <- c(
  "ihme_u5m_2015_log", "post_u5m_log", "post_imr_log", 
  "dpt1", "bcg1", "mcv1", 
  "mal.mort", "mal.inc", "itn.acc"
)

# ===============================================================================
# MIXED EFFECTS MODELS WITH AGE STRATIFICATION
# ===============================================================================

cat("Running mixed effects models with age stratification...\n")

# Define model formulas
formula_glmer_full <- deaths ~ trt*covariate + offset(log(person_time)) + (1|community_id_factor)
formula_glmer_default <- deaths ~ trt + offset(log(person_time)) + (1|community_id_factor)
formula_glmer_partial <- deaths ~ trt + covariate + offset(log(person_time)) + (1|community_id_factor)

# Add age splines to original dataset
data_summary_age_spline <- data_analysis |> 
  dplyr::mutate(ihme_u5m_2015_log = log(ihme_u5m_2015)) |>
  dplyr::left_join(df.mod |> 
                   dplyr::select(community_id, trial_phase, country, post_u5mr, post_u5mr_log),
                   by = c("community_id", "trial_phase", "country")) |>  
  sum_deaths_and_person_time(c("trial_phase", "country", "community_id", "trt", "age_start", 
                               cov_subset)) |>
  dplyr::ungroup() |>
  dplyr::filter(age_start >= 1 & age_start <= 59, !is.na(trt), !is.na(ihme_u5m_2015_log)) %>%
  {
    df <- .
    bind_cols(df,
              tibble::as_tibble(splines::ns(df$age_start, 
                                           knots = c(5.5, 11.5, 23.5, 41.5),
                                           Boundary.knots = c(1,59)) |>
                               as.matrix()) |>
              purrr::set_names(paste0("age_ns_", 1:5)) |>
              dplyr::mutate(dplyr::across(starts_with("age_ns"), ~as.numeric(.x))))
  }

# Function to add splines to formulas
add_splines_to_formula <- function(form) {
  update(form, . ~ . + age_ns_1 + age_ns_2 + age_ns_3 + age_ns_4 + age_ns_5)
}

# Set up parallel processing
plan(multisession, workers = min(8, length(cov_subset)))

# Create analysis grid
me_age_indices <- tidyr::expand_grid(
  covariate_name = cov_subset,
  data = data_summary_age_spline |> dplyr::group_by(country, trial_phase) |> dplyr::group_split(),
  tibble::tibble(
    formula = list(
      model_full = formula_glmer_full,
      model_partial = formula_glmer_partial,
      model_default = formula_glmer_default
    ) |> purrr::map(add_splines_to_formula),
    formula_name = c("model_full", "model_partial", "model_default")
  )
)

# Function to fit mixed effects models
fit_me_model <- function(.covariate_name, .data, .formula, .formula_name) {
  .data[["covariate"]] <- .data[[.covariate_name]]
  .data <- .data |>
    dplyr::filter(!is.na(trt), !is.na(covariate)) |>
    dplyr::mutate(community_id_factor = factor(community_id))
  
  tryCatch({
    lme4::glmer(.formula, data = .data, family = poisson())
  }, error = function(e) {
    warning("Model fitting failed for ", .covariate_name, ": ", e$message)
    return(NULL)
  })
}

# Fit all models in parallel
plan(multisession, workers = parallel::detectCores() - 1)

mixed_effects_age <- furrr::future_pmap(
  with(me_age_indices, list(covariate_name, data, formula, formula_name)),
  fit_me_model,
  .progress = TRUE
)

# Process results
mixed_effects_age <- me_age_indices |>
  dplyr::mutate(
    country = purrr::map_chr(data, ~.x$country[1]),
    trial_phase = purrr::map_chr(data, ~.x$trial_phase[1])
  ) |>
  dplyr::select(-data, -formula) |>
  dplyr::mutate(
    covariate_title = make_long_title(covariate_name),
    model = mixed_effects_age
  ) |>
  tidyr::pivot_wider(names_from = formula_name, values_from = model)

# Remove intermediate objects
rm(me_age_indices)

# Run ANOVA comparisons
mixed_effects_age <- mixed_effects_age |>
  dplyr::rowwise() |>
  dplyr::mutate(
    anova_res = list({
      if (!is.null(model_full) && !is.null(model_partial) && !is.null(model_default)) {
        tryCatch({
          anova(model_full, model_partial, model_default)
        }, error = function(e) {
          warning("ANOVA failed: ", e$message)
          return(NULL)
        })
      } else {
        NULL
      }
    })
  ) |>
  dplyr::ungroup()

# ===============================================================================
# PROCESS ANOVA RESULTS
# ===============================================================================

cat("Processing ANOVA results...\n")

# Extract ANOVA results
mixed_effects_age_anova <- mixed_effects_age |>
  dplyr::select(country, trial_phase, covariate_name, covariate_title, anova_res) |>
  dplyr::filter(!purrr::map_lgl(anova_res, is.null)) |>
  dplyr::mutate(anova_res = purrr::map(anova_res, ~broom::tidy(.x))) |>
  tidyr::unnest(anova_res) |>
  dplyr::group_by(country, trial_phase, covariate_name, covariate_title) |>
  dplyr::mutate(
    AIC_diff = dplyr::lag(AIC) - AIC,
    BIC_diff = dplyr::lag(BIC) - BIC,
    chisq = statistic,
    p.value = p.value,
    contrast = c("Base", "Covariate vs. Base", "Cov-Trt Interaction vs. Cov. only")
  ) |>
  dplyr::filter(!is.na(AIC_diff)) |>
  dplyr::ungroup() |>
  dplyr::mutate(contrast = factor(contrast, levels = c("Covariate vs. Base",
                                                       "Cov-Trt Interaction vs. Cov. only")))

# Function to bin p-values for visualization
sigbins <- function(df) {
  df |> dplyr::mutate(sig_bin = dplyr::case_when(
    p.value < 0.01  ~ "<0.01",
    p.value < 0.05  ~ "<0.05",
    p.value < 0.10  ~ "<0.10",
    p.value < 0.15  ~ "<0.15",
    TRUE            ~ "≥0.15"
  ))
}

# Create heatmap data
heatmap_data1 <- mixed_effects_age_anova |>
  dplyr::filter(contrast == "Covariate vs. Base") |>
  dplyr::mutate(
    trial_study = paste(country, trial_phase, sep = " – "),
    covariate_short = make_short_title(covariate_name)
  ) |> 
  sigbins() |> 
  dplyr::mutate(type = "Covariate")

heatmap_data2 <- mixed_effects_age_anova |>
  dplyr::filter(contrast == "Cov-Trt Interaction vs. Cov. only") |>
  dplyr::mutate(
    trial_study = paste(country, trial_phase, sep = " – "),
    covariate_short = make_short_title(covariate_name)
  ) |> 
  sigbins() |> 
  dplyr::mutate(type = "Cov-Trt Interaction")

# ===============================================================================
# CREATE HEATMAP VISUALIZATION
# ===============================================================================

cat("Creating heatmap visualization...\n")

# Color palette for significance levels
sig_colors <- c(
  "<0.01"  = "#08519C",  # strong blue
  "<0.05"  = "#3182BD",  # medium-dark blue
  "<0.10"  = "#6BAED6",  # medium blue
  "<0.15"  = "#C6DBEF",  # light blue
  "≥0.15"  = "white"
)

# Combine heatmap data
heatmaps <- rbind(heatmap_data1, heatmap_data2) |>
  dplyr::mutate(
    type = factor(type, levels = c("Covariate", "Cov-Trt Interaction")), 
    covariate_short = factor(covariate_short, 
                            levels = c("IHME U5M (log)", "Posterior U5M (log)", 
                                       "Posterior IMR (log)",
                                     "DPT1 unmet", "BCG unmet", "MCV1 unmet", 
                                     "PF mortality", "PF incidence", 
                                     "ITN access"))
  )

# Function to create heatmap
get_heatmaps <- function(df) {
  ggplot(df, aes(x = covariate_short, y = forcats::fct_rev(trial_study), fill = sig_bin)) +
    geom_tile(color = "grey") +
    geom_text(aes(label = sprintf("%.2f", p.value)), 
              size = 2.5, color = "grey60", alpha = 0.5) +
    scale_fill_manual(
      values = sig_colors,
      drop = FALSE,
      name = "P-value",
      guide = guide_legend(byrow = TRUE, nrow = 1)
    ) +
    labs(
      x = "Covariate",
      y = "Trial – Country",
      subtitle = "Statistical Significance by Covariate"
    ) +
    ggthemes::theme_hc(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
      panel.grid = element_blank(), 
      legend.position = "top"
    ) + 
    facet_wrap(~type)
}

# Create and save heatmap
pvals <- get_heatmaps(heatmaps)

output_dir <- "../../results/figs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

png(paste0(output_dir, "/contextual_pvalue_heatmap.png"),  
    width = 3200, height = 1500, res = 300)
print(pvals)
dev.off()

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

cat("Saving results...\n")

# Create results directory
if (!dir.exists("../../results")) {
  dir.create("../../results", recursive = TRUE)
}

# Save analysis data with contextual variables
save(data_analysis, file = "../../results/contextual_analysis_data.rda")

# Save mixed effects results
save(mixed_effects_age, file = "../../results/mixed_effects_age_models.rda")

# Save ANOVA results
save(mixed_effects_age_anova, heatmaps, file = "../../results/contextual_anova_results.rda")

# Save correlation matrix
save(cor_matrix, full_matrix, file = "../../results/correlation_analysis.rda")

# Print summary
cat("\n=== CONTEXTUAL ANALYSIS SUMMARY ===\n")
cat("Covariates analyzed: ", length(covariate_names), "\n")
cat("Trial-countries: ", length(unique(paste(data_analysis$trial_phase, data_analysis$country))), "\n")
cat("Mixed effects models fitted: ", sum(!is.na(mixed_effects_age$model_full)), "\n")
cat("Significant covariate effects (p<0.05): ", 
    sum(heatmap_data1$p.value < 0.05, na.rm = TRUE), "/", nrow(heatmap_data1), "\n")
cat("Significant interactions (p<0.05): ", 
    sum(heatmap_data2$p.value < 0.05, na.rm = TRUE), "/", nrow(heatmap_data2), "\n")

cat("\nFiles saved:\n")
cat("  - ../../results/contextual_analysis_data.rda\n")
cat("  - ../../results/mixed_effects_age_models.rda\n")
cat("  - ../../results/contextual_anova_results.rda\n")
cat("  - ../../results/correlation_analysis.rda\n")
cat("  - ../../results/figs/contextual_pvalue_heatmap.png\n")
cat("  - ../../results/figs/correlation_matrix.png\n")
cat("  - Individual scatter plots in ../../results/figs/\n")

message("Contextual analysis completed successfully!")




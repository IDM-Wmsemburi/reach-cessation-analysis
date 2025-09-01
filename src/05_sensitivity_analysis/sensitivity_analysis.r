# ===============================================================================
# REACH Secondary Analysis: Enhanced Threshold Sensitivity Analysis 
# ===============================================================================
# Purpose: Derive age- and context-specific cessation thresholds using Stan & INLA
# Full posterior distributions + sensitivity analysis dropping locations
# Author: William Msemburi (enhanced with posterior densities)
# Date: August 2025
# Updated version with systematic exclusions only
# ===============================================================================

rm(list = ls())

library(tidyverse)
library(INLA)
library(ggridges)
library(posterior)
library(patchwork)
library(RColorBrewer)

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

# Load existing Stan models and data
results_dir <- "../../results"
load(paste0(results_dir, "/threshold_models.rda"))

# ===============================================================================
# INLA FITTING FUNCTIONS
# ===============================================================================

fit_inla_threshold <- function(df, mortality_col, exclude_trials = character()) {
  # Filter data
  if (length(exclude_trials) > 0) {
    df <- df %>% 
      mutate(trial_country = paste(trial_phase, country)) %>%
      filter(!(trial_country %in% exclude_trials))
  }
  
  # Fit INLA model
  dat <- df %>%
    mutate(
      log_mort = .data[[mortality_col]],
      log_off = log(person_time)
    )
  
  fit <- INLA::inla(
    deaths ~ 1 + log_mort + trt + log_mort:trt + f(group_id, model = "iid"),
    family = "poisson",
    data = dat,
    offset = log_off,
    control.compute = list(config = TRUE),
    num.threads = 2
  )
  
  return(list(
    fit = fit, 
    n_clusters = length(unique(dat$community_id)), 
    d_deaths = sum(dat$deaths), 
    py_tot = sum(dat$person_time),
    data = dat
  ))
}

extract_inla_coefs <- function(inla_fit, n_samp = 3000) {
  samps <- INLA::inla.posterior.sample(n = n_samp, result = inla_fit)
  
  # Get one sample to examine latent structure  
  s1 <- samps[[1]]$latent
  lat_names <- rownames(s1)
  
  # Helper to escape regex special characters
  esc <- function(x) gsub("([\\^$.|?*+(){}\\[\\]\\\\])","\\\\\\1", x)
  
  # Map each coefficient name to its latent name
  map1 <- function(nm) {
    # Try exact match first
    if (nm %in% lat_names) return(nm)
    # Try common INLA pattern with ":1"
    if (paste0(nm, ":1") %in% lat_names) return(paste0(nm, ":1"))
    # Try any ":k" pattern
    hit <- grep(paste0("^", esc(nm), ":[0-9]+$"), lat_names, value = TRUE)
    if (length(hit)) return(hit[1])
    stop("Couldn't find latent name for '", nm, "'.")
  }
  
  # Map coefficient names
  tryCatch({
    lat_map <- c(
      b0   = map1("(Intercept)"),
      b_m  = map1("log_mort"),
      b_t  = map1("trt"),
      b_int = map1("log_mort:trt")
    )
    
    # Extract coefficients from each sample
    coefs <- do.call(rbind, lapply(samps, function(s) {
      as.numeric(s$latent[unname(lat_map), 1])
    }))
    
    colnames(coefs) <- c("b0", "b_m", "b_t", "b_int")
    return(coefs)
    
  }, error = function(e) {
    cat("Error mapping INLA coefficients:", e$message, "\n")
    cat("Available latent names:", paste(head(lat_names, 10), collapse = ", "), "...\n")
    return(NULL)
  })
}

# ===============================================================================
# THRESHOLD EXTRACTION FUNCTIONS
# ===============================================================================

extract_analytical <- function(coefs) {
  if (is.null(coefs) || ncol(coefs) < 4) {
    return(c(median = NA, lower = NA, upper = NA))
  }
  
  # Only drop if denominator is 0 or result is not finite
  denom_valid <- coefs[,"b_int"] != 0
  thresh <- exp(-coefs[,"b_t"] / coefs[,"b_int"])
  valid_thresh <- thresh[denom_valid & is.finite(thresh)]
  
  if (length(valid_thresh) == 0) {
    return(c(median = NA, lower = NA, upper = NA))
  }
  
  c(
    median = median(valid_thresh) * 1000,
    lower = quantile(valid_thresh, 0.025) * 1000,
    upper = quantile(valid_thresh, 0.975) * 1000
  )
}

extract_crossing <- function(coefs, log_seq) {
  if (is.null(coefs) || ncol(coefs) < 4) {
    return(c(median = NA, lower = NA, upper = NA))
  }
  
  # Find crossings for each MCMC sample
  crossings <- apply(coefs, 1, function(row) {
    n_pred <- length(log_seq)
    eta_base <- row["b0"] + row["b_m"] * log_seq + log(365.25)
    eta_trt <- eta_base + row["b_t"] + row["b_int"] * log_seq
    
    # Find where treatment effect changes sign
    for (i in 1:(n_pred-1)) {
      d0 <- eta_trt[i] - eta_base[i]
      d1 <- eta_trt[i+1] - eta_base[i+1]
      if (abs(d0) < 1e-12) return(exp(log_seq[i]))
      if ((d0 < 0 && d1 > 0) || (d0 > 0 && d1 < 0)) {
        t_cross <- d0 / (d0 - d1)
        log_cross <- log_seq[i] + t_cross * (log_seq[i+1] - log_seq[i])
        return(exp(log_cross))
      }
    }
    return(NA)
  })
  
  valid_crossings <- crossings[is.finite(crossings)]
  
  if (length(valid_crossings) == 0) {
    return(c(median = NA, lower = NA, upper = NA))
  }
  
  c(
    median = median(valid_crossings) * 1000,
    lower = quantile(valid_crossings, 0.025) * 1000,
    upper = quantile(valid_crossings, 0.975) * 1000
  )
}

# ===============================================================================
# POSTERIOR DRAWS EXTRACTION
# ===============================================================================

extract_posterior_draws <- function() {
  
  cat("Extracting full posterior draws...\n")
  
  # Check what Stan objects are available
  variable_models_work <- tryCatch({
    test_u5 <- fit_u5mr$draws("beta_trt", format = "matrix")[1:5, 1]
    test_im <- fit_imr$draws("beta_trt", format = "matrix")[1:5, 1]
    TRUE
  }, error = function(e) {
    cat("Variable models not accessible:", e$message, "\n")
    FALSE
  })
  
  fixed_models_work <- tryCatch({
    if (exists("fit_u5mr_fixed") && exists("fit_imr_fixed")) {
      test_u5 <- fit_u5mr_fixed$draws("beta_trt", format = "matrix")[1:5, 1]
      test_im <- fit_imr_fixed$draws("beta_trt", format = "matrix")[1:5, 1]
      TRUE
    } else {
      FALSE
    }
  }, error = function(e) {
    cat("Fixed models not accessible:", e$message, "\n")
    FALSE
  })
  
  draws_list <- list()
  
  # Stan Variable models
  if (variable_models_work) {
    cat("  Extracting Stan Variable draws...\n")
    
    for (model in c("u5mr", "imr")) {
      fit <- if (model == "u5mr") fit_u5mr else fit_imr
      
      # Crossing-based draws
      crossing_draws <- posterior::as_draws_df(fit$draws(c("threshold", "has_cross"))) %>%
        filter(has_cross == 1, is.finite(threshold)) %>%
        mutate(
          threshold_1000 = threshold * 1000,
          scenario = "STAN Random Baseline",
          model = toupper(model),
          type = "crossing",
          draw = row_number()
        ) %>%
        select(scenario, model, type, draw, threshold_1000)
      
      # Analytical draws
      coef_draws <- posterior::as_draws_df(fit$draws(c("beta_trt", "beta_trt_mortality")))
      
      # Only drop if denominator is 0 or result is not finite
      analytical_draws <- coef_draws %>%
        filter(beta_trt_mortality != 0) %>%  # No division by zero
        mutate(
          threshold_raw = exp(-beta_trt / beta_trt_mortality),
          threshold_1000 = threshold_raw * 1000
        ) %>%
        filter(is.finite(threshold_1000)) %>%  # Only drop infinite/NaN
        mutate(
          scenario = "STAN Random Baseline",
          model = toupper(model),
          type = "analytical",
          draw = row_number()
        ) %>%
        select(scenario, model, type, draw, threshold_1000)
      
      draws_list <- append(draws_list, list(crossing_draws, analytical_draws))
      cat("    Extracted", nrow(crossing_draws), "crossing and", nrow(analytical_draws), "analytical draws for", toupper(model), "\n")
    }
  }
  
  # Stan Fixed models
  if (fixed_models_work) {
    cat("  Extracting Stan Fixed draws...\n")
    
    for (model in c("u5mr", "imr")) {
      fit <- if (model == "u5mr") fit_u5mr_fixed else fit_imr_fixed
      
      # Crossing-based draws
      crossing_draws <- posterior::as_draws_df(fit$draws(c("threshold", "has_cross"))) %>%
        filter(has_cross == 1, is.finite(threshold)) %>%
        mutate(
          threshold_1000 = threshold * 1000,
          scenario = "STAN Fixed Baseline",
          model = toupper(model),
          type = "crossing",
          draw = row_number()
        ) %>%
        select(scenario, model, type, draw, threshold_1000)
      
      # Analytical draws
      coef_draws <- posterior::as_draws_df(fit$draws(c("beta_trt", "beta_trt_mortality")))
      
      # Only drop if denominator is 0 or result is not finite
      analytical_draws <- coef_draws %>%
        filter(beta_trt_mortality != 0) %>%  # No division by zero
        mutate(
          threshold_raw = exp(-beta_trt / beta_trt_mortality),
          threshold_1000 = threshold_raw * 1000
        ) %>%
        filter(is.finite(threshold_1000)) %>%  # Only drop infinite/NaN
        mutate(
          scenario = "STAN Fixed Baseline",
          model = toupper(model),
          type = "analytical",
          draw = row_number()
        ) %>%
        select(scenario, model, type, draw, threshold_1000)
      
      draws_list <- append(draws_list, list(crossing_draws, analytical_draws))
      cat("    Extracted", nrow(crossing_draws), "crossing and", nrow(analytical_draws), "analytical draws for", toupper(model), "\n")
    }
  }
  
  # INLA models
  cat("  Extracting INLA draws...\n")
  
  inla_scenarios <- list(
    list(name = "INLA All", exclude = character()),
    list(name = "INLA -Niger MORDOR", exclude = "MORDOR I/II Niger"),
    list(name = "INLA -Niger AVENIR", exclude = "AVENIR Niger"),
    list(name = "INLA -Tanzania MORDOR", exclude = "MORDOR I Tanzania"),
    list(name = "INLA -Malawi MORDOR", exclude = "MORDOR I Malawi"),
    list(name = "INLA -Burkina CHAT", exclude = "CHAT Burkina Faso")
  )
  
  for (scenario in inla_scenarios) {
    cat("    Processing", scenario$name, "...\n")
    
    for (model in c("u5mr", "imr")) {
      tryCatch({
        mortality_col <- if (model == "u5mr") "post_u5mr_log" else "post_imr_log"
        log_seq <- if (model == "u5mr") {
          seq(log(0.035), log(0.235), length.out = 100)
        } else {
          seq(log(0.020), log(0.110), length.out = 100)
        }
        
        inla_result <- fit_inla_threshold(df, mortality_col, scenario$exclude)
        coefs <- extract_inla_coefs(inla_result$fit, n_samp = 3000)
        
        if (!is.null(coefs)) {
          # Analytical draws
          denom_valid <- coefs[,"b_int"] != 0
          analytical_raw <- exp(-coefs[,"b_t"] / coefs[,"b_int"])
          analytical_valid <- analytical_raw[denom_valid & is.finite(analytical_raw)]
          
          if (length(analytical_valid) > 0) {
            analytical_draws <- tibble(
              scenario = scenario$name,
              model = toupper(model),
              type = "analytical", 
              draw = seq_along(analytical_valid),
              threshold_1000 = analytical_valid * 1000
            )
            draws_list <- append(draws_list, list(analytical_draws))
          }
          
          # Crossing draws
          crossings <- apply(coefs, 1, function(row) {
            n_pred <- length(log_seq)
            eta_base <- row["b0"] + row["b_m"] * log_seq + log(365.25)
            eta_trt <- eta_base + row["b_t"] + row["b_int"] * log_seq
            
            for (i in 1:(n_pred-1)) {
              d0 <- eta_trt[i] - eta_base[i]
              d1 <- eta_trt[i+1] - eta_base[i+1]
              if (abs(d0) < 1e-12) return(exp(log_seq[i]))
              if ((d0 < 0 && d1 > 0) || (d0 > 0 && d1 < 0)) {
                t_cross <- d0 / (d0 - d1)
                log_cross <- log_seq[i] + t_cross * (log_seq[i+1] - log_seq[i])
                return(exp(log_cross))
              }
            }
            return(NA)
          })
          
          crossing_valid <- crossings[is.finite(crossings)]
          if (length(crossing_valid) > 0) {
            crossing_draws <- tibble(
              scenario = scenario$name,
              model = toupper(model),
              type = "crossing",
              draw = seq_along(crossing_valid),
              threshold_1000 = crossing_valid * 1000
            )
            draws_list <- append(draws_list, list(crossing_draws))
          }
          
          cat("      Extracted", length(analytical_valid), "analytical and", 
              length(crossing_valid), "crossing draws for", toupper(model), "\n")
        }
      }, error = function(e) {
        cat("      Error extracting draws for", scenario$name, model, ":", e$message, "\n")
      })
    }
  }
  
  return(bind_rows(draws_list))
}

# ===============================================================================
# SUMMARY STATISTICS EXTRACTION
# ===============================================================================

run_summary_scenarios <- function() {
  
  cat("Computing summary statistics...\n")
  
  # Check what Stan objects are available and working
  variable_models_work <- tryCatch({
    test_u5 <- fit_u5mr$draws("beta_trt", format = "matrix")[1:5, 1]
    test_im <- fit_imr$draws("beta_trt", format = "matrix")[1:5, 1]
    TRUE
  }, error = function(e) {
    cat("Variable models not accessible:", e$message, "\n")
    FALSE
  })
  
  fixed_models_work <- tryCatch({
    if (exists("fit_u5mr_fixed") && exists("fit_imr_fixed")) {
      test_u5 <- fit_u5mr_fixed$draws("beta_trt", format = "matrix")[1:5, 1]
      test_im <- fit_imr_fixed$draws("beta_trt", format = "matrix")[1:5, 1]
      TRUE
    } else {
      FALSE
    }
  }, error = function(e) {
    cat("Fixed models not accessible:", e$message, "\n")
    FALSE
  })
  
  # Define scenarios based on what works
  scenarios <- list()
  
  if (variable_models_work) {
    scenarios <- append(scenarios, list(list(name = "STAN Random Baseline", type = "stan", model = "variable")))
    cat("Including Stan Variable models\n")
  }
  
  if (fixed_models_work) {
    scenarios <- append(scenarios, list(list(name = "STAN Fixed Baseline", type = "stan", model = "fixed")))
    cat("Including Stan Fixed models\n")
  }
  
  # Add INLA scenarios
  inla_scenarios <- list(
    list(name = "INLA All", type = "inla", exclude = character()),
    list(name = "INLA -Niger MORDOR", type = "inla", exclude = "MORDOR I/II Niger"),
    list(name = "INLA -Niger AVENIR", type = "inla", exclude = "AVENIR Niger"),
    list(name = "INLA -Tanzania MORDOR", type = "inla", exclude = "MORDOR I Tanzania"),
    list(name = "INLA -Malawi MORDOR", type = "inla", exclude = "MORDOR I Malawi"),
    list(name = "INLA -Burkina CHAT", type = "inla", exclude = "CHAT Burkina Faso")
  )
  scenarios <- append(scenarios, inla_scenarios)
  
  # Prediction grids
  log_seq_u5 <- seq(log(0.035), log(0.235), length.out = 100)
  log_seq_im <- seq(log(0.020), log(0.110), length.out = 100)
  
  # Storage
  results <- list()
  
  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    cat("Processing:", scenario$name, "\n")
    
    if (scenario$type == "stan") {
      # Extract from existing Stan fits
      tryCatch({
        if (scenario$model == "variable") {
          fits <- list(u5mr = fit_u5mr, imr = fit_imr)
        } else {
          fits <- list(u5mr = fit_u5mr_fixed, imr = fit_imr_fixed)
        }
        
        # Extract Stan thresholds
        for (model in c("u5mr", "imr")) {
          tryCatch({
            fit <- fits[[model]]
            
            # Crossing thresholds
            crossing_draws <- posterior::as_draws_df(fit$draws(c("threshold", "has_cross")))
            crossing_draws$threshold[!is.finite(crossing_draws$threshold)] <- NA
            valid_crossings <- crossing_draws %>% 
              filter(has_cross == 1, is.finite(threshold)) %>% 
              pull(threshold)
            
            crossing_thresh <- if (length(valid_crossings) > 0) {
              c(median = median(valid_crossings) * 1000,
                lower = quantile(valid_crossings, 0.025) * 1000,
                upper = quantile(valid_crossings, 0.975) * 1000)
            } else {
              c(median = NA, lower = NA, upper = NA)
            }
            
            # Analytical thresholds
            coef_draws <- posterior::as_draws_df(fit$draws(c("beta_trt", "beta_trt_mortality")))
            
            # Only drop if denominator is 0 or result is not finite  
            analytical_valid <- coef_draws %>%
              filter(beta_trt_mortality != 0) %>%  # No division by zero
              mutate(threshold_raw = exp(-beta_trt / beta_trt_mortality)) %>%
              filter(is.finite(threshold_raw)) %>%  # Only drop infinite/NaN
              pull(threshold_raw)
            
            analytical_thresh <- if (length(analytical_valid) > 0) {
              c(median = median(analytical_valid) * 1000,
                lower = quantile(analytical_valid, 0.025) * 1000,
                upper = quantile(analytical_valid, 0.975) * 1000)
            } else {
              c(median = NA, lower = NA, upper = NA)
            }
            
            results[[length(results) + 1]] <- data.frame(
              scenario = scenario$name,
              model = toupper(model),
              type = "analytical",
              median = analytical_thresh[1],
              lower = analytical_thresh[2], 
              upper = analytical_thresh[3],
              n_clusters = length(unique(df$community_id)),
              total_deaths = sum(df$deaths), 
              total_person_time = sum(df$person_time)
            )
            
            results[[length(results) + 1]] <- data.frame(
              scenario = scenario$name,
              model = toupper(model),
              type = "crossing", 
              median = crossing_thresh[1],
              lower = crossing_thresh[2],
              upper = crossing_thresh[3],
              n_clusters = length(unique(df$community_id)),
              total_deaths = sum(df$deaths), 
              total_person_time = sum(df$person_time)
            )
            
          }, error = function(e) {
            cat("  Error processing", toupper(model), "model:", e$message, "\n")
          })
        }
        
      }, error = function(e) {
        cat("  Error with", scenario$name, ":", e$message, "\n")
      })
      
    } else {
      # Run INLA scenarios
      for (model in c("u5mr", "imr")) {
        tryCatch({
          mortality_col <- if (model == "u5mr") "post_u5mr_log" else "post_imr_log"
          log_seq <- if (model == "u5mr") log_seq_u5 else log_seq_im
          
          inla_result <- fit_inla_threshold(df, mortality_col, scenario$exclude)
          coefs <- extract_inla_coefs(inla_result$fit)
          
          analytical_thresh <- extract_analytical(coefs)
          crossing_thresh <- extract_crossing(coefs, log_seq)
          
          results[[length(results) + 1]] <- data.frame(
            scenario = scenario$name,
            model = toupper(model),
            type = "analytical",
            median = analytical_thresh[1],
            lower = analytical_thresh[2],
            upper = analytical_thresh[3],
            n_clusters = inla_result$n_clusters,
            total_deaths = inla_result$d_deaths, 
            total_person_time = inla_result$py_tot
          )
          
          results[[length(results) + 1]] <- data.frame(
            scenario = scenario$name,
            model = toupper(model),
            type = "crossing",
            median = crossing_thresh[1],
            lower = crossing_thresh[2],
            upper = crossing_thresh[3],
            n_clusters = inla_result$n_clusters,
            total_deaths = inla_result$d_deaths, 
            total_person_time = inla_result$py_tot
          )
          
        }, error = function(e) {
          cat("  Error processing", toupper(model), "INLA model:", e$message, "\n")
        })
      }
    }
  }
  
  return(bind_rows(results))
}

# ===============================================================================
# ENHANCED VISUALIZATION
# ===============================================================================

create_density_plots <- function(posterior_draws, summary_results) {
  
  create_density_plot <- function(threshold_type) {
    # Filter posterior draws - only basic filtering for visualization
    plot_draws <- posterior_draws %>%
      filter(type == threshold_type) %>%
      # Only filter for clear display purposes, not arbitrary exclusions
      filter(threshold_1000 >= 0, threshold_1000 <= 100)  # Sensible display range
    
    # Get summary stats for overlay
    plot_summaries <- summary_results %>%
      filter(type == threshold_type, !is.na(median))
    
    # Order scenarios by U5MR median threshold (highest to lowest)
    u5mr_medians <- plot_summaries %>%
      filter(model == "U5MR") %>%
      arrange(desc(median)) %>%
      pull(scenario)
    
    # Filter to available scenarios and apply U5MR-based ordering
    available_scenarios <- intersect(u5mr_medians, unique(plot_draws$scenario))
    
    plot_draws <- plot_draws %>%
      filter(scenario %in% available_scenarios) %>%
      mutate(scenario = factor(scenario, levels = available_scenarios))
    
    plot_summaries <- plot_summaries %>%
      filter(scenario %in% available_scenarios) %>%
      mutate(scenario = factor(scenario, levels = available_scenarios))
    
    # Create median labels for display
    median_labels <- plot_summaries %>%
      mutate(
        label = paste0(round(median, 1), " (", round(lower), "-", round(upper), ")"),
        x_pos = median,
        y_pos = as.numeric(scenario) + ifelse(model == "U5MR", -0.22, -0.27)
      )
    
    # Create the plot - both models on same row
    p <- ggplot(plot_draws, aes(x = threshold_1000, y = scenario, fill = model)) +
      
      # Posterior density ridges - both models on same y-position
      geom_density_ridges(
        alpha = 0.6, 
        scale = 0.8, 
        rel_min_height = 0.01,
        bandwidth = 3
      ) +
      
      # Overlay summary statistics with slight vertical offset to avoid overlap
      geom_point(
        data = plot_summaries, 
        aes(x = median, y = scenario, color = model), 
        size = 2.5, 
        alpha = 1,
        position = position_nudge(y = ifelse(plot_summaries$model == "U5MR", -0.1, -0.15)),
        inherit.aes = FALSE
      ) +
      
      geom_errorbarh(
        data = plot_summaries,
        aes(xmin = lower, xmax = upper, y = scenario, color = model),
        height = 0.05, 
        alpha = 0.8,
        linewidth = 0.8,
        position = position_nudge(y = ifelse(plot_summaries$model == "U5MR", -0.1, -0.15)),
        inherit.aes = FALSE
      ) +
      
      # Add median value labels
      geom_text(
        data = median_labels,
        aes(x = x_pos, y = y_pos, label = label, color = model),
        size = 3,
        fontface = "bold",
        hjust = 0.5,
        inherit.aes = FALSE
      ) +
      
      scale_fill_manual(values = c("U5MR" = "steelblue", "IMR" = "darkred"),
                        name = "Mortality\nMeasure") +
      scale_color_manual(values = c("U5MR" = "navy", "IMR" = "firebrick"),
                         name = "Mortality\nMeasure") +
      
      # Set sensible x-axis limits
      xlim(0, 100) +
      
      labs(
        title = paste(str_to_title(threshold_type), "Threshold Distributions"),
        subtitle = "Density ridges: full posterior distributions | Points + bars: medians with 95% CIs",
        x = "Threshold (per 1,000 live births)",
        y = "Analysis Scenario"
      ) +
      
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40")
      ) +
      
      guides(
        fill = guide_legend(override.aes = list(alpha = 0.7)),
        color = guide_legend(override.aes = list(size = 2))
      )
    
    return(p)
  }
  
  list(
    analytical = create_density_plot("analytical"),
    crossing = create_density_plot("crossing")
  )
}

create_summary_tables <- function(results_df) {
  create_table <- function(threshold_type) {
    results_df %>%
      filter(type == threshold_type, !is.na(median)) %>%
      select(scenario, model, median, lower, upper, n_clusters, total_deaths) %>%
      mutate(
        ci_95 = paste0(round(median, 1), " (", round(lower, 1), "-", round(upper, 1), ")")
      ) %>%
      select(scenario, model, ci_95, n_clusters, total_deaths) %>%
      pivot_wider(
        names_from = model, 
        values_from = c(ci_95, n_clusters, total_deaths),
        names_sep = "_"
      ) %>%
      select(
        scenario, 
        n_clusters_U5MR, total_deaths_U5MR, ci_95_U5MR,
        n_clusters_IMR, total_deaths_IMR, ci_95_IMR
      )
  }
  
  list(
    analytical = create_table("analytical"),
    crossing = create_table("crossing")
  )
}

# ===============================================================================
# ROBUSTNESS ASSESSMENT
# ===============================================================================

assess_robustness <- function(results_df, posterior_draws) {
  # Coefficient of variation across scenarios
  cv_analysis <- results_df %>%
    filter(!is.na(median)) %>%
    group_by(model, type) %>%
    summarise(
      mean_threshold = mean(median, na.rm = TRUE),
      sd_threshold = sd(median, na.rm = TRUE),
      cv = sd_threshold / mean_threshold,
      n_scenarios = n(),
      range_width = max(median, na.rm = TRUE) - min(median, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      robustness = case_when(
        cv < 0.1 ~ "Very Robust",
        cv < 0.2 ~ "Robust", 
        cv < 0.3 ~ "Moderately Robust",
        TRUE ~ "Not Robust"
      )
    )
  
  # Posterior uncertainty analysis
  uncertainty_analysis <- posterior_draws %>%
    group_by(scenario, model, type) %>%
    summarise(
      n_draws = n(),
      mean_threshold = mean(threshold_1000, na.rm = TRUE),
      sd_threshold = sd(threshold_1000, na.rm = TRUE),
      cv_posterior = sd_threshold / mean_threshold,
      .groups = "drop"
    ) %>%
    group_by(model, type) %>%
    summarise(
      avg_posterior_cv = mean(cv_posterior, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    scenario_robustness = cv_analysis,
    posterior_uncertainty = uncertainty_analysis
  )
}

# ===============================================================================
# MAIN ANALYSIS EXECUTION
# ===============================================================================

run_enhanced_analysis <- function() {
  
  cat("=== ENHANCED THRESHOLD SENSITIVITY ANALYSIS ===\n")
  cat("Starting comprehensive analysis...\n\n")
  
  # Extract full posterior distributions
  posterior_draws <- extract_posterior_draws()
  
  # Get summary statistics
  summary_results <- run_summary_scenarios()
  
  # Create enhanced visualizations
  density_plots <- create_density_plots(posterior_draws, summary_results)
  
  # Create summary tables
  summary_tables <- create_summary_tables(summary_results)
  
  # Assess robustness
  robustness <- assess_robustness(summary_results, posterior_draws)
  
  # Print key results
  cat("\n=== SUMMARY STATISTICS ===\n")
  print(summary_results %>% filter(!is.na(median)) %>% head(10))
  
  cat("\n=== ROBUSTNESS ASSESSMENT ===\n")
  print(robustness$scenario_robustness)
  
  # Save all results
  save(
    posterior_draws, 
    summary_results, 
    density_plots, 
    summary_tables, 
    robustness,
    file = paste0(results_dir, "/enhanced_threshold_analysis.rda")
  )
  
  # Export tables
  readr::write_csv(summary_tables$analytical, 
                   paste0(results_dir, "/analytical_thresholds_enhanced.csv"))
  readr::write_csv(summary_tables$crossing, 
                   paste0(results_dir, "/crossing_thresholds_enhanced.csv"))
  readr::write_csv(summary_results,
                   paste0(results_dir, "/all_threshold_results.csv"))
  readr::write_csv(robustness$scenario_robustness,
                   paste0(results_dir, "/threshold_robustness.csv"))
  
  # Save enhanced plots
  ggsave(paste0(results_dir, "/figs/analytical_threshold_densities.pdf"), 
         density_plots$analytical, width = 8, height = 7)
  ggsave(paste0(results_dir, "/figs/crossing_threshold_densities.pdf"), 
         density_plots$crossing, width = 8, height = 7)
  
  # Combined plot
  combined_plot <- density_plots$analytical / density_plots$crossing +
    plot_annotation(
      title = "Comprehensive Threshold Sensitivity Analysis",
      subtitle = "Posterior distributions across Stan and INLA models with various exclusion scenarios",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  ggsave(paste0(results_dir, "/figs/combined_threshold_analysis.pdf"), 
         combined_plot, width = 16, height = 12)
  
  cat("\n=== FILES CREATED ===\n")
  cat("Data files:\n")
  cat("  - enhanced_threshold_analysis.rda (all objects)\n")
  cat("  - analytical_thresholds_enhanced.csv\n")
  cat("  - crossing_thresholds_enhanced.csv\n") 
  cat("  - all_threshold_results.csv\n")
  cat("  - threshold_robustness.csv\n")
  cat("\nPlot files:\n")
  cat("  - figs/analytical_threshold_densities.pdf\n")
  cat("  - figs/crossing_threshold_densities.pdf\n")
  cat("  - figs/combined_threshold_analysis.pdf\n")
  
  return(list(
    posterior_draws = posterior_draws,
    summary_results = summary_results,
    plots = density_plots,
    tables = summary_tables,
    robustness = robustness
  ))
}

# ===============================================================================
# RUN THE ANALYSIS
# ===============================================================================

# Execute the enhanced analysis
enhanced_results <- run_enhanced_analysis()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Enhanced threshold sensitivity analysis completed successfully!\n")

# Print final summary
cat("\nFinal recommendations based on robustness:\n")
robust_models <- enhanced_results$robustness$scenario_robustness %>%
  filter(robustness %in% c("Robust", "Very Robust")) %>%
  arrange(cv)

if (nrow(robust_models) > 0) {
  cat("Most robust threshold estimates:\n")
  print(robust_models)
} else {
  cat("Consider using median values across scenarios due to moderate robustness\n")
}
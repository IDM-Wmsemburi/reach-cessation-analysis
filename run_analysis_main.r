# ===============================================================================
# REACH Secondary Analysis: Main Analysis Runner
# ===============================================================================
# Purpose: Coordinate execution of the complete REACH cessation analysis pipeline
# Author: William Msemburi
# Date: August 2025
# ===============================================================================

# Clear environment and set options
rm(list = ls())
options(warn = -1)
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

# ===============================================================================
# COMMAND LINE ARGUMENT PARSING
# ===============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract stage argument if provided
stage_arg <- NULL
if (length(args) > 0) {
  stage_matches <- grep("--stage=", args, value = TRUE)
  if (length(stage_matches) > 0) {
    stage_arg <- sub("--stage=", "", stage_matches[1])
  }
}

# Define available stages
available_stages <- c("prep", "baseline", "context", "thresholds", "all")

# Default to running all stages
run_stage <- if (is.null(stage_arg) || !stage_arg %in% available_stages) "all" else stage_arg

# ===============================================================================
# SETUP AND CONFIGURATION
# ===============================================================================

cat("=== REACH SECONDARY ANALYSIS PIPELINE ===\n")
cat("Starting analysis pipeline...\n")
cat("Stage:", run_stage, "\n")
cat("Start time:", format(Sys.time()), "\n\n")

# Record start time
pipeline_start <- Sys.time()

project_root <- getwd()

cat("Project root:", project_root, "\n")

# Load core utilities
source("src/utils/reach_analysis_utils.R")

# Create required directories
required_dirs <- c("results", "results/figs", "models/stan")
for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat("Created directory:", dir, "\n")
  }
}

# ===============================================================================
# PIPELINE STAGE DEFINITIONS
# ===============================================================================

#' Run data preparation stage
#' @return TRUE if successful, FALSE otherwise
run_prep_stage <- function() {
  cat("\n=== STAGE 1: DATA PREPARATION ===\n")
  
  stage_start <- Sys.time()
  
  tryCatch({
    cat("1.1 Extracting DHS hazards...\n")
    source("src/01_data_preparation/prepare_dhs_hazards.R")
    
    cat("1.2 Deriving trial priors...\n")
    source("src/01_data_preparation/derive_trial_priors.R")
    
    stage_end <- Sys.time()
    elapsed <- round(as.numeric(stage_end - stage_start, units = "mins"), 2)
    cat("Data preparation completed in", elapsed, "minutes\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR in data preparation stage:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

#' Run baseline mortality estimation stage
#' @return TRUE if successful, FALSE otherwise
run_baseline_stage <- function() {
  cat("\n=== STAGE 2: BASELINE MORTALITY ESTIMATION ===\n")
  
  stage_start <- Sys.time()
  
  tryCatch({
    cat("2.1 Running Bayesian HSGP spatial model...\n")
    source("src/02_baseline_mortality/bayesian_hsgp_model.R")
    
    stage_end <- Sys.time()
    elapsed <- round(as.numeric(stage_end - stage_start, units = "mins"), 2)
    cat("Baseline mortality estimation completed in", elapsed, "minutes\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR in baseline mortality stage:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

#' Run contextual analysis stage
#' @return TRUE if successful, FALSE otherwise
run_context_stage <- function() {
  cat("\n=== STAGE 3: CONTEXTUAL ANALYSIS ===\n")
  
  stage_start <- Sys.time()
  
  tryCatch({
    cat("3.1 Analyzing contextual heterogeneity...\n")
    source("src/03_contextual_analysis/contextual_heterogeneity.R")
    
    stage_end <- Sys.time()
    elapsed <- round(as.numeric(stage_end - stage_start, units = "mins"), 2)
    cat("Contextual analysis completed in", elapsed, "minutes\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR in contextual analysis stage:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

#' Run threshold estimation stage  
#' @return TRUE if successful, FALSE otherwise
run_thresholds_stage <- function() {
  cat("\n=== STAGE 4: THRESHOLD ESTIMATION ===\n")
  
  stage_start <- Sys.time()
  
  tryCatch({
    cat("4.1 Estimating cessation thresholds...\n")
    source("src/04_threshold_estimation/threshold_analysis.R")
    
    stage_end <- Sys.time()
    elapsed <- round(as.numeric(stage_end - stage_start, units = "mins"), 2)
    cat("Threshold estimation completed in", elapsed, "minutes\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR in threshold estimation stage:", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# ===============================================================================
# DEPENDENCY CHECKING
# ===============================================================================

#' Check if required input files exist for a given stage
#' @param stage Stage name to check
#' @return TRUE if dependencies exist, FALSE otherwise
check_stage_dependencies <- function(stage) {
  
  dependencies <- list(
    prep = c(
      # DHS data files should exist
      "data/DHS" # directory should exist
    ),
    baseline = c(
      "results/dhs_hazard_data.rda",
      "results/dhs_h0_pooled_model.rds", 
      "results/trial_hazard_priors.rds",
      "data/Clean/data_crt_geospatial.rda"
    ),
    context = c(
      "data/Clean/data_crt_geospatial.rda",
      "results/baseline_mortality_estimates.rda"
    ),
    thresholds = c(
      "data/Clean/data_crt_geospatial.rda", 
      "results/baseline_mortality_estimates.rda"
    )
  )
  
  if (!stage %in% names(dependencies)) {
    return(TRUE) # No dependencies defined
  }
  
  required_files <- dependencies[[stage]]
  missing_files <- c()
  
  for (file in required_files) {
    if (!file.exists(file) && !dir.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    cat("Missing dependencies for stage '", stage, "':\n")
    for (file in missing_files) {
      cat("  -", file, "\n")
    }
    return(FALSE)
  }
  
  return(TRUE)
}

# ===============================================================================
# MAIN PIPELINE EXECUTION
# ===============================================================================

#' Execute the analysis pipeline
#' @param stages Character vector of stages to run
#' @return Named logical vector indicating success/failure of each stage
run_pipeline <- function(stages = "all") {
  
  if ("all" %in% stages) {
    stages <- c("prep", "baseline", "context", "thresholds")
  }
  
  results <- c()
  
  for (stage in stages) {
    
    # Check dependencies
    if (!check_stage_dependencies(stage)) {
      cat("Skipping stage '", stage, "' due to missing dependencies.\n")
      results[stage] <- FALSE
      next
    }
    
    # Run the stage
    success <- switch(stage,
      "prep" = run_prep_stage(),
      "baseline" = run_baseline_stage(), 
      "context" = run_context_stage(),
      "thresholds" = run_thresholds_stage(),
      {
        cat("Unknown stage:", stage, "\n")
        FALSE
      }
    )
    
    results[stage] <- success
    
    if (!success) {
      cat("Stage '", stage, "' failed. Stopping pipeline.\n")
      break
    }
  }
  
  return(results)
}

# ===============================================================================
# RESULTS SUMMARY FUNCTION
# ===============================================================================

#' Generate a summary of pipeline results
#' @return Invisible TRUE
generate_pipeline_summary <- function() {
  cat("\n=== PIPELINE SUMMARY ===\n")
  
  # Check which result files exist
  result_files <- c(
    "DHS hazards" = "results/dhs_hazard_data.rda",
    "DHS H0 model" = "results/dhs_h0_pooled_model.rds",
    "Trial priors" = "results/trial_hazard_priors.rds", 
    "Baseline mortality" = "results/baseline_mortality_estimates.rda",
    "Contextual data" = "results/contextual_analysis_data.rda",
    "Threshold results" = "results/threshold_analysis_results.rda"
  )
  
  cat("Generated result files:\n")
  for (name in names(result_files)) {
    file <- result_files[[name]]
    status <- if (file.exists(file)) "✓" else "✗"
    cat("  ", status, name, "\n")
  }
  
  # Count figures
  figs_dir <- "results/figs"
  if (dir.exists(figs_dir)) {
    fig_count <- length(list.files(figs_dir, pattern = "\\.(png|pdf)$"))
    cat("\nGenerated figures:", fig_count, "\n")
  }
  
  # Check for key outputs
  key_outputs <- c(
    "results/threshold_results.csv",
    "results/figs/threshold_analysis_plots.png",
    "results/figs/contextual_pvalue_heatmap.png"
  )
  
  cat("\nKey outputs:\n")
  for (output in key_outputs) {
    status <- if (file.exists(output)) "✓" else "✗"
    cat("  ", status, basename(output), "\n")
  }
  
  invisible(TRUE)
}

# ===============================================================================
# EXECUTE PIPELINE
# ===============================================================================

# Main execution
tryCatch({
  
  # Load required packages for utilities
  load_reach_packages()
  
  # Print analysis summary function
  print_reach_utils_summary()
  
  # Run the requested stage(s)
  results <- run_pipeline(run_stage)
  
  # Calculate total elapsed time  
  pipeline_end <- Sys.time()
  total_elapsed <- round(as.numeric(pipeline_end - pipeline_start, units = "mins"), 2)
  
  # Print results summary
  cat("\n=== EXECUTION RESULTS ===\n")
  for (stage in names(results)) {
    status <- if (results[[stage]]) "SUCCESS" else "FAILED"
    cat(stage, ":", status, "\n")
  }
  
  cat("\nTotal execution time:", total_elapsed, "minutes\n")
  
  # Generate summary of outputs
  generate_pipeline_summary()
  
  # Final message
  if (all(results)) {
    cat("\n🎉 Analysis pipeline completed successfully!\n")
  } else {
    cat("\n⚠️ Analysis pipeline completed with errors.\n")
  }
  
  cat("\nEnd time:", format(Sys.time()), "\n")
  
}, error = function(e) {
  cat("\n❌ PIPELINE ERROR:", conditionMessage(e), "\n")
  cat("Please check the error messages above and ensure all dependencies are installed.\n")
})

# ===============================================================================
# USAGE EXAMPLES AND HELP
# ===============================================================================

if (interactive()) {
  cat("\n=== USAGE EXAMPLES ===\n")
  cat("# Run complete pipeline:\n")
  cat("source('run_analysis.R')\n\n")
  
  cat("# Run specific stages:\n")
  cat("source('run_analysis.R'); run_pipeline('prep')\n")
  cat("source('run_analysis.R'); run_pipeline(c('baseline', 'context'))\n\n")
  
  cat("# Command line usage:\n")
  cat("Rscript run_analysis.R --stage=all\n")
  cat("Rscript run_analysis.R --stage=thresholds\n\n")
  
  cat("# Check dependencies:\n")
  cat("check_stage_dependencies('baseline')\n\n")
}
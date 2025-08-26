# ===============================================================================
# REACH Secondary Analysis: CmdStan Configuration
# ===============================================================================
# Purpose: Configure CmdStan for consistent usage across analysis scripts
# Author: William Msemburi
# Date: August 2025
# ===============================================================================

# Load required packages
if (!require("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr package is required but not installed. Install with: install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))")
}

#' Configure CmdStan for the current system
#' @param force_recheck Force rechecking CmdStan installation
#' @return TRUE if successfully configured, FALSE otherwise
configure_cmdstan <- function(force_recheck = FALSE) {
  
  cat("Configuring CmdStan...\n")
  
  # Check system type for WSL configuration
  if (Sys.info()["sysname"] == "Windows") {
    cat("Detected Windows system\n")
    
    # Check if WSL is available and preferred
    if (Sys.getenv("CMDSTANR_USE_WSL") == "TRUE" || 
        !nzchar(Sys.which("make")) || 
        !dir.exists("C:/tools/msys64")) {
      
      cat("Configuring for WSL usage\n")
      Sys.setenv(CMDSTANR_USE_WSL = "TRUE")
      
    } else {
      cat("Configuring for native Windows usage\n")
    }
  }
  
  # Check if CmdStan is installed and working
  if (force_recheck || !cmdstanr::cmdstan_version(error_on_NA = FALSE) %in% c(NA, "")) {
    
    tryCatch({
      cmdstan_version <- cmdstanr::cmdstan_version()
      cat("Found CmdStan version:", cmdstan_version, "\n")
      
      # Test compilation with a simple model
      test_model_code <- "
        data {
          int<lower=0> N;
        }
        parameters {
          real theta;
        }
        model {
          theta ~ normal(0, 1);
        }
      "
      
      temp_file <- tempfile(fileext = ".stan")
      writeLines(test_model_code, temp_file)
      
      test_model <- cmdstanr::cmdstan_model(temp_file, compile = TRUE)
      cat("CmdStan compilation test successful\n")
      
      # Clean up
      unlink(temp_file)
      
      return(TRUE)
      
    }, error = function(e) {
      cat("CmdStan test failed:", conditionMessage(e), "\n")
      return(FALSE)
    })
    
  } else {
    
    cat("CmdStan not found or needs installation\n")
    
    # Attempt to install CmdStan
    tryCatch({
      cat("Attempting to install CmdStan...\n")
      cmdstanr::install_cmdstan()
      
      # Verify installation
      cmdstan_version <- cmdstanr::cmdstan_version()
      cat("CmdStan installed successfully, version:", cmdstan_version, "\n")
      return(TRUE)
      
    }, error = function(e) {
      cat("CmdStan installation failed:", conditionMessage(e), "\n")
      cat("Please install CmdStan manually using cmdstanr::install_cmdstan()\n")
      return(FALSE)
    })
  }
}

#' Get appropriate model directory for Stan files
#' @param create_if_missing Create directory if it doesn't exist
#' @return Path to model directory
get_model_directory <- function(create_if_missing = TRUE) {
  
  # Determine project root
  if (file.exists("run_analysis.R")) {
    project_root <- getwd()
  } else if (file.exists("../run_analysis.R")) {
    project_root <- normalizePath("..")
  } else if (file.exists("../../run_analysis.R")) {
    project_root <- normalizePath("../..")
  } else {
    project_root <- getwd()
    warning("Could not determine project root, using current directory")
  }
  
  model_dir <- file.path(project_root, "models", "stan")
  
  if (create_if_missing && !dir.exists(model_dir)) {
    dir.create(model_dir, recursive = TRUE)
    cat("Created model directory:", model_dir, "\n")
  }
  
  return(model_dir)
}

#' Set up parallel processing options for Stan
#' @param chains Number of chains to run
#' @param cores Number of CPU cores to use (NULL for auto-detect)
#' @param threads_per_chain Threads per chain for within-chain parallelization
#' @return List of parallel processing settings
setup_parallel_stan <- function(chains = 4, cores = NULL, threads_per_chain = 1) {
  
  # Auto-detect cores if not specified
  if (is.null(cores)) {
    available_cores <- parallel::detectCores()
    cores <- min(chains, available_cores - 1, na.rm = TRUE)
    cores <- max(cores, 1)
  }
  
  cat("Parallel Stan setup:\n")
  cat("  Chains:", chains, "\n")
  cat("  Parallel chains:", cores, "\n")  
  cat("  Threads per chain:", threads_per_chain, "\n")
  
  # Set up parallel options
  options(mc.cores = cores)
  
  return(list(
    chains = chains,
    parallel_chains = cores,
    threads_per_chain = threads_per_chain
  ))
}

#' Get recommended Stan sampling settings
#' @param model_complexity Character: "simple", "medium", "complex"
#' @return Named list of Stan sampling parameters
get_stan_sampling_settings <- function(model_complexity = "medium") {
  
  settings <- switch(model_complexity,
    "simple" = list(
      iter_warmup = 500,
      iter_sampling = 1000,
      adapt_delta = 0.8,
      max_treedepth = 10
    ),
    "medium" = list(
      iter_warmup = 1000,
      iter_sampling = 1500,
      adapt_delta = 0.9,
      max_treedepth = 12
    ),
    "complex" = list(
      iter_warmup = 2000,
      iter_sampling = 2000,
      adapt_delta = 0.95,
      max_treedepth = 15
    ),
    {
      warning("Unknown model complexity, using medium settings")
      list(
        iter_warmup = 1000,
        iter_sampling = 1500,
        adapt_delta = 0.9,
        max_treedepth = 12
      )
    }
  )
  
  return(settings)
}

#' Print CmdStan system information
#' @return Invisible TRUE
print_cmdstan_info <- function() {
  
  cat("=== CMDSTAN SYSTEM INFO ===\n")
  cat("Operating System:", Sys.info()["sysname"], "\n")
  cat("R Version:", R.version.string, "\n")
  
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    cat("cmdstanr version:", as.character(utils::packageVersion("cmdstanr")), "\n")
    
    tryCatch({
      cmdstan_version <- cmdstanr::cmdstan_version()
      cat("CmdStan version:", cmdstan_version, "\n")
      cat("CmdStan path:", cmdstanr::cmdstan_path(), "\n")
    }, error = function(e) {
      cat("CmdStan: Not found or not working\n")
    })
  } else {
    cat("cmdstanr: Not installed\n")
  }
  
  cat("Available cores:", parallel::detectCores(), "\n")
  
  # WSL info for Windows
  if (Sys.info()["sysname"] == "Windows") {
    cat("WSL usage:", Sys.getenv("CMDSTANR_USE_WSL", "FALSE"), "\n")
  }
  
  cat("============================\n")
  
  invisible(TRUE)
}

#' Check Stan model syntax without compilation
#' @param stan_code Character string containing Stan code
#' @param model_name Optional model name for temporary file
#' @return TRUE if syntax is valid, FALSE otherwise
check_stan_syntax <- function(stan_code, model_name = "temp_model") {
  
  tryCatch({
    temp_file <- tempfile(pattern = model_name, fileext = ".stan")
    writeLines(stan_code, temp_file)
    
    # Try to parse the model (this checks syntax without compiling)
    model <- cmdstanr::cmdstan_model(temp_file, compile = FALSE)
    
    # Clean up
    unlink(temp_file)
    
    cat("Stan syntax check passed for", model_name, "\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("Stan syntax error in", model_name, ":", conditionMessage(e), "\n")
    
    # Clean up
    if (exists("temp_file")) unlink(temp_file)
    
    return(FALSE)
  })
}

# ===============================================================================
# AUTO-CONFIGURATION
# ===============================================================================

# Automatically configure CmdStan when this file is sourced
if (!exists("CMDSTAN_CONFIGURED") || !CMDSTAN_CONFIGURED) {
  
  cat("Running automatic CmdStan configuration...\n")
  
  config_success <- configure_cmdstan(force_recheck = FALSE)
  
  if (config_success) {
    cat("CmdStan configured successfully\n")
    CMDSTAN_CONFIGURED <- TRUE
  } else {
    cat("CmdStan configuration failed - some analyses may not work\n")
    CMDSTAN_CONFIGURED <- FALSE
  }
  
  # Print system info if in interactive mode
  if (interactive()) {
    print_cmdstan_info()
  }
}
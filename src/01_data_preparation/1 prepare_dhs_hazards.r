# ===============================================================================
# REACH Secondary Analysis: DHS Hazard Extraction
# ===============================================================================
# Purpose: Extract age-specific mortality hazards from DHS birth histories
# Author: William Msemburi 
# Date: August 2025
# ===============================================================================

# Clear environment and set options
rm(list = ls())
options(warn = -1, gsubfn.engine = "R")

# Load required libraries
libraries <- c("SUMMER", "dplyr", "tidyr", 
               "INLA", "survey", "ggplot2",  
               "haven", "labelled", "data.table",
               "mgcv", "ggthemes")
invisible(lapply(libraries, library, character.only = TRUE))

# Set working directory to script location
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

# Load utility functions
source("../utils/reach_analysis_utils.R")

# ===============================================================================
# DHS HAZARD EXTRACTION FUNCTION
# ===============================================================================

#' Extract direct hazards from DHS data with survey weighting
#' @param mod.data DHS survey data processed by SUMMER::getBirths
#' @param regionVar Region variable name (default: "region")
#' @param timeVar Time variable name (default: "years")
#' @param clusterVar Cluster variable formula (default: "~cluster")
#' @param ageVar Age variable name (default: "age")
#' @param Ntrials Number of trials variable (default: "total")
#' @param weightsVar Weights variable (default: "v005")
#' @param CI Confidence interval level (default: 0.95)
#' @return data.table with hazard estimates by region, time, and age group
getDirectHazards <- function(mod.data, 
                            regionVar = "region",
                            timeVar   = "years",
                            clusterVar =  "~cluster",
                            ageVar = "age", 
                            Ntrials = "total",
                            weightsVar = "v005",
                            CI = 0.95) {
  
  require(data.table)
  require(survey)
  require(MASS)
  
  isids <- unique(mod.data$isid)
  results.all <- list()
  
  # Process each survey separately
  for (isi in isids) {
    print(paste("Processing DHS", isi))
    births <- mod.data |> filter(isid == isi)
    years <- sort(unique(births$years))
    
    # Set up variable names
    births$died <- births$Y
    births$region0 <- births[, regionVar]
    births$weights0 <- births[, weightsVar]
    births$time0 <- births[, timeVar]
    births$age0 <- births[, ageVar]
    
    # Age group recoding with interval lengths
    if (sum(c(0, 1, 12, 24, 36, 48) %in% births$age0) == 6) {
      births$age0[births$age0 == 0] <- "0"
      births$age0[births$age0 == 1] <- "1-11"
      births$age0[births$age0 == 12] <- "12-23"
      births$age0[births$age0 == 24] <- "24-35"
      births$age0[births$age0 == 36] <- "36-47"
      births$age0[births$age0 == 48] <- "48-59"
      ns <- c(1, 11, 12, 12, 12, 12)
    }
    
    # Calculate interval lengths for each age group
    labels <- as.character(unique(births$age0))
    ns <- rep(1, length(labels))
    for (i in 1:length(labels)) {
      if (labels[i] == "0") {
        ns[i] <- 1
        next
      }
      tmp <- as.numeric(strsplit(as.character(labels[i]), "-")[[1]])
      ns[i] <- tmp[2] - tmp[1] + 1
    }
    
    births$strata <- factor(births$strata)
    births$n <- births[, Ntrials]
    births$died <- births$died/births$n
    
    # Set up survey design
    options(survey.lonely.psu = "adjust")
    my.svydesign <- survey::svydesign(
      ids = stats::formula(clusterVar), 
      strata = ~strata, 
      nest = TRUE, 
      weights = ~weights0, 
      data = births
    )
    
    # Define regions
    regions_list <- as.character(sort(names(table(births$region0))[as.vector(table(births$region0) != 0)]))
    regions_num <- 1:length(regions_list)
    regions_list <- c("All", regions_list)
    regions_num <- c(0, regions_num)
    
    # Function to fit hazard model for each region-time combination
    region.time.HT.withNA <- function(which.area, which.time) {
      if (which.area == "All") {
        tmp <- subset(my.svydesign, (time0 == which.time))
      } else {
        tmp <- subset(my.svydesign, (time0 == which.time & region0 == which.area))
      }
      
      if (dim(tmp)[1] == 0) {
        return(rep(NA, 5))
      } else if (sum(tmp$variables$died) == 0) {
        message(paste0(which.area, " ", which.time, " has no deaths, set to NA"))
        return(rep(NA, 5))
      } else if (length(unique(tmp$variables$age0)) > 1) {
        
        # Fit GLM with survey design
        if (is.null(Ntrials)) {
          glm.ob <- survey::svyglm(died ~ (-1) + factor(age0), 
                                   design = tmp, family = stats::quasibinomial, 
                                   maxit = 50)
        } else {
          glm.ob <- survey::svyglm(died ~ (-1) + factor(age0), 
                                   design = tmp, family = stats::quasibinomial, 
                                   maxit = 50, weights = tmp$variables$n)
        }
        
        # Handle missing age groups
        if (dim(summary(glm.ob)$coef)[1] < length(labels)) {
          bins.nodata <- length(labels) - dim(summary(glm.ob)$coef)[1]
          if (bins.nodata >= length(ns)/2) {
            message(paste0(which.area, " ", which.time, 
                          " has no observation in more than half of the age bins, set to NA"))
            return(rep(NA, 5))
          }
          ns.comb <- ns
          ns.comb[length(ns) - bins.nodata] <- sum(ns[c(length(ns) - bins.nodata):length(ns)])
          message(paste0(which.area, " ", which.time, " has no observations in ", 
                        bins.nodata, " age bins. They are collapsed with previous bins"))
          return(get.est.withNA(glm.ob, labels, ns.comb))
        }
        return(get.est.withNA(glm.ob, labels, which.area, which.time))
      }
    }
    
    # Function to extract estimates and calculate U5MR
    get.est.withNA <- function(glm.ob, labels_input, region_val, year_val) {
      K <- length(labels)
      V <- matrix(0, K, K)
      betas <- rep(NA, K)
      labels <- paste("factor(age0)", labels, sep = "")
      colnames(V) <- rownames(V) <- labels
      names(betas) <- labels
      
      V2 <- stats::vcov(glm.ob)
      V[rownames(V2), colnames(V2)] <- V2
      betas2 <- coef(glm.ob)
      labels <- names(betas2)
      age_group <- gsub("factor\\(age0\\)", "", labels)
      
      # Dynamically derive ns for this subset of age groups
      ns <- sapply(age_group, function(lab) {
        if (lab == "0") return(1)
        bounds <- as.numeric(strsplit(lab, "-")[[1]])
        return(bounds[2] - bounds[1] + 1)
      })
      
      betas[names(betas2)] <- betas2
      probs <- plogis(betas)  # Use plogis instead of expit
      
      # Monte Carlo integration for U5MR
      set.seed(123)
      B <- 2000
      
      beta_draws <- MASS::mvrnorm(B, mu = betas, Sigma = V)
      probs_draws <- plogis(beta_draws)  # Use plogis
      u5mr_draws <- apply(probs_draws, 1, function(pr) 1 - prod((1-pr)^ns, na.rm = TRUE))
      
      mean.est <- mean(u5mr_draws, na.rm = TRUE)
      sd.est <- sd(u5mr_draws, na.rm = TRUE)
      CI <- quantile(u5mr_draws, c(0.025, 0.975), na.rm = TRUE)
      
      output_df <- data.table(
        region = region_val,
        years = year_val,
        age_group = age_group,
        ns = ns,
        probs = probs,
        hs = -log(1 - probs),
        u5m.m = mean.est, 
        u5m.lo = CI[1], 
        u5m.hi = CI[2], 
        u5m.sd = sd.est
      )
      
      return(output_df)
    }
    
    # Create all region-year combinations
    combo_df <- expand.grid(
      region = regions_list,
      years = years,
      stringsAsFactors = FALSE
    )
    
    # Apply function to all combinations
    results <- purrr::pmap_dfr(
      combo_df,
      function(region, years) {
        region.time.HT.withNA(region, years)
      }
    )
    
    results <- data.table(isid = isi, results)
    results.all[[isi]] <- results
  }
  
  return(data.table::rbindlist(results.all))
}

# ===============================================================================
# MAIN ANALYSIS PIPELINE
# ===============================================================================

cat("Starting DHS hazard extraction...\n")

# Initialize storage
dhs.ests <- list()
dhs.files <- list.files("../../data/DHS", full.names = TRUE, pattern = "\\.dta$", ignore.case = TRUE)

if (length(dhs.files) == 0) {
  stop("No DHS .dta files found in data/DHS directory")
}

# Process each DHS file
for (files.br in dhs.files) {
  print("########################################################")
  df <- haven::read_dta(files.br, n_max = 1)
  
  # Handle different state variable names
  if ("sstate" %in% names(df)) {
    df$v024 <- df$sstate
  }
  
  survey_year <- df$v007[1]
  sid <- df$v000[1]
  cat("File:", files.br, "\n")
  cat("Survey code:", sid, "\n")
  cat("Survey year:", survey_year, "\n")
  
  # Set time window
  beg.year <- survey_year - 25  # 25 years before survey
  end.year <- survey_year       # survey year
  country <- sid
  
  message("Reading birth recode and extracting births")
  
  # Use SUMMER::getBirths to process the data
  dat.tmp <- SUMMER::getBirths(
    filepath = files.br,
    surveyyear = survey_year,
    year.cut = seq(beg.year, end.year, by = 5),
    compact = TRUE,
    variables = c("caseid", "v001", "v002", "v004", "v005",
                 "v021", "v022", "v023", "v024",
                 "v025", "v139", "bidx"),
    compact.by = c("v001", "v024", "v025", "v005")
  ) |> 
  dplyr::mutate(isid = sid, survey = survey_year)
  
  dhs.ests[[paste(sid, survey_year)]] <- dat.tmp
  message("DONE!")
}

# Combine all DHS estimates
dhs.ests.df <- dplyr::bind_rows(dhs.ests)

# Standardize column names
dhs.ests.df <- dhs.ests.df[, c("v001", "age", "time", "total", "died", "v005", 
                               "strata", "v025", "isid", "v024", "survey")]
colnames(dhs.ests.df) <- c("cluster", "age", "years", "total",
                          "Y", "v005", "strata", "urban", "isid", "region", "survey")

# Process year coding
mod.data <- dhs.ests.df |> 
  dplyr::mutate(
    years = as.numeric(substring(paste(years), 1, 2)), 
    years = ifelse(years < 25, years + 2002, years + 1902)
  )

# Extract hazards (excluding problematic surveys if needed)
hazards.df <- getDirectHazards(mod.data |> dplyr::filter(isid != "MW4"))

# ===============================================================================
# POST-PROCESSING AND VISUALIZATION
# ===============================================================================

# Create country mappings
nat.u5mr <- hazards.df |>
  dplyr::mutate(id = substr(isid, 1, 2)) |>
  dplyr::mutate(country = dplyr::case_when(
    id == "BF" ~ "Burkina Faso", 
    id == "MW" ~ "Malawi", 
    id == "NG" ~ "Nigeria", 
    id == "NI" ~ "Niger", 
    id == "TZ" ~ "Tanzania"
  )) |> 
  dplyr::filter(region == "All") |> 
  dplyr::select(country, years, u5m.m, u5m.lo, u5m.hi, u5m.sd) |>
  dplyr::distinct() |> 
  dplyr::arrange(years)

# Create visualization of trends
create_u5mr_trend_plot <- function(data, output_path) {
  if (!dir.exists(dirname(output_path))) {
    dir.create(dirname(output_path), recursive = TRUE)
  }
  
  png(output_path, width = 2500, height = 1500, res = 300)
  
  p <- ggplot2::ggplot(
    data = subset(data, country != "Nigeria" & years > 1999), 
    ggplot2::aes(x = years, y = 1e3 * u5m.m, 
                ymax = 1e3 * u5m.hi, ymin = 1e3 * u5m.lo, 
                fill = country, color = country)
  ) + 
  ggplot2::geom_errorbar() + 
  ggplot2::geom_point(pch = 21, size = 2.5, color = "black") +  
  ggthemes::theme_hc() +
  ggplot2::facet_wrap(~country) + 
  ggplot2::ylim(0, 220) + 
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(
    title = "Direct U5MR estimates using DHS",
    x = "Year",
    y = "U5MR per 1,000 live births"
  )
  
  print(p)
  dev.off()
}

save(hazards.df, dhs.ests.df, file = "../../results/dhs_hazard_data.rda")

# Create trend plot
create_u5mr_trend_plot(nat.u5mr, "../../results/figs/dhs_u5mr_trends.png")

cat("DHS hazard extraction completed successfully!\n")
cat("Results saved to: ../../results/dhs_hazard_data.rda\n")
cat("Trend plot saved to: ../../results/figs/dhs_u5mr_trends.png\n")

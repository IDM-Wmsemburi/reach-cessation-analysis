# REACH Secondary Analysis: Azithromycin MDA Cessation Framework

Evidence-based cessation thresholds for azithromycin mass drug administration based on secondary analysis of REACH trial consortium data (AVENIR, CHAT, MORDOR I/II).

## Overview

This repository implements a Bayesian framework for developing cessation criteria for azithromycin mass drug administration (MDA) programs. The analysis integrates baseline mortality estimation, contextual factors, and probabilistic threshold derivation using data from five cluster-randomized trials across sub-Saharan Africa.

## Analysis Pipeline

The analysis follows a four-stage pipeline:

### 1. Data Preparation (`src/01_data_preparation/`)
- **DHS Hazard Extraction** (`1 prepare_dhs_hazards.r`): Extract age-specific mortality hazards from DHS birth histories using survey-weighted methods
- **Trial Prior Derivation** (`2 derive_trial_priors.r`): Generate data-driven priors from REACH trial placebo data using quasi-Poisson regression with cluster-robust standard errors

**Key Outputs:**
- `results/dhs_hazard_data.rda`: Country-specific neonatal and post-neonatal hazard estimates
- `results/trial_hazard_priors.rds`: Trial-country specific hazard priors for Bayesian updating
- `results/dhs_h0_pooled_model.rds`: Pooled regression linking neonatal to post-neonatal hazards

### 2. Baseline Mortality Estimation (`src/02_baseline_mortality/`)
- **Bayesian HSGP Model** (`bayesian_hsgp_model.r`): Estimate community-level baseline U5MR and IMR using spatial hierarchical models

**Method Details:**
- Hilbert Space Gaussian Process (HSGP) for computational efficiency
- Regularized Horseshoe for cluster-level deviations
- Age-split priors derived from DHS cumulative hazards 
- DHS neonatal bridge regression with full uncertainty propagation

**Key Outputs:**
- `results/posterior_hsgp_fit.rda`: Full posterior samples from Stan model
- `results/baseline_mortality_estimates.rda`: Community-level U5MR and IMR estimates
- `results/age_beta_priors_from_dhs.rda`: Beta priors for age splitting 

### 3. Contextual Analysis (`src/03_contextual_analysis/`)
- **Contextual Heterogeneity** (`contextual_heterogeneity.r`): Explore vaccination coverage, malaria burden, and other contextual factors using mixed-effects models

**Covariates Analyzed:**
- Vaccination: DPT1, DPT3, BCG, MCV1, Polio3 unmet need
- Malaria: ITN access/use, treatment coverage, incidence, mortality  
- Mortality: IHME U5M estimates, posterior U5M estimates, posterior IMR estimates

**Key Outputs:**
- `results/contextual_analysis_data.rda`: Merged contextual data with mortality estimates
- `results/mixed_effects_age_models.rda`: Age-stratified mixed effects models
- `results/contextual_anova_results.rda`: Statistical significance testing
- `results/figs/contextual_pvalue_heatmap.png`: P-value heatmap showing covariate and interaction effects

### 4. Threshold Estimation (`src/04_threshold_estimation/`)
- **Threshold Analysis** (`threshold_analysis.r`): Derive mortality-specific cessation thresholds using crossing-based probabilistic estimation

**Method Details:**
- Separate U5MR and IMR models
- Stan-based crossing detection where treatment effects vanish
- Hierarchical trial-country structure with mortality × treatment interactions
- Full uncertainty quantification via posterior sampling

**Key Outputs:**
- `results/threshold_models.rda`: Fitted U5MR and IMR threshold models
- `results/threshold_results.csv`: Threshold point estimates and credible intervals
- `results/threshold_coefficients.csv`: Model coefficients for all parameters
- `results/figs/threshold_analysis_plots.png`: Summary threshold and coefficient plots
- `results/figs/threshold_mortality_plots.png`: Observed vs predicted mortality relationships

## Data Sources

### Primary Analysis Data
* **AVENIR** (Niger): **2,158** communities, **619,228** children
* **CHAT** (Burkina Faso): **285** communities, **237,434** children
* **MORDOR I/II** (Niger): **594** communities, **400,111** children
* **MORDOR I** (Tanzania): **613** communities, **131,095** children
* **MORDOR I** (Malawi): **304** communities, **240,384** children

*(Total across settings: 3,954 communities, 1,628,252 children.)*

### Auxiliary Data
- **DHS birth histories**: 17 surveys across 5 countries for neonatal-postneonatal hazard relationships
- **IHME vaccination data**: MCV1, DTP1/3, BCG, Polio3 coverage rasters (2019)
- **Malaria indicators**: ITN use/access, treatment coverage, P. falciparum incidence/mortality (2019)

## Quick Start

### Prerequisites
- R ≥ 4.2.0
- Stan ≥ 2.30.0 via cmdstanr
- Required spatial packages: `sf`, `raster`
- Bayesian packages: `posterior`, `bayesplot`

## Repository Structure

```
├── src/                          # Analysis source code
│   ├── 01_data_preparation/      # DHS extraction and prior derivation
│   ├── 02_baseline_mortality/    # Spatial mortality modeling
│   ├── 03_contextual_analysis/   # Vaccination and context analysis  
│   ├── 04_threshold_estimation/  # Threshold derivation
│   └── utils/                    # Utility functions and configuration
├── models/stan/                  # Stan model files
├── data/                         # Data files (structure only - data not shared)
├── results/                      # Analysis outputs and figures
├── manuscript/                   # LaTeX manuscript and presentation
```

## Methods Summary

### Baseline Mortality Estimation
- **Spatial Model**: HSGP approximation with Matérn-3/2 covariance for computational scalability
- **Hierarchical Structure**: Trial-country intercepts + spatial smoothing + regularized horseshoe local deviations
- **Age Decomposition**: Beta priors derived from DHS hazard relationships
- **Neonatal Bridge**: Pooled regression with country fixed effects to translate post-neonatal to full U5MR

### Treatment Effect Modeling
- **Mortality Interactions**: Baseline mortality × treatment effects allow threshold identification
- **Hierarchical Structure**: Trial-country random effects account for between-study heterogeneity
- **Uncertainty Propagation**: Full Bayesian treatment preserves uncertainty through all analysis stages

### Threshold Derivation
- **Crossing Detection**: Probabilistic identification of mortality levels where treatment effects vanish
- **Stan Implementation**: Root-finding within generated quantities for robust threshold calculation
- **Coverage Assessment**: Proportion of posterior draws yielding valid crossings quantifies certainty

## Citation

If using this code or methodology, please cite:

> Msemburi, W. et al. Evidence-based cessation thresholds for azithromycin mass drug administration. *Manuscript in preparation* (2025).

## Contact

William Msemburi - william.msemburi@gatesfoundation.org 
Project Repository: https://github.com/IDM-Wmsemburi/reach-cessation-analysis

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
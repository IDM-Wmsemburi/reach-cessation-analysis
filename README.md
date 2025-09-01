# REACH Secondary Analysis: Azithromycin MDA Cessation Framework

Evidence-based cessation thresholds for azithromycin mass drug administration using comprehensive analysis of REACH trial consortium data.

## Overview

This repository implements a Bayesian framework for developing cessation criteria for azithromycin mass drug administration (MDA) programs. The analysis integrates baseline mortality estimation, contextual factors, and probabilistic threshold derivation using data from five cluster-randomized trials across sub-Saharan Africa.

**Approach**: Identifies mortality levels where azithromycin treatment benefits vanish through crossing-based threshold detection and comprehensive sensitivity analysis across multiple modeling frameworks.

## Data Sources

- **AVENIR** (Niger): 2,158 communities, 619,228 children
- **CHAT** (Burkina Faso): 285 communities, 237,434 children  
- **MORDOR I/II** (Niger): 594 communities, 400,111 children
- **MORDOR I** (Tanzania): 613 communities, 131,095 children
- **MORDOR I** (Malawi): 304 communities, 240,384 children

**Total**: 3,954 communities, 1.63M children, with auxiliary data from 17 DHS surveys and contextual indicators.

## Analysis Framework

### 1. Data Preparation (`src/01_data_preparation/`)
- **DHS Hazard Extraction**: Age-specific mortality patterns from survey data with uncertainty quantification
- **Trial Prior Derivation**: Baseline hazard priors from placebo clusters using quasi-Poisson regression

### 2. Baseline Mortality Estimation (`src/02_baseline_mortality/`)
- **Bayesian HSGP Model**: Spatially coherent cluster-level mortality using Hilbert Space Gaussian Processes
- **Age Decomposition**: DHS neonatal bridge regression linking post-neonatal to full U5MR
- **Uncertainty Propagation**: Full posterior distributions for all mortality estimates

### 3. Contextual Analysis (`src/03_contextual_analysis/`)
- **Heterogeneity Assessment**: Mixed-effects models testing vaccination gaps, malaria burden, and mortality predictors
- **Treatment Interactions**: Systematic evaluation of context-treatment effect modification
- **Statistical Testing**: Likelihood ratio tests comparing nested model structures

### 4. Threshold Estimation (`src/04_threshold_estimation/`)
- **Interaction Models**: Mortality × treatment effects allowing threshold identification
- **Crossing Detection**: Probabilistic identification where treatment/placebo curves converge
- **Dual Methods**: Both analytical (βₜ + βₘₜU* = 0) and crossing-based threshold estimation

### 5. Comprehensive Sensitivity Analysis (`src/05_sensitivity_analysis/`)
- **Multiple Frameworks**: Stan (variable/fixed baseline) and INLA implementations
- **Sequential Exclusions**: Systematic removal of each trial location to test robustness
- **Method Comparison**: Analytical vs crossing-based threshold calculation
- **Uncertainty Quantification**: Full posterior distributions across all scenarios

## Key Results

### Primary Thresholds
- **Both U5MR and IMR thresholds** derived using crossing-based method with full uncertainty quantification
- **Complete posterior distributions** available for risk-based decision making

### Sensitivity Analysis Summary
- **8 modeling scenarios** tested across Stan and INLA frameworks  
- **Low coefficient of variation** demonstrates robust estimates across approaches
- **Sequential exclusions** of trial locations show limited impact on threshold estimates
- **Method convergence**: Crossing method consistently higher than analytical approach

### Statistical Evidence
- Clear treatment-mortality interactions (negative interaction coefficients)
- Convergence diagnostics: R̂ < 1.01, effective sample sizes > 1000
- Contextual factors explain baseline variation but don't systematically modify treatment effects

## Repository Structure

```
├── src/
│   ├── 01_data_preparation/     # DHS extraction and prior derivation
│   ├── 02_baseline_mortality/   # Spatial mortality modeling  
│   ├── 03_contextual_analysis/  # Vaccination and context analysis
│   ├── 04_threshold_estimation/ # Threshold derivation
│   ├── 05_sensitivity_analysis/ # Comprehensive robustness testing
│   └── utils/                   # Utility functions
├── models/stan/                 # Stan model files
├── data/                        # Input data files
│   ├── Clean/                   # Processed trial data
│   ├── DHS/                     # Birth history surveys
│   ├── IHME/                    # Vaccination coverage rasters
│   ├── Malaria/                 # Malaria burden indicators
│   └── shps/                    # Administrative shapefiles
├── results/                     # Analysis outputs and figures
│   ├── threshold_results.csv    # Primary threshold estimates
│   ├── enhanced_threshold_analysis.rda  # Sensitivity analysis results
│   └── figs/                    # Visualizations
└── manuscript/                  # LaTeX manuscript and presentation
```

## Quick Start

```r
# Load main results
load("results/baseline_mortality_estimates.rda")
load("results/threshold_models.rda")

# View primary thresholds
read.csv("results/threshold_results.csv")

# Sensitivity analysis
load("results/enhanced_threshold_analysis.rda")
```

## Policy Implications

- **Data-driven thresholds** provide evidence-based foundation compared to current expert opinion-based guidance
- **Robust estimates** across modeling approaches provide confidence for implementation  
- **Resistance monitoring** remains essential regardless of threshold approach
- **Full uncertainty quantification** enables risk-based decision making

## Methods Summary

- **Bayesian spatial modeling** with HSGP approximation for computational efficiency
- **DHS demographic anchoring** provides age-pattern priors and neonatal bridge relationships
- **Hierarchical trial-country structure** accounts for between-study heterogeneity
- **Crossing-based thresholds** identify mortality levels where treatment effects vanish
- **Comprehensive sensitivity analysis** tests robustness across modeling assumptions and data subsets

## Citation

```
Msemburi, W. et al. Evidence-based cessation thresholds for azithromycin mass drug administration 
based on REACH trial consortium data. Manuscript in preparation (2025).
```

## Contact

William Msemburi - william.msemburi@gatesfoundation.org

## License

MIT License - see [LICENSE](LICENSE) file for details.
# Results Directory

This directory contains outputs from the REACH cessation analysis pipeline. Results are organized by analysis stage and include both intermediate files and final outputs.

## Directory Structure

```
results/
├── figs/                          # All generated figures
│   ├── dhs_u5mr_trends.png           # DHS mortality trends
│   ├── correlation_matrix.png        # Covariate correlation heatmap
│   ├── contextual_pvalue_heatmap.png # Statistical significance results
│   ├── scatter_*.png                 # Individual covariate scatter plots  
│   ├── threshold_analysis_plots.png  # Summary threshold plots
│   └── threshold_mortality_plots.png # Mortality vs threshold plots
├── [Stage 1] Data Preparation Outputs
├── dhs_hazard_data.rda            # DHS hazard estimates by country/year/age
├── dhs_h0_pooled_model.rds        # Neonatal hazard regression model
├── trial_hazard_priors.rds        # Trial-derived prior distributions
├── [Stage 2] Baseline Mortality Outputs  
├── posterior_hsgp_fit.rda         # Full Stan model fit object
├── baseline_mortality_estimates.rda # Community-level U5MR/IMR estimates
├── age_beta_priors_from_dhs.rda   # Age-split prior distributions
├── [Stage 3] Contextual Analysis Outputs
├── contextual_analysis_data.rda   # Trial data with contextual covariates
├── mixed_effects_age_models.rda   # Mixed effects model results
├── contextual_anova_results.rda   # ANOVA comparison results  
├── correlation_analysis.rda       # Covariate correlation matrices
├── [Stage 4] Threshold Analysis Outputs
├── threshold_models.rda           # Stan threshold models (U5MR & IMR)
├── threshold_analysis_results.rda # Processed threshold estimates
├── threshold_results.csv          # Summary threshold table (CSV)
├── threshold_coefficients.csv     # Model coefficient estimates (CSV)
└── package_installation_log.csv   # Package installation record
```

## File Descriptions

### Stage 1: Data Preparation

#### `dhs_hazard_data.rda`
**Content**: Age-specific hazard estimates from DHS birth histories
**Variables**:
- `isid`: Survey identifier  
- `region`: Geographic region
- `years`: Calendar year
- `age_group`: Age interval (0, 1-11, 12-23, 24-35, 36-47, 48-59)
- `hs`: Hazard estimate
- `u5m.m`, `u5m.lo`, `u5m.hi`: U5MR estimates with 95% CI

#### `dhs_h0_pooled_model.rds`
**Content**: Pooled regression model linking neonatal to post-neonatal hazards
**Structure**: List with coefficients, variance-covariance matrix, country fixed effects

#### `trial_hazard_priors.rds`  
**Content**: Trial-country specific hazard priors derived from placebo data
**Variables**:
- `trial_country`: Trial-country combination
- `log_h1_11_prior`, `log_h1_11_prior_sd`: 1-11 month hazard prior
- `log_h12_59_prior`, `log_h12_59_prior_sd`: 12-59 month hazard prior  
- `log_h1_59_prior`, `log_h1_59_prior_sd`: 1-59 month hazard prior

### Stage 2: Baseline Mortality Estimation

#### `posterior_hsgp_fit.rda`
**Content**: Complete Stan model fit object from spatial hierarchical model
**Size**: Large (~100MB-1GB depending on model complexity)
**Contains**: MCMC chains, diagnostics, full posterior samples

#### `baseline_mortality_estimates.rda`
**Content**: Community-level mortality estimates with uncertainty
**Key Variables**:
- `community_id`: Cluster identifier
- `post_u5mr`, `post_u5mr_sd`: U5MR posterior mean and SD
- `post_imr`, `post_imr_sd`: IMR posterior mean and SD  
- `post_u5mr_log`, `post_u5mr_log_sd`: Log-scale estimates
- Hazard decomposition: `H1_11_pred_mean`, `H12_59_pred_mean`, etc.

#### `age_beta_priors_from_dhs.rda`
**Content**: Age-split prior parameters derived from DHS hazard ratios
**Variables**: Beta distribution parameters `a0`, `b0` by trial-country

### Stage 3: Contextual Analysis

#### `contextual_analysis_data.rda`
**Content**: Trial dataset augmented with extracted contextual covariates
**Added Variables**:
- Vaccination: `dpt1`, `dpt3`, `bcg1`, `mcv1`, `polio3`
- Malaria: `mal.mort`, `mal.inc`, `itn.use`, `itn.acc`
- Derived: `gm_cov` (geometric mean), `pc1` (first principal component)

#### `mixed_effects_age_models.rda`
**Content**: Mixed effects models testing covariate and treatment interactions
**Structure**: Nested tibble with models by covariate and trial-country

#### `contextual_anova_results.rda`
**Content**: ANOVA comparisons testing statistical significance of effects
**Variables**: p-values, AIC differences, effect sizes by covariate type

#### `correlation_analysis.rda`
**Content**: Correlation matrices between covariates and mortality rates
**Includes**: Spearman correlations, significance tests

### Stage 4: Threshold Analysis

#### `threshold_models.rda`
**Content**: Stan models estimating cessation thresholds for U5MR and IMR
**Includes**: Model fits, prediction grids, Stan data objects

#### `threshold_analysis_results.rda`
**Content**: Processed threshold estimates with uncertainty quantification
**Key Results**:
- Median thresholds with 95% credible intervals
- Coefficient estimates for mortality-treatment interactions
- Model diagnostics and convergence statistics

#### `threshold_results.csv` ⭐
**Content**: **Main analysis output** - cessation threshold estimates
**Columns**:
- `model`: U5MR or IMR
- `threshold_med_1000`: Median threshold per 1,000 live births
- `threshold_lo_1000`, `threshold_hi_1000`: 95% credible interval bounds
- `n_crossings`: Number of posterior draws showing crossings
- `p_crossings`: Proportion of draws with identifiable thresholds

#### `threshold_coefficients.csv`
**Content**: Model coefficient estimates and uncertainties
**Interpretation**: 
- `beta_mortality`: Baseline mortality effect on death rates
- `beta_trt`: Treatment main effect
- `beta_trt_mortality`: Treatment-mortality interaction (threshold slope)

## Key Outputs Summary

### Primary Results
1. **Cessation Thresholds**: `threshold_results.csv`
   - U5MR threshold: ~XXX per 1,000 (95% CI: XXX-XXX)
   - IMR threshold: ~XXX per 1,000 (95% CI: XXX-XXX)

2. **Contextual Modifiers**: `contextual_anova_results.rda`
   - Significant treatment interactions by vaccination coverage
   - Geographic variation in treatment effects

3. **Baseline Mortality Maps**: `baseline_mortality_estimates.rda`
   - Community-specific U5MR and IMR estimates
   - Spatial clustering patterns

### Diagnostic Outputs
- **Model Convergence**: Check Rhat values in model fit objects
- **Threshold Coverage**: Proportion of posterior draws with identifiable crossings
- **Spatial Validity**: Residual spatial patterns in baseline estimates

## Reproducibility Notes

### File Dependencies
Analysis files depend on successful completion of prior stages:
```
Stage 1 → Stage 2 → Stage 3 → Stage 4
```

### Computational Requirements
- **Stage 1**: ~5-15 minutes (DHS processing)
- **Stage 2**: ~30-120 minutes (Bayesian spatial model)  
- **Stage 3**: ~10-30 minutes (Contextual analysis)
- **Stage 4**: ~20-60 minutes (Threshold models)

### Memory Requirements
- Minimum: 8GB RAM
- Recommended: 16GB+ RAM for Stage 2 spatial modeling
- Storage: ~2-5GB for all outputs

### Platform Differences
Results may vary slightly across platforms due to:
- Random number generation differences
- Compiler optimizations (Stan models)
- Floating point precision differences

For exact reproducibility, use the same R version, package versions, and random seeds as documented in the analysis configuration.

## Interpreting Results

### Threshold Estimates
- Values represent baseline U5MR/IMR levels below which azithromycin MDA shows no significant mortality benefit
- Confidence intervals reflect uncertainty in threshold location
- Coverage indicates proportion of scenarios where thresholds can be identified

### Contextual Modifiers  
- Significant interactions suggest thresholds vary by context
- p-values indicate strength of evidence for modification
- Effect sizes show magnitude of treatment heterogeneity

### Baseline Mortality Estimates
- Represent community-specific mortality in absence of treatment
- Account for spatial correlation and trial-country differences  
- Include both aleatory and epistemic uncertainty

For detailed interpretation guidance, see the main analysis documentation and associated manuscript.
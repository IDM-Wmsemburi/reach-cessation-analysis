
data {
  int<lower=1> C;
  int<lower=1> K;
  int<lower=1> Kc;
  array[C] int<lower=1,upper=K>  group_id;
  array[C] int<lower=1,upper=Kc> country_id;
  int<lower=0> C_obs;
  array[C_obs] int<lower=1,upper=C> idx_obs;

  array[C] int<lower=0> y;
  vector[C] logE;
  array[C] int<lower=0,upper=1> is_obs;

  // HSGP
  int<lower=1> M_hsgp;
  matrix[C, M_hsgp] PHI;
  vector<lower=0>[M_hsgp] lambda_hsgp; // omega^2

  int<lower=1> n_groups;
  array[n_groups] int<lower=1> group_basis_start;
  array[n_groups] int<lower=1> group_basis_size;

  real<lower=0> slab_scale;

  // DHS neonatal bridge
  matrix[Kc,3] ab_mean;
  array[Kc] matrix[3,3] L_ab;
  vector<lower=0>[Kc] sigma0_hat;       // fixed residual SD

  // Group mean priors (per-year log-rate)
  array[K] real mu_logLambda_prior;
  array[K] real<lower=0> sd_logLambda_prior;

  // Age-split prior (DHS-only)
  array[K] real<lower=0> a0;
  array[K] real<lower=0> b0;
}
parameters {
  vector[K] mu_g;
  real<lower=0.1, upper=5.0> sigma;          // GP marginal SD
  real<lower=0.1, upper=10.0> rho;           // GP range
  vector[M_hsgp] beta_hsgp;                  // HSGP weights

  vector[C] z_v;
  vector<lower=0>[C] lambda;
  vector<lower=0>[K] tau_g;

  vector<lower=0,upper=1>[K] p_g;            // age split (sampled)
}
transformed parameters {
  vector[C] u;
  vector[C] v;
  vector[C] log_lambda_total;
  vector[C] log_lambda_1_11;
  vector[C] log_lambda_12_59;

  // Matérn 3/2 in 2D: spectral density S(w) ∝ kappa^3 (kappa^2 + w^2)^(-2.5)
  real<lower=0> kappa = sqrt(3.0) / rho;
  {
    vector[M_hsgp] sqrt_eigenvalues;
    vector[M_hsgp] beta_scaled;
    for (m in 1:M_hsgp) {
      real sqrt_spd = sigma
        * sqrt( 3.0 * pi()
                * pow(kappa, 3)
                * pow(kappa^2 + lambda_hsgp[m], -2.5) );
      sqrt_eigenvalues[m] = sqrt_spd;
      beta_scaled[m] = sqrt_eigenvalues[m] * beta_hsgp[m];
    }
    u = PHI * beta_scaled;
  }
  
  for (i in 1:C) {
    real tau_lam = tau_g[group_id[i]] * lambda[i];
    real c2 = square(slab_scale);
    real lt2 = square(slab_scale * lambda[i]) / (c2 + square(tau_lam));
    real lambda_tilde = sqrt(fmax(lt2, 1e-10));
    v[i] = tau_g[group_id[i]] * lambda_tilde * z_v[i];
    log_lambda_total[i] = mu_g[group_id[i]] + u[i] + v[i];
  }
  
  // Age-specific log rates via sampled p_g
  for (i in 1:C) {
    int g = group_id[i];
    real p = p_g[g];
    log_lambda_1_11[i]  = log_lambda_total[i] + log(p)   + log(59.0/11.0);
    log_lambda_12_59[i] = log_lambda_total[i] + log1m(p) + log(59.0/48.0);
  }
}
model {
  // group means (from DHS H1_59 priors)
  for (g in 1:K)
    mu_g[g] ~ normal(mu_logLambda_prior[g], sd_logLambda_prior[g]);

  // age split: *only* DHS-implied Beta prior
  p_g ~ beta(a0, b0);

  // GP + horseshoe
  sigma    ~ normal(0, 1);
  rho      ~ lognormal(0, 1);
  beta_hsgp ~ normal(0, 1);
  z_v      ~ normal(0, 1);
  lambda   ~ student_t(3, 0, 1);
  tau_g    ~ student_t(3, 0, 1);

  // likelihood (placebo clusters only)
  y[idx_obs] ~ poisson_log(logE[idx_obs] + log_lambda_total[idx_obs]);
}
generated quantities {
  array[Kc] vector[3] ab_rng;             // country-specific coef draws
  array[C] real H1_11_true;
  array[C] real H12_59_true;
  array[C] real H1_total_true;
  array[C] real H1_11_pred;
  array[C] real H12_59_pred;
  array[C] real H1_total_pred;
  array[C] real H0_pred;
  array[C] real imr_0_11;
  array[C] real u5mr_0_59;

  // Draw country-level neonatal coefficients; residual SD fixed to sigma0_hat
  for (k in 1:Kc) {
    vector[3] mu_k = to_vector(ab_mean[k]');
    ab_rng[k] = multi_normal_cholesky_rng(mu_k, L_ab[k]);
  }

  // Predictive cluster hazards and mortality
  real c2 = square(slab_scale);
  for (i in 1:C) {
    int g  = group_id[i];
    int kc = country_id[i];

    // latent hazards (fitted)
    real logH1_total_true = log_lambda_total[i] + log(59.0/12.0);
    real logH1_11_true    = log_lambda_1_11[i] + log(11.0/12.0);
    real logH12_59_true   = log_lambda_12_59[i] + log(48.0/12.0);
    H1_total_true[i] = exp(logH1_total_true);
    H1_11_true[i]    = exp(logH1_11_true);
    H12_59_true[i]   = exp(logH12_59_true);

    // predictive new horseshoe innovation
    real lambda_new = sqrt(square(student_t_rng(3, 0, 1)));
    real z_new = normal_rng(0, 1);
    real tau_lam_new = tau_g[g] * lambda_new;
    real lt2_new = square(slab_scale * lambda_new) / (c2 + square(tau_lam_new));
    real lambda_tilde_new = sqrt(lt2_new);
    real v_new = tau_g[g] * lambda_tilde_new * z_new;

    real log_lambda_total_pred = mu_g[g] + u[i] + v_new;
    real p = p_g[g];
    real log_lambda_1_11_pred  = log_lambda_total_pred + log(p)   + log(59.0/11.0);
    real log_lambda_12_59_pred = log_lambda_total_pred + log1m(p) + log(59.0/48.0);

    real logH1_total_pred = log_lambda_total_pred + log(59.0/12.0);
    real logH1_11_pred    = log_lambda_1_11_pred + log(11.0/12.0);
    real logH12_59_pred   = log_lambda_12_59_pred + log(48.0/12.0);
    H1_total_pred[i] = exp(logH1_total_pred);
    H1_11_pred[i]    = exp(logH1_11_pred);
    H12_59_pred[i]   = exp(logH12_59_pred);

    // Neonatal hazard via DHS bridge
    {
      vector[3] ab = ab_rng[kc];
      real mu_logH0 = ab[1]
                    + ab[2] * log(H1_11_pred[i])
                    + ab[3] * log(H12_59_pred[i]);
      real logH0    = normal_rng(mu_logH0, sigma0_hat[kc]);  // fixed residual SD
      H0_pred[i]    = exp(logH0);
    }

    // Mortality rates
    imr_0_11[i]  = 1 - exp(-(H0_pred[i] + H1_11_pred[i]));
    u5mr_0_59[i] = 1 - exp(-(H0_pred[i] + H1_11_pred[i] + H12_59_pred[i]));
  }
}


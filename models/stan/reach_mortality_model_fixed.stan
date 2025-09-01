
functions {
  real crossing_from_preds(int Np, array[,] real pred, vector log_x) {
    real neg_inf = negative_infinity();
    for (u in 1:(Np-1)) {
      real d0 = log(pred[u,2])   - log(pred[u,1]);
      real d1 = log(pred[u+1,2]) - log(pred[u+1,1]);
      if (abs(d0) < 1e-12) return exp(log_x[u]);
      if ((d0 < 0 && d1 > 0) || (d0 > 0 && d1 < 0)) {
        real t  = d0 / (d0 - d1);
        real lx = log_x[u] + t * (log_x[u+1] - log_x[u]);
        return exp(lx);
      }
    }
    return neg_inf;
  }
}
data {
  int<lower=1> N, J, N_pred;
  array[N] int<lower=0> deaths;
  vector<lower=0>[N] person_time;
  array[N] int<lower=0,upper=1> trt;
  array[N] int<lower=1,upper=J> group;
  vector[N] log_mortality_mean;
  vector<lower=0>[N] log_mortality_sd;
  vector[N_pred] log_mortality_pred;
}
parameters {
  // Hierarchical group effects
  vector[J] alpha_group_raw;
  real mu_alpha_group;
  real<lower=0> sigma_group;
  
  // Model coefficients
  real beta_mortality;
  real beta_trt;
  real beta_trt_mortality;
}
transformed parameters {
  vector[J] alpha_group = mu_alpha_group + sigma_group * alpha_group_raw;
}
model {
  // Priors
  alpha_group_raw ~ std_normal();
  mu_alpha_group ~ normal(0, 1);
  sigma_group ~ normal(0, 1);
  
  beta_mortality ~ normal(0, 1);
  beta_trt ~ normal(0, 0.5);
  beta_trt_mortality ~ normal(0, 0.5);
  
  // Likelihood
  vector[N] mu;
  for (i in 1:N) {
    mu[i] = alpha_group[group[i]] + 
            beta_mortality * log_mortality_mean[i] +
            trt[i] * (beta_trt + beta_trt_mortality * log_mortality_mean[i]) +
            log(person_time[i]);
  }
  deaths ~ poisson_log(mu);
}
generated quantities {
  // Predictions
  array[N_pred, 2] real pred_lambda;
  
  // Crossing-based threshold
  real threshold;
  int has_cross;
  
  real mean_alpha = mean(alpha_group);
  
  // Generate predictions
  for (u in 1:N_pred) {
    real mu_base = mean_alpha +
                   beta_mortality * log_mortality_pred[u] +
                   log(365.25);
    real mu_trt = mu_base +
                  beta_trt +
                  beta_trt_mortality * log_mortality_pred[u];
    
    pred_lambda[u, 1] = exp(mu_base);  // placebo
    pred_lambda[u, 2] = exp(mu_trt);   // treatment
  }
  
  // Calculate crossing threshold
  threshold = crossing_from_preds(N_pred, pred_lambda, log_mortality_pred);
  has_cross = (threshold != negative_infinity());
}


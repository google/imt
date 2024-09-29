/*
# Copyright 2024 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
*/

data {
  // Data dimensions
  int<lower=0> N;        // Number of observations
  int<lower=0, upper=1> treatment[N]; // Treatment indicator (0 or 1)
  int<lower=1> K;        // Number of predictors (columns in X)
  // Data
  matrix[N, K] X;        // Design matrix
  vector<lower=0>[N] y;  // Observed outcome (non-negative)
  // Prior parameters
  real mean_alpha_logit;
  real<lower=0> sd_alpha_logit;
  vector[K] mean_beta_logit;
  vector<lower=0>[K] sd_beta_logit;
  real tau_mean_logit;
  real<lower=0> tau_sd_logit;
  real mean_tau;
  real<lower=0> sigma_tau;
  // Flag for running estimation (0: no, 1: yes)
  int<lower=0, upper=1> run_estimation;
}

transformed data {
  matrix[N, K] X_std;       // Standardized design matrix
  int N_positive;           // Number of positive y values
  vector[N] log_y;          // Log-transformed y values (including zeros)
  real log_y_mean;
  real log_y_sd;

  // Standardize the design matrix
  for (k in 1:K) {
    X_std[, k] = (X[, k] - mean(X[, k])) / sd(X[, k]);
  }
  // Count positive values and log-transform y
  N_positive = 0;
  for (n in 1:N) {
    if (y[n] > 0) {
      N_positive += 1;
      log_y[n] = log(y[n]);
    } else {
      log_y[n] = negative_infinity();
    }
  }
  // Calculate mean and sd of log-transformed positive y values
  {
    vector[N_positive] positive_log_y;
    int pos_index = 1;
    for (n in 1:N) {
      if (y[n] > 0) {
        positive_log_y[pos_index] = log_y[n];
        pos_index += 1;
      }
    }
    log_y_mean = mean(positive_log_y);
    log_y_sd = sd(positive_log_y);
  }
}

parameters {
  // Zero-inflation (logit) parameters
  real alpha_zero;
  vector[K] beta_zero;
  real tau_zero;

  // Log-normal parameters
  real alpha_lnorm;
  vector[K] beta_lnorm;
  real tau_lnorm;
  real<lower=0> sigma_lnorm;
}

model {
  vector[N] theta_zero;
  vector[N] mu_lnorm;
  // Priors for zero-inflation parameters
  alpha_zero ~ normal(mean_alpha_logit, sd_alpha_logit);
  beta_zero ~ normal(mean_beta_logit, sd_beta_logit);
  tau_zero ~ normal(tau_mean_logit, tau_sd_logit);
  // Priors for log-normal parameters
  alpha_lnorm ~ normal(log_y_mean, 1);
  beta_lnorm ~ normal(0, 0.5);
  tau_lnorm ~ normal(mean_tau, sigma_tau);
  sigma_lnorm ~ normal(0, 0.5);
  // Linear predictors
  theta_zero = alpha_zero + X_std * beta_zero + tau_zero * to_vector(treatment);
  mu_lnorm = alpha_lnorm + X_std * beta_lnorm + tau_lnorm * to_vector(treatment);
  if (run_estimation==1) {
    // Likelihood
    for (n in 1:N) {
      if (y[n] == 0) {
        target += bernoulli_logit_lpmf(1 | theta_zero[n]);
      } else {
        target += bernoulli_logit_lpmf(0 | theta_zero[n]);
        target += normal_lpdf((log_y[n] - log_y_mean) / log_y_sd | mu_lnorm[n], sigma_lnorm);
      }
    }
  }
}

generated quantities {
  vector[N] y_pred;
  vector[N] y0; // Potential outcome under control
  vector[N] y1; // Potential outcome under treatment
  real ATE;
  real tau_prob_zero;

  int<lower=0, upper=1> y0_zero[N];
  int<lower=0, upper=1> y1_zero[N];

  for (n in 1:N) {
    // Linear predictors for both treatment and control
    real theta_zero_0 = alpha_zero + X_std[n] * beta_zero;
    real theta_zero_1 = theta_zero_0 + tau_zero;
    real mu_lnorm_0 = alpha_lnorm + X_std[n] * beta_lnorm;
    real mu_lnorm_1 = mu_lnorm_0 + tau_lnorm;

    // Predictions for actual treatment
    if (bernoulli_logit_rng(theta_zero_0 + tau_zero * treatment[n]) == 1) {
      y_pred[n] = 0;
    } else {
      y_pred[n] = exp(normal_rng(mu_lnorm_0 + tau_lnorm * treatment[n], sigma_lnorm) * log_y_sd + log_y_mean);
    }

    // Potential outcomes
    if (bernoulli_logit_rng(theta_zero_0) == 1) {
      y0[n] = 0;
    } else {
      y0[n] = exp(normal_rng(mu_lnorm_0, sigma_lnorm) * log_y_sd + log_y_mean);
    }

    if (bernoulli_logit_rng(theta_zero_1) == 1) {
      y1[n] = 0;
    } else {
      y1[n] = exp(normal_rng(mu_lnorm_1, sigma_lnorm) * log_y_sd + log_y_mean);
    }

    // Potential outcomes for zero inflation
    y0_zero[n] = bernoulli_logit_rng(theta_zero_0);
    y1_zero[n] = bernoulli_logit_rng(theta_zero_1);
  }
  ATE = mean(y1 - y0);
  tau_prob_zero = mean(y1_zero) - mean(y0_zero);
}

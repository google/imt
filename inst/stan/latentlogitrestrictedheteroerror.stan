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
  int<lower=1> C;             // Number of calls (units)
  int<lower=0> K;             // Number of covariates
  matrix[C, K] X;             // Covariate matrix
  vector[C] treat;            // Treatment indicator (0 or 1)

  // Satisfaction Data
  array[C] int<lower=0> k;    // Number of 'dissatisfied' (1) ratings for each call
  array[C] int<lower=1> N;    // Total number of ratings for each call
  
  // Difficulty Data
  vector<lower=0, upper=1>[C] d_obs;  // Observed difficulty for each call

  // Priors for the latent logistic regression
  real mean_alpha;
  real<lower=0> sd_alpha;
  vector[K] mean_beta;
  vector[K] sd_beta;
  real tau_mean;
  real<lower=0> tau_sd;

  // Flag for running estimation
  int<lower=0, upper=1> run_estimation;

  int<lower=0, upper=1> run_covariates; // Flag for running covariates
  int<lower=0, upper=1> run_treatment;  // Flag for running treatment
}

transformed data {
  matrix[C, K] X_std;
  vector[K] mean_X;
  vector[K] sd_X;
  for (k_i in 1:K) {
    mean_X[k_i] = mean(X[, k_i]);
    sd_X[k_i] = sd(X[, k_i]);
    if (sd_X[k_i] > 1e-6) {
      X_std[, k_i] = (X[, k_i] - mean_X[k_i]) / sd_X[k_i];
    } else {
      X_std[, k_i] = rep_vector(0.0, C);
    }
  }
}

parameters {
  // Latent dissatisfaction regression parameters
  real alpha;
  vector[K] beta;
  real tau;

  // Parameters for the error rate functions
  // logit(epsilon) = alpha_err + beta_err * d_obs
  real alpha_0;
  real<lower=0> beta_0; // Slope for FPR
  real alpha_1;
  real<lower=0> beta_1; // Slope for FNR
}

model {
  real log_prob_if_dissatisfied;
  real log_prob_if_satisfied;
  real epsilon_0_c;
  real epsilon_1_c;

  // Priors
  alpha ~ normal(mean_alpha, sd_alpha);
  
  if (run_covariates == 1) {
    beta ~ normal(mean_beta, sd_beta);
  } else {
    beta ~ normal(0, 0.1); // Pin unused beta to 0
  }

  if (run_treatment == 1) {
    tau ~ normal(tau_mean, tau_sd);
  } else {
    tau ~ normal(0, 0.1); // Pin unused tau to 0
  }

  alpha_0 ~ normal(-1.5, 1);
  beta_0 ~ lognormal(0, 0.5);
  alpha_1 ~ normal(-1.5, 1);
  beta_1 ~ lognormal(0, 0.5);

  // Likelihood
  if (run_estimation == 1) {
    vector[C] theta_c = rep_vector(alpha, C);
    if (run_covariates == 1) {
      theta_c += X_std * beta;
    }
    if (run_treatment == 1) {
      theta_c += tau * treat;
    }
    for (c in 1:C) {
      // Calculate call-specific error rates directly from observed difficulty (capped at 0.5)
      epsilon_0_c = inv_logit(alpha_0 + beta_0 * d_obs[c]) / 2.0;
      epsilon_1_c = inv_logit(alpha_1 + beta_1 * d_obs[c]) / 2.0;

      log_prob_if_dissatisfied = bernoulli_logit_lpmf(1 | theta_c[c]) +
                                      binomial_lpmf(k[c] | N[c], 1 - epsilon_1_c);
      log_prob_if_satisfied = bernoulli_logit_lpmf(0 | theta_c[c]) +
                                   binomial_lpmf(k[c] | N[c], epsilon_0_c);
      
      target += log_sum_exp(log_prob_if_dissatisfied, log_prob_if_satisfied);
    }
  }
}

generated quantities {
  real eta;
  vector[C] prob_dissatisfied;
  matrix[C, 2] individualized_errors;
  
  // Calculate theta for the observed treatment assignment
  vector[C] theta_c_obs = rep_vector(alpha, C);
  
  // Calculate potential outcomes for the ATE (eta)
  vector[C] theta_c_treated;
  vector[C] theta_c_control;

  real log_prob_if_dissatisfied;
  real log_prob_if_satisfied;

  real epsilon_0_c;
  real epsilon_1_c;

  if (run_covariates == 1) {
    theta_c_obs += X_std * beta;
  }
  if (run_treatment == 1) {
    theta_c_obs += tau * treat;
  }

  // Build counterfactuals and eta
  if (run_treatment == 1) {
    // Build control (no treatment)
    theta_c_control = rep_vector(alpha, C);
    if (run_covariates == 1) {
      theta_c_control += X_std * beta;
    }
    
    // Build treated
    theta_c_treated = theta_c_control + tau;
    
    // Calculate ATE
    eta = mean(inv_logit(theta_c_treated)) - mean(inv_logit(theta_c_control));
    
  } else {
    // If no treatment, ATE is 0 and counterfactuals are just the observed
    eta = 0.0;
    theta_c_control = theta_c_obs;
    theta_c_treated = theta_c_obs;
  }

  for (c in 1:C) {
    epsilon_0_c = inv_logit(alpha_0 + beta_0 * d_obs[c]) / 2.0; // Using the fix from above
    epsilon_1_c = inv_logit(alpha_1 + beta_1 * d_obs[c]) / 2.0; // Using the fix from above

    individualized_errors[c, 1] = epsilon_0_c;
    individualized_errors[c, 2] = epsilon_1_c;

    // Use theta_c_obs for calculating posterior probability of the latent state
    log_prob_if_dissatisfied = bernoulli_logit_lpmf(1 | theta_c_obs[c]) +
                                    binomial_lpmf(k[c] | N[c], 1 - epsilon_1_c);
    log_prob_if_satisfied = bernoulli_logit_lpmf(0 | theta_c_obs[c]) +
                                 binomial_lpmf(k[c] | N[c], epsilon_0_c);
                                 
    prob_dissatisfied[c] = exp(log_prob_if_dissatisfied - log_sum_exp(log_prob_if_dissatisfied, log_prob_if_satisfied));
  }
}

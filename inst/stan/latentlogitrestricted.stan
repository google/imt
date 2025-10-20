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
  array[C] int<lower=0> k;    // Number of 'dissatisfied' (1) ratings for each call
  array[C] int<lower=1> N;    // Total number of ratings for each call

  // Priors for the latent logistic regression
  real mean_alpha;            // Prior mean for intercept
  real<lower=0> sd_alpha;     // Prior SD for intercept
  vector[K] mean_beta;        // Prior mean for covariates
  vector[K] sd_beta;          // Prior SD for covariates
  real tau_mean;              // Prior mean for treatment effect
  real<lower=0> tau_sd;       // Prior SD for treatment effect

  // Priors for the measurement error rates (default: Beta(1,1))
  real<lower=0> epsilon0_alpha; // Beta prior alpha for False Positive Rate
  real<lower=0> epsilon0_beta;  // Beta prior beta for False Positive Rate
  real<lower=0> epsilon1_alpha; // Beta prior alpha for False Negative Rate
  real<lower=0> epsilon1_beta;  // Beta prior beta for False Negative Rate

  // Flag for running estimation (0: prior only, 1: full)
  int<lower=0, upper=1> run_estimation;
}

transformed data {
  matrix[C, K] X_std;  // Standardized covariates
  vector[K] mean_X;    // Means of covariates
  vector[K] sd_X;      // Standard deviations of covariates
  for (k_i in 1:K) {
    mean_X[k_i] = mean(X[, k_i]);
    sd_X[k_i] = sd(X[, k_i]);
    // Handle constant covariates
    if (sd_X[k_i] > 1e-6) {
      X_std[, k_i] = (X[, k_i] - mean_X[k_i]) / sd_X[k_i];
    } else {
      X_std[, k_i] = rep_vector(0.0, C);
    }
  }
}

parameters {
  real alpha;                     // Intercept for latent logistic model
  vector[K] beta;                 // Coefficients for covariates
  real tau;                       // Treatment effect on latent dissatisfaction
  real<lower=0, upper=0.5> epsilon_0; // False Positive Rate: P(R=1 | D=0)
  real<lower=0, upper=0.5> epsilon_1; // False Negative Rate: P(R=0 | D=1)
}

model {
  // --- Priors ---
  alpha ~ normal(mean_alpha, sd_alpha);
  beta ~ normal(mean_beta, sd_beta);
  tau ~ normal(tau_mean, tau_sd);

  epsilon_0 ~ beta(epsilon0_alpha, epsilon0_beta);
  epsilon_1 ~ beta(epsilon1_alpha, epsilon1_beta);

  // --- Likelihood (Mixture Model) ---
  if (run_estimation == 1) {
    // Linear predictor for the latent probability of dissatisfaction
    vector[C] theta_c = alpha + X_std * beta + tau * treat;

    // Log-likelihood for each call
    for (c in 1:C) {
      // Scenario 1: Call is truly dissatisfied (D_c = 1)
      // Prob of this scenario: logit^{-1}(theta_c[c])
      // Prob of data given scenario: k[c] ~ Binomial(N[c], 1 - epsilon_1)
      real log_prob_if_dissatisfied = bernoulli_logit_lpmf(1 | theta_c[c]) +
                                      binomial_lpmf(k[c] | N[c], 1 - epsilon_1);

      // Scenario 2: Call is truly satisfied (D_c = 0)
      // Prob of this scenario: logit^{-1}(-theta_c[c])
      // Prob of data given scenario: k[c] ~ Binomial(N[c], epsilon_0)
      real log_prob_if_satisfied = bernoulli_logit_lpmf(0 | theta_c[c]) +
                                    binomial_lpmf(k[c] | N[c], epsilon_0);

      // Add the marginal log-likelihood for call 'c' to the target
      // This sums the probabilities of the two mutually exclusive scenarios
      target += log_sum_exp(log_prob_if_dissatisfied, log_prob_if_satisfied);
    }
  }
}

generated quantities {
  // Individual-level posterior probability of dissatisfaction
  // P(D_c=1 | R_c, params)
  vector[C] prob_dissatisfied;

  // Average Treatment Effect (ATE) on the latent probability of dissatisfaction
  real eta;

  // Linear predictors
  vector[C] theta_c_treated = alpha + X_std * beta + tau * 1;
  vector[C] theta_c_control = alpha + X_std * beta;
  vector[C] theta_c_observed = alpha + X_std * beta + tau * treat;

  for (c in 1:C) {
    // Calculate P(D_c=1 | R_c, params) using Bayes' rule
    real log_prob_if_dissatisfied = bernoulli_logit_lpmf(1 | theta_c_observed[c]) +
                                    binomial_lpmf(k[c] | N[c], 1 - epsilon_1);
    real log_prob_if_satisfied = bernoulli_logit_lpmf(0 | theta_c_observed[c]) +
                                  binomial_lpmf(k[c] | N[c], epsilon_0);
    
    // The above probabilities need to be standardized due to the scaling of p(R_c)
    // prob = exp(A) / (exp(A) + exp(B)) = 1 / (1 + exp(B - A)) = inv_logit(A - B)
    prob_dissatisfied[c] = inv_logit(log_prob_if_dissatisfied - log_prob_if_satisfied);
  }

  // Calculate ATE on the probability of dissatisfaction
  eta = mean(inv_logit(theta_c_treated)) - mean(inv_logit(theta_c_control));
}

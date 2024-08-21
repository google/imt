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
  int N;                          // Number of observations
  int<lower = 0, upper = 1> y[N]; // Outcome (binary 0 or 1)
  int K;                          // Number of covariates
  matrix[N, K] X;                 // Model matrix (contains predictor values)
  real mean_alpha;
  real sd_alpha;
  vector[K] mean_beta;
  vector[K] sd_beta;
  vector[N] treat;      // Treatment indicator
  real tau_mean;        // Prior mean for treatment effect
  real<lower=0> tau_sd; // Prior standard deviation for treatment effect
  // Flag for running estimation (0: no, 1: yes)
  int<lower=0, upper=1> run_estimation;
}

transformed data {
  matrix[N, K] X_std;  // Standardized covariates
  vector[K] mean_X;    // Means of covariates
  vector[K] sd_X;      // Standard deviations of covariates
  for (k in 1:K) {
    mean_X[k] = mean(X[,k]);
    sd_X[k] = sd(X[,k]);
    X_std[,k] = (X[,k] - mean_X[k]) / sd_X[k];
  }
}

parameters {
  real alpha;         // Intercept parameter
  vector[K] beta;     // Regression parameters (coefficients for each covariate)
  real tau;         // treatment effect parameter
}

model {
  vector[N] theta;        // Linear predictor (stores calculated probabilities)
  // Priors
  alpha ~ normal(mean_alpha, sd_alpha);
  for (k in 1:K) {
    beta[k] ~ normal(mean_beta[k], sd_beta[k]);
  }
  tau ~ normal(tau_mean, tau_sd);
  theta = alpha + X_std * beta + tau * treat;     // Calculate linear predictor
  if (run_estimation==1) {
  // Likelihood with logit link function:
      y ~ bernoulli_logit(theta); // Models the outcome 'y' as following a
            //  Bernoulli distribution with probabilities determined
            // by the logit of 'theta'
  }
}

generated quantities {
  vector[N] y_sim;
  vector[N] y0;
  vector[N] y1;
  real eta;
  real mean_y_sim;
  for (i in 1:N) {
    y_sim[i] = bernoulli_logit_rng(alpha + tau * treat[i] +
      dot_product(X_std[i], beta));
    y0[i] = bernoulli_logit_rng(alpha + dot_product(X_std[i], beta));
    y1[i] = bernoulli_logit_rng(alpha + dot_product(X_std[i], beta) + tau);
    }
  eta = mean(y1) - mean(y0);
  mean_y_sim = mean(y_sim);
}

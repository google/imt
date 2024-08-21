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
  int N;                // Number of observations
  int K;                // Number of covariates
  int y[N];             // Outcome variable
  matrix[N, K] X;       // Model matrix
  vector[N] treat;      // Treatment indicator
  real tau_mean;        // Prior mean for treatment effect
  real<lower=0> tau_sd; // Prior standard deviation for treatment effect
}

transformed data {
  matrix[N, K] X_std;    // Standardized covariates
  vector[K] mean_X;      // Means of covariates
  vector[K] sd_X;        // Standard deviations of covariates
  for (k in 1:K) {
    mean_X[k] = mean(X[,k]);
    sd_X[k] = sd(X[,k]);
    X_std[,k] = (X[,k] - mean_X[k]) / sd_X[k];
  }
}

parameters {
  real alpha;         // Intercept
  vector[K] beta;     // Regression parameters
  real log_phi;       // Logarithm of overdispersion parameter
  real tau;           // Treatment effect
}

transformed parameters {
  real phi = exp(log_phi);  // Overdispersion parameter
}

model {
  vector[N] mu0;      // Linear predictor for control group
  vector[N] mu1;      // Linear predictor for treatment group
  tau ~ normal(tau_mean, tau_sd);
  log_phi ~ std_normal();
  log_phi ~ std_normal();
  alpha ~ std_normal();
  beta ~ std_normal();
  mu0 = exp(X * beta);  // Using log link
  mu1 = mu0 * (1 + tau);
  for (n in 1:N) {
    if (treat[n] == 1) {
      y[n] ~ neg_binomial_2(mu1[n], phi);
    } else {
      y[n] ~ neg_binomial_2(mu0[n], phi);
    }
  }
}

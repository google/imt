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
  // Number of observations
  int<lower=1> N;
  // Outcome variable
  vector[N] y;
  // Treatment indicator
  vector[N] treat;
  // Prior mean for treatment effect
  real eta_mean;
  // Prior standard deviation for treatment effect
  real<lower=0> eta_sd;
  // Flag for running estimation (0: no, 1: yes)
  int<lower=0, upper=1> run_estimation;
}

transformed data {
  // Mean of outcome
  real mean_y = mean(y);
  // Standard deviation of outcome
  real sd_y = sd(y);
  // Standardized data
  vector[N] y_std = (y - mean_y) / sd_y;
}

parameters {
  // standard deviation
  real<lower=0> sigma;
  // Treatment effect
  real eta;
  // Intercept
  real alpha;
}

model {
  // Priors
  alpha ~ std_normal();
  eta ~ normal(eta_mean, eta_sd);
  sigma ~ std_normal();

  // Likelihood
  if (run_estimation==1) {
    target += normal_lpdf(y_std | eta*treat + alpha, sigma);
  }
}

generated quantities {
  // Simulated outcome
  vector[N] y_sim;
  // Outcome if treated
  vector[N] y_one;
  // Outcome if untreated
  vector[N] y_zero;
  // Rescaled treatment effect
  real eta_rsc = eta*sd_y;

  for (i in 1:N) {
    y_sim[i] = normal_rng(eta*treat[i] + alpha , sigma) * sd_y + mean_y;
    y_one[i] = normal_rng(eta + alpha, sigma) * sd_y + mean_y;
    y_zero[i] = normal_rng(alpha, sigma) * sd_y + mean_y;
  }
}

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

// This is a simple Bayesian meta analysis model without covariates based on
// the clasic 8 schools example (Rubin 1981,
// and Gelman, Carlin, Stern, and Rubin 2003 Section 5.5).

data {
  int<lower=0> J;         // number of studies
  real y[J];              // estimated lift
  real<lower=0> sigma[J]; // standard error of lift estimates
  int<lower=0, upper=1> run_estimation;
  real<lower=0> sd_mu;
  real mean_mu;
}

parameters {
  real mu;                // population lift
  real<lower=0> tau;      // standard deviation in lift
  vector[J] eta;          // unscaled deviation from mu by study
}

transformed parameters {
  vector[J] theta = mu + tau * eta;        // study lift
}

model {
  mu ~ normal(mean_mu, sd_mu);
  tau ~ std_normal();
  eta ~ normal(0,sd_mu);
  if (run_estimation==1) {
    y ~ normal(theta, sigma);
  }
}

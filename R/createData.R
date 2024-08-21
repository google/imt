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

#' Converts a dataframe into a named list to provide data to a Stan model
#'
#' @param y The name of the outcome variable in the data frame
#' @param x A vector of names of all covariates to be used in the data frame
#' @param treatment The name of the treatment indicator variable in the data frame
#' @param data The data frame to be used
#' @param eta_mean The prior mean to be used for estimating the treatment effect
#' @param eta_sd The prior standard deviation to be used for estimating the treatment effect
#' @param run_estimation Whether to run the estimation, or merely draw data from the priors
#'
#' @return Returns \cr\code{stan_data} a named list providing the data for the stan model
#' @export
createData <-
  function(data,
           y,
           treatment,
           x = NULL,
           eta_mean = 0,
           eta_sd = 1,
           run_estimation = 1) {
    cleaned_data <- cleanData(
      data = data,
      y = y,
      treatment = treatment,
      x = x
    )
    # Stan data list
    stan_data <- list(
      N = cleaned_data$N,
      y = cleaned_data$Y,
      treat = cleaned_data$treat_vec,
      K = cleaned_data$K,
      X = cleaned_data$X,
      eta_mean = eta_mean,
      eta_sd = eta_sd,
      run_estimation = run_estimation
    )
    return(stan_data)
  }

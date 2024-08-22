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

#' Fits Stan model.
#' @param stan_data A named list providing the data for the Stan model.
#' @param model A character string specifying the model type. 
#'   Must be either "blm" (Bayesian linear model) or "bnb" (Bayesian negative binomial model).
#'   Defaults to "blm".
#' @param ... Additional options to be passed through to `rstan::sampling`.
#' @return The complete StanFit object from the fitted Stan model.
#' @export
#'
fit <- function(stan_data, model = "blm", ...) {
  # Check model type
  if (!model %in% c("blm", "bnb")) {
    stop("Invalid model type. Please choose either \"blm\" or \"bnb\".")
  }
  if (model == "blm") {
    message("Bayesian linear model")
    if (stan_data$K == 0) {
      model <- stanmodels$blmnox
    } else {
      model <- stanmodels$blm
    }
  } else {
    message("Bayesian negative binomial model")
    model <- stanmodels$NB
  }
  # Fit the Stan model using the provided data and optional arguments.
  out <- rstan::sampling(model, data = stan_data, ...)
  # Return the StanFit object.
  return(out)
}

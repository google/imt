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

##' Calculate credible interval from MCMC draws
##'
##' This function calculates a credible interval of the specified
##' width from a vector of MCMC draws.
##'
##' @param draws A numeric vector containing MCMC draws.
##' @param width A numeric value between 0 and 1 specifying the desired
##'     width of the credible interval.
##'
##' @return A named list containing three elements:
##'   - width: The specified width of the credible interval.
##'   - lower_bound: The lower bound of the credible interval.
##'   - upper_bound: The upper bound of the credible interval.
##'
##' @details The function calculates the lower and upper bounds of the
##'     credible interval using the quantile function based on the
##'     specified width.
##'
##' @examples
##' # Generate example draws
##' draws <- rnorm(1000)
##'
##' # Calculate 95% credible interval
##' credibleInterval(draws, width = 0.95)
##'
##' @export
##'
credibleInterval <- function(draws, width) {
  # Check argument validity
  if (!is.numeric(draws) || !is.vector(draws)) {
    stop("'draws' must be a numeric vector")
  }
  if (!is.numeric(width) || width < 0 || width > 1) {
    stop("'width' must be a numeric value between 0 and 1")
  }
  # Calculate lower and upper bounds
  lower_bound <- stats::quantile(draws, (1 - width) / 2)
  upper_bound <- stats::quantile(draws, 1 - (1 - width) / 2)
  # Return named list
  return(list(
    width = width,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  ))
}

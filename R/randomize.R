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

utils::globalVariables(c(
  "group", "variables", "std_diff", "balance", "xmin",
  "xmax", "balance"
))
#' @title Randomization Class for Treatment Assignment
#' @docType class
#' @description This class provides methods to randomly assign treatments
#' to a dataset while ensuring baseline covariate balance. It can
#' handle both simple and stratified randomization.
#' @export
#'
#' @field version The version of the \code{im} package used for randomization.
#' @field data The data frame with the assigned treatment.
#' @field seed The random seed used for reproducibility.
#' @field balance_summary A summary (or list of summaries) of the balance
#' assessment after randomization.
#' @field balance_plot A plot (or list of plots) of the balance assessment
#' after randomization.
#'
randomizer <- R6::R6Class(
  classname = "randomized",
  private = list(
    ..version = NULL,
    ..data = NULL,
    ..seed = NULL,
    ..balance_summary = NULL, # Can now hold a list for stratified randomization
    ..balance_plot = NULL # Can now hold a list for stratified randomization
  ),
  active = list(
    version = function() {
      return(private$..version)
    },
    data = function() {
      return(private$..data)
    },
    seed = function() {
      return(private$..seed)
    },
    balance_summary = function() {
      return(private$..balance_summary)
    },
    balance_plot = function() {
      return(imt::balancePlot(data = private$..balance_summary))
    }
  ),
  public = list(
    #' @description Initialize a new Randomizer object.
    #' @param data The input data frame.
    #' @param seed (Optional) An integer to set the random seed.
    #' @param max_attempts (Optional) Maximum number of randomization attempts.
    #' @param variables A vector of covariate names to check for balance.
    #' @param standard The desired level of baseline equivalence.
    #' Must be one of "Not Concerned",
    #' "Concerned", or "Very Concerned". Default is "Not Concerned".
    #' ("Not Concerned", "Concerned", or "Very Concerned").
    #' @param group_by (Optional) A character vector of column names to
    #' stratify randomization.
    #' @return A new \code{randomizer} object.
    initialize = function(data, variables, standard = "Not Concerned",
                          seed = NULL, max_attempts = 100, group_by = NULL) {
      private$..version <- packageVersion("imt")
      randomized <- randomize(
        data = data,
        seed = seed, max_attempts = max_attempts,
        variables = variables, standard = standard,
        pr_treated = 0.5, group_by = group_by
      )
      private$..data <- randomized$data
      private$..seed <- randomized$seed
      private$..balance_summary <- randomized$balance

      # Check for stratified randomization and create plots accordingly
      if (inherits(randomized$balance, "list")) {
        private$..balance_plot <- lapply(randomized$balance, imt::balancePlot)
      } else {
        private$..balance_plot <- imt::balancePlot(data = randomized$balance)
      }
      return(invisible(self))
    }
  )
)

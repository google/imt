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

#' @title Bayesian Negative Binomial Model Factory
#' @docType class
#' @export

negativeBinomial <- R6::R6Class(
  classname = "bnb",
  private = list(
    ..mcmc_checks = NULL,
    ..version = NULL,
    ..stanfit = NULL,
    ..stan_data = NULL,
    ..tau_draws = NULL,
    ..tau_mean = NULL,
    ..tau_sd = NULL,
    ..credible_interval = NULL
  ),
  active = list(
    #' @field version Package version used to fit the model
    version = function() {
      return(private$..version)
    },
    #' @field mcmChecks MCMC diagnostics
    mcmChecks = function() {
      return(private$..mcmc_checks)
    },
    #' @field credible_interval Credible interval for the treatment effect
    credible_interval = function() {
      return(private$..credible_interval)
    },
    #' @field tau_draws Posterior draws for the treatment effect
    tau_draws = function() {
      return(private$..tau_draws)
    }
  ),
  public = list(
    #' @description
    #' Create a new Bayesian Negative Binomial Model object.
    #'
    #' @param y Name of the outcome variable in the data frame
    #' @param x Vector of names of all covariates in the data frame
    #' @param treatment Name of the treatment indicator variable in the data frame
    #' @param data Data frame to be used
    #' @param tau_mean Prior mean for the treatment effect estimation
    #' @param tau_sd Prior standard deviation for the treatment effect estimation
    #' @param run_estimation Integer flag to control whether estimation is run (1) or not (0)
    #' @param seed Seed for Stan fitting
    #' @param ... Additional arguments for Stan
    #' @return invisible
    initialize = function(data, y, x, treatment, tau_mean, tau_sd,
                          run_estimation = 1, seed = 1982, ...) {
      private$..version <- packageVersion("imt")
      private$..tau_mean <- tau_mean
      private$..tau_sd <- tau_sd
      Y <- data |> pull({{ y }})
      if (any(Y < 0)) {
        stop("The outcome variable should be non-negative")
      }
      if (!is.integer(Y)) {
        stop("The outcome variable should Integers")
      }
      cleaned_data <- cleanData(
        data = data,
        y = y,
        treatment = treatment,
        x = x
      )
      # Fit model
      private$..stan_data <- list(
        N = cleaned_data$N,
        K = cleaned_data$K,
        X = cleaned_data$X,
        y = cleaned_data$Y,
        treat = cleaned_data$treat_vec,
        tau_mean = tau_mean,
        tau_sd = tau_sd
      )

      message("Fitting model to the data")
      private$..stanfit <- fit(private$..stan_data, model = "bnb", ...)
      private$..mcmc_checks <- mcmcChecks$new(
        fit = private$..stanfit,
        pars = c("tau", "phi")
      )
      private$..tau_draws <-
        as.data.frame(private$..stanfit, pars = "tau") |>
        dplyr::pull(tau)
      return(invisible())
    },
    #' @description
    #' Plot MCMC trace for the eta and sigma parameters.
    #' @param ... Additional arguments for Stan
    #' @return A ggplot object.
    tracePlot = function(...) {
      return(
        bayesplot::mcmc_trace(private$..stanfit,
          pars = c("tau", "phi"),
          facet_args = list(nrow = 2), ...
        )
      )
    },
    #' @description
    #' Calculate posterior probability of effect exceeding a threshold
    #'
    #' This function calculates the posterior probability of the effect
    #' being larger or smaller than a specified threshold.
    #'
    #' @param threshold A numeric value specifying the threshold.
    #' @param gt A logical value indicating whether to calculate the probability
    #'     of the effect being greater than (TRUE) or less than (FALSE)
    #'     the threshold.
    #'
    #' @return A character string summarizing the estimated probability.
    #'
    #' @details This function uses the private$..tau_draws internal variable
    #'     to obtain draws from the posterior distribution of the effect size.
    #'     Based on the specified arguments, the function calculates the
    #'     proportion of draws exceeding/falling below the threshold and
    #'     returns a formatted statement describing the estimated probability.
    #'
    posteriorProb = function(threshold = 0, gt = TRUE) {
      if (gt) {
        pr <- mean(private$..tau_draws > threshold)
        moreless <- "more"
      } else {
        pr <- mean(private$..tau_draws < threshold)
        moreless <- "less"
      }
      statement <- glue::glue(
        "Given the data, we estimate that the ",
        "probability that the effect is {moreless} ",
        "than {threshold} is {scales::percent(pr)}."
      )
      return(statement)
    },
    #' Calculate point estimate of the effect
    #'
    #' This R6 method calculates the point estimate of the effect size
    #' based on the posterior draws of the tau parameter.
    #'
    #' @param median Logical value. If TRUE (default), the median of
    #'     the tau draws is returned. If FALSE, the mean is returned.
    #'
    #' @return A numeric value representing the point estimate.
    #'
    #' @details This method uses the private$..tau_draws internal variable
    #'     which contains MCMC draws of the tau parameter representing the
    #'     effect size. Based on the specified median argument, the method
    #'     calculates and returns either the median or the mean of the draws.
    pointEstimate = function(median = TRUE) {
      if (median) {
        return(median(private$..tau_draws))
      } else {
        return(mean(private$..tau_draws))
      }
    },
    #' Calculates credible interval for the effect of the intervention
    #'
    #' This R6 method calculates and returns a formatted statement summarizing
    #' the credible interval of a specified width for the effect of the intervention.
    #'
    #' @param width Numeric value between 0 and 1 representing the desired
    #'     width of the credible interval (e.g., 0.95 for a 95% credible interval).
    #' @param round Integer value indicating the number of decimal places to round
    #'     the lower and upper bounds of the credible interval.
    #'
    #' @return A character string with the following information:
    #'   - The probability associated with the specified width
    #'   - The lower and upper bounds of the credible interval, rounded to the
    #'      specified number of decimal places
    #'
    #' @details This method uses the private$..tau_draws internal variable
    #'     containing MCMC draws of the tau parameter representing the effect
    #'     size. It calculates the credible interval, stores it internally, and
    #'     returns a formatted statement summarizing the findings.
    #'
    credibleInterval = function(width = 0.75, round = 0) {
      private$..credible_interval <- credibleInterval(
        draws = private$..tau_draws, width = width
      )
      statement <- glue::glue(
        "Given the data, we estimate that there is a ",
        "{scales::percent(width)} probability that the effect is between ",
        "{round(private$..credible_interval$lower_bound, round)} and ",
        "{round(private$..credible_interval$upper_bound, round)}."
      )
      return(statement)
    },
    #' @description
    #' Plots impact's prior and posterior distributions.
    #'
    #' For more details see [vizdraws::vizdraws()].
    #' @param ... other arguments passed to vizdraws.
    #' @return An interactive plot of the prior and posterior distributions.
    vizdraws = function(...) {
      prior <-
        glue::glue("N({private$..tau_mean}, {private$..tau_sd})")
      p <- vizdraws::vizdraws(
        prior = prior,
        posterior = private$..tau_draws, ...
      )
      return(p)
    }
  )
)

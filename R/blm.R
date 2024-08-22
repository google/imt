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

#' @title Bayesian Linear Model Factory
#' @docType class
#' @export
#' @field version im package version used to fit model
#' @field eta_draws Posterior draws for the treatment effect
#' @field mcmChecks MCMC diagnostics
#' @field credible_interval Credible interval for the treatment effect

blm <- R6::R6Class(
  classname = "blm",
  private = list(
    ..mcmc_checks = NULL, # a mcmc_checks object
    ..version = NULL,
    ..stanfit = NULL,
    ..eta_draws = NULL,
    ..stan_data = NULL,
    ..true_eta = NULL,
    ..credible_interval = NULL,
    ..eta_prior_mean = NULL,
    ..eta_prior_sd = NULL,
    ..y_sd = NULL
  ),
  active = list(
    #' @description Get the package version
    version = function() {
      return(private$..version)
    },
    #' @description Get the posterior draws for eta
    eta_draws = function() {
      return(private$..eta_draws)
    },
    #' @description Get the MCMC diagnostics
    mcmChecks = function() {
      return(private$..mcmc_checks)
    },
    #' @description Get the credible interval
    credible_interval = function() {
      return(private$..credible_interval)
    }
  ),
  public = list(
    #' @description
    #' Create a new Bayesian Linear Model object.
    #'
    #' @param data Data frame to be used
    #' @param y Name of the outcome variable in the data frame
    #' @param x Vector of names of all covariates in the data frame
    #' @param treatment Name of the treatment indicator variable in the data frame
    #' @param eta_mean Prior mean for the treatment effect estimation
    #' @param eta_sd Prior standard deviation for the treatment effect estimation
    #' @param generate_fake_data Flag for generating fake data
    #' @param seed Seed for Stan fitting
    #' @param ... Additional arguments for Stan
    #' @return invisible
    initialize = function(data, y, x, treatment, eta_mean, eta_sd,
                          generate_fake_data = 0, seed = 1982, ...) {
      private$..version <- packageVersion("imt")
      stan_data <- createData(
        data = data, y = y, treatment = treatment,
        x = x, eta_mean = eta_mean, eta_sd = eta_sd,
        run_estimation = 1 - generate_fake_data
      )
      private$..y_sd <- sd(stan_data$y)
      private$..eta_prior_sd <- eta_sd * private$..y_sd
      private$..eta_prior_mean <- eta_mean * private$..y_sd
      # Fit model
      if (generate_fake_data == 1) {
        message("Generating fake data")
        # Generate fake data
        sim_out <- fit(stan_data, seed = seed, iter = 100)
        fake_data_matrix <- as.data.frame(sim_out, pars = "y_sim")
        draw <- sample.int(100, 1)
        y_sim <- as.numeric(fake_data_matrix[draw, ])
        analysis_df <- data |> dplyr::mutate(
          {{ y }} := y_sim,
          {{ treatment }} := as.logical({{ treatment }} == 1)
        )
        stan_data <- createData(
          data = analysis_df, y = y, treatment = treatment,
          x = x, eta_mean = eta_mean, eta_sd = eta_sd, run_estimation = 1
        )
        private$..true_eta <- as.matrix(sim_out, pars = "eta_rsc")[draw]
      } # end of fake data generating code
      private$..stan_data <- stan_data
      message("Fitting model to the data")
      private$..stanfit <- fit(stan_data, ...)
      private$..mcmc_checks <- mcmcChecks$new(
        fit = private$..stanfit,
        pars = c("eta", "sigma")
      )
      private$..eta_draws <-
        as.data.frame(private$..stanfit, pars = "eta_rsc") |>
        dplyr::pull(eta_rsc)
      return(invisible())
    },
    #' @description
    #' This method compares the empirical distribution of the data 'y'
    #' to the distributions of simulated/replicated data 'yrep' from the
    #' posterior predictive distribution. This is done by creating a density
    #' overlay plot using the `bayesplot::ppc_dens_overlay` function.
    #' @param n Number of posterior draws to use for the overlay
    #' @return ggplot2 visualization
    ppcDensOverlay = function(n = 50) {
      y_sim <- as.data.frame(private$..stanfit, pars = "y_sim")
      return(
        bayesplot::ppc_dens_overlay(
          y = private$..stan_data$y,
          yrep = as.matrix(y_sim[sample.int(nrow(y_sim), n), ])
        )
      )
    },
    #' @description
    #' Plot MCMC trace for the eta and sigma parameters.
    #' @param ... Additional arguments for Stan
    #' @return A ggplot object.
    tracePlot = function(...) {
      return(
        bayesplot::mcmc_trace(private$..stanfit,
          pars = c("eta", "sigma"),
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
    #' @details This function uses the private$..eta_draws internal variable
    #'     to obtain draws from the posterior distribution of the effect size.
    #'     Based on the specified arguments, the function calculates the
    #'     proportion of draws exceeding/falling below the threshold and
    #'     returns a formatted statement describing the estimated probability.
    #'
    posteriorProb = function(threshold = 0, gt = TRUE) {
      if (gt) {
        pr <- mean(private$..eta_draws > threshold)
        moreless <- "more"
      } else {
        pr <- mean(private$..eta_draws < threshold)
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
    #' based on the posterior draws of the eta parameter.
    #'
    #' @param median Logical value. If TRUE (default), the median of
    #'     the eta draws is returned. If FALSE, the mean is returned.
    #'
    #' @return A numeric value representing the point estimate.
    #'
    #' @details This method uses the private$..eta_draws internal variable
    #'     which contains MCMC draws of the eta parameter representing the
    #'     effect size. Based on the specified median argument, the method
    #'     calculates and returns either the median or the mean of the draws.
    pointEstimate = function(median = TRUE) {
      if (median) {
        return(median(private$..eta_draws))
      } else {
        return(mean(private$..eta_draws))
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
    #' @details This method uses the private$..eta_draws internal variable
    #'     containing MCMC draws of the eta parameter representing the effect
    #'     size. It calculates the credible interval, stores it internally, and
    #'     returns a formatted statement summarizing the findings.
    #'
    credibleInterval = function(width = 0.75, round = 0) {
      private$..credible_interval <- credibleInterval(
        draws = private$..eta_draws, width
      )
      statement <- glue::glue(
        "Given the data, we estimate that there is a ",
        "{scales::percent(width)} probability that the effect is between ",
        "{round(private$..credible_interval$lower_bound, round)} and ",
        "{round(private$..credible_interval$upper_bound, round)}."
      )
      return(statement)
    },
    #' Calculate and format probability statement based on prior
    #'
    #' This method calculates the probability that the effect is greater than or
    #' less than a given threshold based on the prior distribution of the effect.
    #' The results are formatted into a statement suitable for reporting.
    #'
    #' @param threshold Numerical threshold for comparison.
    #' @param gt Logical indicating whether to calculate probability greater than
    #'     or less than the threshold.
    #' @return A character string containing the formatted probability statement.
    priorProb = function(threshold = 0, gt = TRUE) {
      pr <- pnorm(
        q = threshold, mean = private$..eta_prior_mean,
        sd = private$..eta_prior_sd, lower.tail = !gt
      )
      if (gt) {
        moreless <- "more"
      } else {
        moreless <- "less"
      }
      statement <- glue::glue(
        "Our prior is that the ",
        "probability that the effect is {moreless} ",
        "than {threshold} is {scales::percent(pr)}."
      )
      return(statement)
    },
    #' Calculates the prior for the effedt of the intervention based on
    #' the hyperpriors.
    #'
    #' This method computes and formats a statement about the probability interval
    #' of the effect based on the prior distribution.
    #'
    #' @param width Desired probability width of the interval (default: 0.75).
    #' @param round Number of decimal places for rounding the bounds (default: 0).
    #' @return A character string containing the formatted probability interval
    #'     statement.
    #' @details
    #' This method calculates the lower and upper bounds of the interval based
    #' on the specified probability width and the prior distribution of the effect.
    #' It then formats the results into a clear and informative statement.
    #'
    #' Note that the probability is checked to be within valid range (0-1) with
    #' consideration of machine precision using .Machine$double.eps.
    priorInterval = function(width = 0.75, round = 0) {
      if (width < 0 + .Machine$double.eps ||
        width > 1 - .Machine$double.eps) {
        stop("Probability (width) must be between 0 and 1.")
      }
      # Calculate half probability for each tail
      tail_prob <- (1 - width) / 2
      # Find z-scores corresponding to tail probabilities
      z_lower <- qnorm(tail_prob)
      z_upper <- qnorm(1 - tail_prob)
      # Calculate thresholds based on mean and standard deviation
      lower_bound <- private$..eta_prior_mean +
        z_lower * private$..eta_prior_sd
      upper_bound <- private$..eta_prior_mean +
        z_upper * private$..eta_prior_sd
      statement <- glue::glue(
        "Our prior is that there is a ",
        "{scales::percent(width)} probability that the effect is between ",
        "{round(lower_bound, round)} and ",
        "{round(upper_bound, round)}."
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
        glue::glue("N({private$..eta_prior_mean}, {private$..eta_prior_sd})")
      p <- vizdraws::vizdraws(
        prior = prior,
        posterior = private$..eta_draws, ...
      )
      return(p)
    }
  )
)

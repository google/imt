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

#' @title Bayesian Logit Model Factory
#' @docType class
#' @export
#' @description A class for creating and managing Bayesian Logit Models
#' @field version im package version used to fit model
#' @field tau_draws Posterior draws for the treatment effect
#' @field mcmChecks MCMC diagnostics
#' @field credible_interval Credible interval for the treatment effect
#' @field prior_eta Prior distribution for eta
#' @field prior_tau Prior distribution for tau
#' @field prior_mean_y Prior distribution for mean y
#' @field eta_draws Posterior draws for eta
#' @field predict_list List of predictions

logit <- R6::R6Class(
  classname = "logit",
  private = list(
    ..mcmc_checks = NULL, # a mcmc_checks object
    ..version = NULL,
    ..stanfit = NULL,
    ..tau_draws = NULL,
    ..stan_data = NULL,
    ..true_tau = NULL,
    ..true_eta = NULL,
    ..credible_interval = NULL,
    ..tau_prior_mean = NULL,
    ..tau_prior_sd = NULL,
    ..eta_draws = NULL,
    ..alpha_draws = NULL,
    ..beta_draws = NULL,
    ..prior_eta = NULL,
    ..prior_tau = NULL,
    ..prior_mean_y = NULL,
    ..var_cols = NULL,
    ..treatment = NULL,
    ..predict_list = NULL,
    ..predictions = NULL
  ),
  active = list(
    #' @description Get the package version
    version = function() {
      return(private$..version)
    },
    #' @description Get the posterior draws for tau
    tau_draws = function() {
      return(private$..tau_draws)
    },
    #' @description Get the MCMC diagnostics
    mcmChecks = function() {
      return(private$..mcmc_checks)
    },
    #' @description Get the credible interval
    credible_interval = function() {
      return(private$..credible_interval)
    },
    #' @description Get the prior for eta
    prior_eta = function() {
      return(private$..prior_eta)
    },
    #' @description Get the prior for tau
    prior_tau = function() {
      return(private$..prior_tau)
    },
    #' @description Get the prior for mean y
    prior_mean_y = function() {
      return(private$..prior_mean_y)
    },
    #' @description Get the posterior draws for eta
    eta_draws = function() {
      return(private$..eta_draws)
    },
    #' @description Get the list of predictions
    predict_list = function() {
      return(private$..predict_list)
    }
  ),
  public = list(
    #' @description
    #' Create a new Bayesian Logit Model object.
    #'
    #' @param data Data frame to be used
    #' @param y Name of the outcome variable in the data frame
    #' @param x Vector of names of all covariates in the data frame
    #' @param treatment Name of the treatment indicator variable in the data frame
    #' @param mean_alpha Prior mean for alpha
    #' @param sd_alpha Prior standard deviation for alpha
    #' @param mean_beta Prior mean for beta
    #' @param sd_beta Prior standard deviation for beta
    #' @param tau_mean Prior mean for the treatment effect estimation
    #' @param tau_sd Prior standard deviation for the treatment effect estimation
    #' @param seed Seed for Stan fitting
    #' @param fit Flag for fitting the data to the model or not
    #' @param ... Additional arguments for Stan
    #' @return invisible
    initialize = function(data, y, x, treatment,
                          mean_alpha,
                          sd_alpha,
                          mean_beta,
                          sd_beta,
                          tau_mean,
                          tau_sd,
                          seed = 1982,
                          fit = TRUE,
                          ...) {
      private$..version <- packageVersion("imt")
      private$..var_cols <- x
      private$..treatment <- treatment
      cleaned_data <- cleanData(
        data = data,
        y = y,
        treatment = treatment,
        x = x,
        binary = TRUE
      )
      stan_data <- list(
        N = cleaned_data$N,
        y = cleaned_data$Y,
        K = cleaned_data$K,
        X = cleaned_data$X,
        mean_alpha = mean_alpha,
        sd_alpha = sd_alpha,
        mean_beta = mean_beta,
        sd_beta = sd_beta,
        treat = cleaned_data$treat_vec,
        tau_mean = tau_mean,
        tau_sd = tau_sd,
        run_estimation = 0
      )
      private$..tau_prior_sd <- tau_sd
      private$..tau_prior_mean <- tau_mean
      private$..stan_data <- stan_data
      # Draw from the prior
      sim_out <- rstan::sampling(stanmodels$logit,
        data = private$..stan_data
      )
      stan_data$run_estimation <- 1
      private$..prior_eta <- as.matrix(sim_out, pars = "eta")[, 1]
      private$..prior_tau <- as.matrix(sim_out, pars = "tau")[, 1]
      private$..prior_mean_y <- as.matrix(sim_out, pars = "mean_y_sim")
      # Fit model
      if (fit) {
        message("Fitting model to the data")
        private$..stanfit <- rstan::sampling(stanmodels$logit,
          data = stan_data, ...
        )
        private$..mcmc_checks <- mcmcChecks$new(
          fit = private$..stanfit,
          pars = "tau"
        )

        private$..tau_draws <-
          as.data.frame(private$..stanfit, pars = "tau") |>
          dplyr::pull(tau)

        private$..eta_draws <-
          as.data.frame(private$..stanfit, pars = "eta") |>
          dplyr::pull(eta)

        private$..alpha_draws <-
          as.data.frame(private$..stanfit, pars = "alpha") |>
          dplyr::pull(alpha)

        private$..beta_draws <-
          as.data.frame(private$..stanfit, pars = "beta")
      }
      return(invisible())
    },

    #' @description
    #' Plot MCMC trace for the eta and sigma parameters.
    #' @param ... Additional arguments for Stan
    #' @return A ggplot object.
    tracePlot = function(...) {
      return(
        bayesplot::mcmc_trace(private$..stanfit,
          pars = "tau", ...
        )
      )
    },

    #' @description
    #' Calculates the posterior of an effect being greater than, less than,
    #' or within a range defined by thresholds.
    #'
    #'
    #' @param x A numeric vector containing draws from the posterior
    # ' distribution of the effect size.
    #' @param a Optional. Lower bound for the threshold.
    #' @param b Optional. Upper bound for the threshold.
    #' @param prior Logical. If TRUE, calculates probabilities based on
    #' the prior distribution.
    #'        If FALSE (default), uses the posterior distribution.
    #'
    #' @return  A character string summarizing the estimated probability
    #'
    #'
    calcProb = function(a = 0, b = NULL, prior = FALSE) {
      # Input validation (same as before)
      if (is.null(a) && is.null(b)) {
        stop("Either 'a' or 'b' must be provided.")
      }
      if (!is.null(a) && !is.null(b) && b <= a) {
        stop("'b' must be greater than 'a'.")
      }
      if (prior) {
        x <- private$..prior_eta
        txt <- "Our prior is "
      } else {
        x <- private$..eta_draws
        txt <- "Given the data, we estimate "
      }
      if (!is.null(a) && is.null(b)) {
        p <- scales::percent(mean(x > a))
        statement <- glue::glue(
          "{txt} that the ",
          "probability that the effect is more than {a}",
          " percentage points is {p}."
        )
      } else if (is.null(a) && !is.null(b)) {
        p <- mean(x < b)
        statement <- glue::glue(
          "{txt} that the probability that the",
          " effect is less than {b} percentage points is {p}."
        )
      } else { # both 'a' and 'b' are present
        p <- mean(x > a & x < b)
        statement <- glue::glue(
          "{txt} that the probability that the effect",
          " is between {a} and {b} percentage points is {p}."
        )
      }
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
        return(median(private$..eta_draws) * 100)
      } else {
        return(mean(private$..eta_draws) * 100)
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
        "{round(private$..credible_interval$lower_bound * 100, round)} and ",
        "{round(private$..credible_interval$upper_bound * 100, round)} ",
        "percentage points."
      )
      return(statement)
    },
    #' @description
    #' Plots impact's prior and posterior distributions.
    #'
    #' @param tau Logical. If TRUE, plot tau instead of eta
    #' @param ... other arguments passed to vizdraws.
    #' @return An interactive plot of the prior and posterior distributions.
    vizdraws = function(tau = FALSE, ...) {
      if (tau) {
        p <- vizdraws::vizdraws(
          prior = private$..prior_tau,
          posterior = private$..tau_draws, ...
        )
      } else {
        p <- vizdraws::vizdraws(
          prior = private$..prior_eta * 100,
          posterior = private$..eta_draws * 100, ...
        )
      }
      return(p)
    },
    #' @description
    #' Plots lollipop chart for the prior and posterior of the impact being
    #' greater or less than a threshold.
    #'
    #' For more details see [vizdraws::lollipops()].
    #' @param threshold cutoff used to calculate the probability. Defaults to
    #' zero percent points
    #' @param ... other arguments passed to vizdraws.
    #' @return A lollipop chart with the prior and posterior probability of
    #' the impact being above or below a threshold.
    lollipop = function(threshold = 0, ...) {
      data <- data.frame(
        Name = "Impact",
        Prior = mean(private$..prior_eta * 100 > threshold),
        Posterior = mean(private$..eta_draws * 100 > threshold)
      )
      p <- vizdraws::lollipops(data, ...)
      return(p)
    },
    #' @description
    #' Plots draws from the prior distribution of the outcome, tau, and impact
    #' in percentage points.
    plotPrior = function() {
      mean_y <- ggplot2::ggplot(
        data = tibble::tibble(draws = private$..prior_mean_y),
        ggplot2::aes(x = draws)
      ) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::scale_x_continuous(labels = scales::percent) +
        ggplot2::annotate("text",
          x = quantile(private$..prior_mean_y, prob = 0.9),
          y = 200, label = paste0(
            "mean: ",
            round(mean(private$..prior_mean_y * 100), 1)
          ),
          hjust = 0
        ) +
        ggplot2::xlab(expression(bar(y))) +
        ggplot2::ylab("N draws") +
        ggplot2::theme_minimal()

      tau <- ggplot2::ggplot(
        data = tibble::tibble(draws = private$..prior_tau),
        ggplot2::aes(x = draws)
      ) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::xlab(expression(tau)) +
        ggplot2::ylab("N draws") +
        ggplot2::theme_minimal()

      eta <- ggplot2::ggplot(
        data = tibble::tibble(draws = private$..prior_eta * 100),
        ggplot2::aes(x = draws)
      ) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::annotate("text",
          x = quantile(private$..prior_eta * 100, prob = 0.9),
          y = 200, label = paste0(
            "mean: ",
            round(mean(private$..prior_eta) * 100, 1)
          ),
          hjust = 0
        ) +
        ggplot2::xlab("Impact in percentage points") +
        ggplot2::ylab("N draws") +
        ggplot2::theme_minimal()

      plots <- ggpubr::ggarrange(eta,
        ggpubr::ggarrange(tau, mean_y, labels = c("tau", "Outcome"), ncol = 2),
        labels = "Impact", nrow = 2
      )

      plots <- ggpubr::annotate_figure(plots,
        top = ggpubr::text_grob("Draws from prior distributions",
          face = "bold", size = 14
        )
      )
      return(plots)
    },

    #' @description
    #' Predict new data
    #'
    #' @param new_data Data frame to be predicted
    #' @param name Group name of the prediction
    #' @param M Number of posterior draws to sample from
    #' @param ... Additional arguments
    #' @return invisible(self)
    predict = function(new_data, name = NULL, M = NULL, ...) {
      # Check new_data contains correct columns
      cols <- c(private$..var_cols, private$..treatment)
      if (!all(cols %in% colnames(new_data))) {
        missing_cols <- cols[which(!(cols %in% colnames(new_data)))]
        statement <- glue::glue(
          "{missing_cols} is missing."
        )
        stop(statement)
      }

      # Check number of posterior draws is not out of range
      if (is.null(M)) {
        M <- length(private$..alpha_draws)
      }
      if (M > length(private$..alpha_draws)) {
        M <- length(private$..alpha_draws)
        warning(
          "Number of posterior draws can not exceed ", M, ". ",
          "Setting number of posterior draws to ", M, ". "
        )
      }

      # Create prediction list if NULL
      if (is.null(private$..predictions)) {
        private$..predictions <- list()
        private$..predict_list <- list()
      }
      # Create group name if NULL
      if (is.null(name)) {
        name <- paste0("pred", length(private$..predictions) + 1)
        warning(
          "No name was supplied, assigning predictions to ", name, "."
        )
      }

      # Transform data
      X_pred <- new_data[, private$..var_cols]
      mean_X <- colMeans((private$..stan_data)$X)
      sd_X <- apply((private$..stan_data)$X, 2, sd)
      K <- length(private$..var_cols)
      for (k in 1:K) {
        X_pred[, k] <- (X_pred[, k] - mean_X[k]) / sd_X[k]
      }

      # Sample posterior predictives
      N <- nrow(X_pred)
      y_sim <- matrix(NA, nrow = N, ncol = M)
      treat <- new_data |> dplyr::pull(private$..treatment)
      y_sim <- purrr::pmap(
        .l = list(
          alpha = private$..alpha_draws,
          tau = private$..tau_draws,
          beta = purrr::array_branch(private$..beta_draws, 1)
        ),
        .f = logitRng,
        X = X_pred, treat = treat, N = N
      )
      y_sim <- t(matrix(unlist(y_sim), ncol = N, byrow = TRUE))

      private$..predict_list[[length(private$..predict_list) + 1]] <- name
      private$..predictions[[name]] <- y_sim
      return(invisible(self))
    },

    #' @description
    #' Get posterior predictive draws
    #' @param name Group name of the prediction
    #' @param ... Additional arguments (not used)
    #' @return Matrix of posterior predictive draws
    getPred = function(name = NULL, ...) {
      # Check if predictions exists
      if (is.null(private$..predictions)) {
        stop("No predictions in the object.")
      }
      # If name is NULL, supply the last predicted item.
      if (is.null(name)) {
        name <- private$..predict_list[[length(private$..predict_list)]]
      }

      # If name is not in the list, throw an error
      if (!(name %in% private$..predict_list)) {
        stop(name, " is not in the prediction list.")
      }

      return(private$..predictions[[name]])
    },

    #' @description
    #' Get point estimate, credible interval and prob summary of predictive draws
    #'
    #' @param name Optional. Group name of the prediction
    #' @param subgroup Optional. A boolean vector to get summary on the conditional group average
    #' @param median Optional. Logical value for using median or mean
    #' @param width Optional. Numeric value for credible interval width
    #' @param round Optional. Integer value for rounding
    #' @param a Optional. Lower bound threshold
    #' @param b Optional. Upper bound threshold
    #' @param ... Additional arguments
    #' @return A character string with summary information
    predSummary = function(name = NULL,
                           subgroup = NULL,
                           median = TRUE,
                           width = 0.75,
                           round = 0,
                           a = NULL,
                           b = NULL,
                           ...) {
      # Check if predictions exists
      if (is.null(private$..predictions)) {
        stop("No predictions in the object.")
      }
      # If name is NULL, supply the last predicted item.
      if (is.null(name)) {
        name <- names(private$..predictions)[length(private$..predictions)]
      }

      # If name is not in the list, throw an error
      if (!(name %in% names(private$..predictions))) {
        stop(name, " is not in the prediction list.")
      }

      # Validate subgroup
      if (!is.null(subgroup)) {
        validate_logical_vector(subgroup, nrow(private$..predictions[[name]]))
      } else {
        subgroup <- rep(TRUE, nrow(private$..predictions[[name]]))
      }

      # Get posterior draws of group average
      mean_y_draws <- colMeans((private$..predictions[[name]])[subgroup, ])
      # Get point estimate of posterior draws
      point_estimate_y <- pointEstimate(x = mean_y_draws, median = median)
      pe_statement <- glue::glue(
        "Given the data, we estimate that for group: {name}, ",
        "the point estimate of the group average is ",
        "{round(point_estimate_y, round)}%."
      )

      # Get credible interval
      credible_interval <- credibleInterval(mean_y_draws, width)
      ci_statement <- glue::glue(
        "With {scales::percent(width)} probability, ",
        "the point estimate is between ",
        "{round(credible_interval$lower_bound * 100, round)} and ",
        "{round(credible_interval$upper_bound * 100, round)} ",
        "percentage points."
      )

      prob_statement <- calcProb(mean_y_draws, a, b,
        group_name = "group average"
      )

      statement <- paste(
        pe_statement,
        ci_statement,
        prob_statement
      )
      return(statement)
    },

    #' @description
    #' Compare the average of the posterior draws of two groups
    #'
    #' @param name1 Group name of the first prediction to be compared
    #' @param name2 Group name of the second prediction to be compared
    #' @param subgroup1 Optional. A boolean vector for the first group
    #' @param subgroup2 Optional. A boolean vector for the second group
    #' @param median Optional. Logical value for using median or mean
    #' @param width Optional. Numeric value for credible interval width
    #' @param round Optional. Integer value for rounding
    #' @param a Optional. Lower bound threshold
    #' @param b Optional. Upper bound threshold
    #' @param ... Additional arguments
    #' @return A character string with comparison summary
    predCompare = function(name1, name2,
                           subgroup1 = NULL,
                           subgroup2 = NULL,
                           median = TRUE,
                           width = 0.75,
                           round = 0,
                           a = NULL,
                           b = NULL,
                           ...) {
      # Check if predictions exists
      if (is.null(private$..predictions)) {
        stop("No predictions in the object.")
      }
      if (length(private$..predictions) < 2) {
        stop("Less than two groups were predicted.")
      }
      # If name is NULL, supply the last predicted item.
      if (is.null(name1) || is.null(name2)) {
        stop("Group names need to be provided.")
      }
      # If name is not in the list, throw an error
      if (!(name1 %in% names(private$..predictions))) {
        stop(name1, " is not in the prediction list.")
      }
      if (!(name2 %in% names(private$..predictions))) {
        stop(name2, " is not in the prediction list.")
      }

      # Validate subgroup1
      if (!is.null(subgroup1)) {
        validate_logical_vector(subgroup1, nrow(private$..predictions[[name1]]))
      } else {
        subgroup1 <- rep(TRUE, nrow(private$..predictions[[name1]]))
      }

      # Validate subgroup2
      if (!is.null(subgroup2)) {
        validate_logical_vector(subgroup2, nrow(private$..predictions[[name2]]))
      } else {
        subgroup2 <- rep(TRUE, nrow(private$..predictions[[name2]]))
      }

      # Get posterior draws of group average
      eta_draws <- colMeans((private$..predictions[[name1]])[subgroup1, ]) -
        colMeans((private$..predictions[[name2]])[subgroup2, ])
      # Get point estimate of posterior draws
      point_estimate_eta <- pointEstimate(x = eta_draws, median = median)
      pe_statement <- glue::glue(
        "Given the data, we estimate that ",
        "the point estimate of the group difference is ",
        "{round(point_estimate_eta, round)}%."
      )

      # Get credible interval
      credible_interval <- credibleInterval(eta_draws, width)
      ci_statement <- glue::glue(
        "With {scales::percent(width)} probability, ",
        "the point estimate is between ",
        "{round(credible_interval$lower_bound * 100, round)} and ",
        "{round(credible_interval$upper_bound * 100, round)} ",
        "percentage points."
      )

      # Calculate prob
      prob_statement <- calcProb(eta_draws, a, b,
        group_name = "group difference"
      )

      statement <- paste(
        pe_statement,
        ci_statement,
        prob_statement
      )
      return(statement)
    }
  )
)

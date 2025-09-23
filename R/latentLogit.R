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

#' @title Bayesian Latent Logit Model for Measurement Error
#' @docType class
#' @export
#' @description
#' A class for creating and managing Bayesian Logit Models that account for
#' measurement error in the binary outcome.
#'
#' @field version im package version used to fit model
#' @field tau_draws Posterior draws for the treatment effect coefficient (tau)
#' @field eta_draws Posterior draws for the Average Treatment Effect (ATE)
#'                  on the latent probability of dissatisfaction.
#' @field epsilon0_draws Posterior draws for the False Positive Rate (epsilon_0)
#' @field epsilon1_draws Posterior draws for the False Negative Rate (epsilon_1)
#' @field prob_dissatisfied_draws A matrix (draws x calls) of posterior
#'                                probabilities for each call being dissatisfied.
#' @field mcmChecks MCMC diagnostics
#' @field credible_interval Credible interval for the treatment effect (ATE)
#' @field prior_eta Prior distribution for eta
#' @field prior_tau Prior distribution for tau

latentLogit <- R6::R6Class(
  classname = "latentLogit",
  private = list(
    ..mcmc_checks = NULL,
    ..version = NULL,
    ..stanfit = NULL,
    ..tau_draws = NULL,
    ..eta_draws = NULL,
    ..epsilon0_draws = NULL,
    ..epsilon1_draws = NULL,
    ..prob_dissatisfied_draws = NULL,
    ..stan_data = NULL,
    ..credible_interval = NULL,
    ..tau_prior_mean = NULL,
    ..tau_prior_sd = NULL,
    ..prior_eta = NULL,
    ..prior_tau = NULL,
    ..var_cols = NULL,
    ..treatment_col = NULL,
    ..call_id_col = NULL,
    ..call_ids = NULL # stores the ordered list of call_ids
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
    #' @description Get the posterior draws for eta (ATE)
    eta_draws = function() {
      return(private$..eta_draws)
    },
    #' @description Get the posterior draws for epsilon_0 (False Positive Rate)
    epsilon0_draws = function() {
      return(private$..epsilon0_draws)
    },
    #' @description Get the posterior draws for epsilon_1 (False Negative Rate)
    epsilon1_draws = function() {
      return(private$..epsilon1_draws)
    },
    #' @description Get the posterior draws for P(Dissatisfied) for each call
    prob_dissatisfied_draws = function() {
      return(private$..prob_dissatisfied_draws)
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
    }
  ),
  public = list(
    #' @description
    #' Create a new Bayesian Latent Logit Model object.
    #'
    #' @param data Data frame in LONG format (one row per AI rating)
    #' @param call_id Name of the column identifying the unit/call
    #' @param y_rating Name of the outcome variable (the 0/1 AI rating)
    #' @param x_covariates Vector of names of all covariates
    #' @param treatment Name of the treatment indicator variable
    #' @param mean_alpha Prior mean for alpha (latent intercept)
    #' @param sd_alpha Prior standard deviation for alpha
    #' @param mean_beta Prior mean for beta (covariates)
    #' @param sd_beta Prior standard deviation for beta
    #' @param tau_mean Prior mean for the treatment effect (tau)
    #' @param tau_sd Prior standard deviation for the treatment effect (tau)
    #' @param epsilon0_alpha Beta prior alpha for False Positive Rate (default 1)
    #' @param epsilon0_beta Beta prior beta for False Positive Rate (default 1)
    #' @param epsilon1_alpha Beta prior alpha for False Negative Rate (default 1)
    #' @param epsilon1_beta Beta prior beta for False Negative Rate (default 1)
    #' @param seed Seed for Stan fitting
    #' @param fit Flag for fitting the data to the model or not
    #' @param ... Additional arguments for rstan::sampling (e.g., chains, iter)
    #' @return invisible
    initialize = function(data,
                          call_id,
                          y_rating,
                          x_covariates,
                          treatment,
                          mean_alpha = -3,
                          sd_alpha = 2,
                          mean_beta = 0,
                          sd_beta = 1,
                          tau_mean = 0,
                          tau_sd = 0.5,
                          epsilon0_alpha = 1,
                          epsilon0_beta = 1,
                          epsilon1_alpha = 1,
                          epsilon1_beta = 1,
                          seed = 1997,
                          restriction = FALSE,
                          fit = TRUE,
                          ...) {
      
      # Store variable names
      private$..version <- packageVersion("imt") # Or your package name
      private$..var_cols <- x_covariates
      private$..treatment_col <- treatment
      private$..call_id_col <- call_id
      
      # Aggregate Data
      # The Stan model needs one row per call, with k and N
      message("Aggregating data by call_id...")
      agg_data <- data |>
        dplyr::group_by(!!dplyr::sym(call_id)) |>
        dplyr::summarize(
          k = sum(!!dplyr::sym(y_rating), na.rm = TRUE),
          N = dplyr::n(),
          # Assumes covariates and treatment are constant for a given call_id
          dplyr::across(
            c(all_of(x_covariates), !!dplyr::sym(treatment)),
            dplyr::first
          )
        ) |>
        dplyr::ungroup()
      
      # Store call_ids in order for later mapping
      private$..call_ids <- agg_data[[call_id]]
      
      # Prepare Stan Data List
      if (length(x_covariates) == 1) {
        X_matrix <- as.matrix(agg_data[, x_covariates])
      } else {
        X_matrix <- agg_data[, x_covariates]
      }
      
      stan_data <- list(
        C = nrow(agg_data),
        K = length(x_covariates),
        X = X_matrix,
        treat = agg_data[[treatment]],
        k = agg_data$k,
        N = agg_data$N,
        mean_alpha = mean_alpha,
        sd_alpha = sd_alpha,
        mean_beta = rep(mean_beta, length(x_covariates)),
        sd_beta = rep(sd_beta, length(x_covariates)),
        tau_mean = tau_mean,
        tau_sd = tau_sd,
        epsilon0_alpha = epsilon0_alpha,
        epsilon0_beta = epsilon0_beta,
        epsilon1_alpha = epsilon1_alpha,
        epsilon1_beta = epsilon1_beta,
        run_estimation = 0 # Start with 0 for prior simulation
      )
      
      private$..tau_prior_sd <- tau_sd
      private$..tau_prior_mean <- tau_mean
      private$..stan_data <- stan_data
      
      # Draw from the prior
      message("Drawing from prior distributions...")
      if (restriction) {
        sim_out <- rstan::sampling(
        imt.models::latentlogitrestricted,
        data = private$..stan_data,
        seed = seed,
        ... 
      )
      } else {
        sim_out <- rstan::sampling(
        imt.models::latentlogit,
        data = private$..stan_data,
        seed = seed,
        ... 
      )
      }

      private$..prior_eta <- rstan::extract(sim_out, pars = "eta")$eta
      private$..prior_tau <- rstan::extract(sim_out, pars = "tau")$tau
      
      # Fit model
      if (fit) {
        message("Fitting model to the data...")
        stan_data$run_estimation <- 1
        if (restriction) {
          private$..stanfit <- rstan::sampling(
          imt.models::latentlogitrestricted,
          data = stan_data,
          seed = seed,
          ...
        )
        } else {
          private$..stanfit <- rstan::sampling(
          imt.models::latentlogit,
          data = stan_data,
          seed = seed,
          ...
        )
        }

        # Extract Draws
        message("Extracting posterior draws...")
        draws <- rstan::extract(private$..stanfit)
        
        private$..mcmc_checks <- mcmcChecks$new(
          fit = private$..stanfit,
          pars = c("tau", "epsilon_0", "epsilon_1")
        )
        
        private$..tau_draws <- as.vector(draws$tau)
        private$..eta_draws <- as.vector(draws$eta)
        private$..epsilon0_draws <- as.vector(draws$epsilon_0)
        private$..epsilon1_draws <- as.vector(draws$epsilon_1)
        
        # prob_dissatisfied is a [draws x C] matrix
        private$..prob_dissatisfied_draws <- draws$prob_dissatisfied
      }
      return(invisible(self))
    },
    
    #' @description
    #' Plot MCMC trace for key parameters.
    #' @param ... Additional arguments for bayesplot::mcmc_trace
    #' @return A ggplot object.
    tracePlot = function(...) {
      return(
        bayesplot::mcmc_trace(
          private$..stanfit,
          pars = c("tau", "eta", "epsilon_0", "epsilon_1", "alpha"),
          ...
        )
      )
    },
    
    #' @description
    #' Calculates the posterior probability of the ATE (eta) being
    #' greater than, less than, or within a range.
    #' @param a Optional. Lower bound for the threshold.
    #' @param b Optional. Upper bound for the threshold.
    #' @param prior Logical. If TRUE, calculates based on the prior.
    #' @return A character string summarizing the estimated probability
    calcProb = function(a = 0, b = NULL, prior = FALSE) {
      if (is.null(a) && is.null(b)) {
        stop("Either 'a' or 'b' must be provided.")
      }
      if (!is.null(a) && !is.null(b) && b <= a) {
        stop("'b' must be greater than 'a'.")
      }
      
      x <- if (prior) private$..prior_eta else private$..eta_draws
      txt <- if (prior) "Our prior is " else "Given the data, we estimate "
      
      # Convert to percentage points for reporting
      x <- x * 100
      
      if (!is.null(a) && is.null(b)) {
        p <- scales::percent(mean(x > a))
        statement <- glue::glue(
          "{txt}that the probability that the effect is more than {a}",
          " percentage points is {p}."
        )
      } else if (is.null(a) && !is.null(b)) {
        p <- scales::percent(mean(x < b))
        statement <- glue::glue(
          "{txt}that the probability that the",
          " effect is less than {b} percentage points is {p}."
        )
      } else { # both 'a' and 'b' are present
        p <- scales::percent(mean(x > a & x < b))
        statement <- glue::glue(
          "{txt}that the probability that the effect",
          " is between {a} and {b} percentage points is {p}."
        )
      }
      return(statement)
    },
    
    #' @description
    #' Calculate point estimate of the ATE (eta) in percentage points.
    #' @param median Logical. If TRUE (default), return median. Else, mean.
    #' @return A numeric value representing the point estimate.
    pointEstimate = function(median = TRUE) {
      draws <- private$..eta_draws * 100 # Convert to percentage points
      if (median) {
        return(median(draws))
      } else {
        return(mean(draws))
      }
    },
    
    #' @description
    #' Calculates credible interval for the ATE (eta).
    #' @param width Numeric value between 0 and 1 (e.g., 0.95).
    #' @param round Integer for rounding decimal places.
    #' @return A character string with the summary.
    credibleInterval = function(width = 0.75, round = 0) {
      private$..credible_interval <- credibleInterval(
        draws = private$..eta_draws, width
      )
      
      statement <- glue::glue(
        "Given the data, we estimate that there is a ",
        "{scales::percent(width)} probability that the ATE is between ",
        "{round(private$..credible_interval$lower_bound * 100, round)} and ",
        "{round(private$..credible_interval$upper_bound * 100, round)} ",
        "percentage points."
      )
      return(statement)
    },
    
    #' @description
    #' Plots prior and posterior distributions for ATE (eta) or coeff (tau).
    #' @param tau Logical. If TRUE, plot tau instead of eta
    #' @param ... other arguments passed to vizdraws::vizdraws.
    #' @return An interactive plot of the prior and posterior distributions.
    vizdraws = function(tau = FALSE, ...) {
      if (tau) {
        p <- vizdraws::vizdraws(
          prior = private$..prior_tau,
          posterior = private$..tau_draws,
          ...
        )
      } else {
        p <- vizdraws::vizdraws(
          prior = private$..prior_eta * 100,
          posterior = private$..eta_draws * 100,
          ...
        )
      }
      return(p)
    },
    
    #' @description
    #' Plots lollipop chart for the prior and posterior of the ATE (eta)
    #' being greater than a threshold.
    #' @param threshold cutoff (in percentage points)
    #' @param ... other arguments passed to vizdraws::lollipops.
    #' @return A lollipop chart.
    lollipop = function(threshold = 0, ...) {
      data <- data.frame(
        Name = "Impact (ATE)",
        Prior = mean(private$..prior_eta * 100 > threshold),
        Posterior = mean(private$..eta_draws * 100 > threshold)
      )
      p <- vizdraws::lollipops(data, ...)
      return(p)
    },
    
    #' @description
    #' Get posterior summary of AI error rates.
    #' @param width Numeric value for credible interval width (e.g., 0.75).
    #' @return A data.frame summarizing epsilon_0 and epsilon_1.
    getErrorRates = function(width = 0.75) {
      prob_lower <- (1 - width) / 2
      prob_upper <- 1 - prob_lower
      
      summarize_draws <- function(draws, name) {
        tibble::tibble(
          parameter = name,
          mean = mean(draws),
          median = median(draws),
          lower_bound = quantile(draws, probs = prob_lower),
          upper_bound = quantile(draws, probs = prob_upper)
        )
      }
      
      e0_summary <- summarize_draws(private$..epsilon0_draws, "Epsilon_0 (False Positive)")
      e1_summary <- summarize_draws(private$..epsilon1_draws, "Epsilon_1 (False Negative)")
      
      return(dplyr::bind_rows(e0_summary, e1_summary))
    },
    
    #' @description
    #' Get posterior summary of dissatisfaction probabilities for each call.
    #' @param width Numeric value for credible interval width (e.g., 0.75).
    #' @return A data.frame with call_id and posterior summaries.
    getDissatisfactionProbabilities = function(width = 0.75) {
      if (is.null(private$..prob_dissatisfied_draws)) {
        stop("Model has not been fit, or no draws were extracted.")
      }
      
      prob_lower <- (1 - width) / 2
      prob_upper <- 1 - prob_lower
      
      # prob_dissatisfied_draws is [draws x C]
      draws_matrix <- private$..prob_dissatisfied_draws
      
      mean_probs <- colMeans(draws_matrix)
      median_probs <- apply(draws_matrix, 2, median)
      lower_bounds <- apply(draws_matrix, 2, quantile, probs = prob_lower)
      upper_bounds <- apply(draws_matrix, 2, quantile, probs = prob_upper)
      
      return(
        data.frame(
          call_id = private$..call_ids,
          mean_prob_dissatisfied = mean_probs,
          median_prob_dissatisfied = median_probs,
          lower_bound = lower_bounds,
          upper_bound = upper_bounds
        )
      )
    },
    
    #' @description
    #' Plots draws from the prior distribution of tau and eta.
    plotPrior = function() {
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
        ggplot2::xlab("Prior ATE (eta) in percentage points") +
        ggplot2::ylab("N draws") +
        ggplot2::theme_minimal()
      
      plots <- ggpubr::ggarrange(tau, eta, labels = c("tau", "eta"), ncol = 2)
      
      plots <- ggpubr::annotate_figure(
        plots,
        top = ggpubr::text_grob("Draws from prior distributions",
                               face = "bold", size = 14)
      )
      return(plots)
    }
  )
)

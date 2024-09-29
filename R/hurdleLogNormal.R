#' @title Bayesian Hurdle Log-Normal Model Factory
#' @docType class
#' @export
#' @description A class for creating and managing Bayesian Hurdle Log-Normal Models
#' @field version imt package version used to fit model
#' @field ATE_draws Posterior draws for the Average Treatment Effect on the positive outcome
#' @field tau_prob_zero_draws Posterior draws for the change in probability of zero outcome due to treatment
#' @field mcmChecks MCMC diagnostics
#' @field credible_interval Credible interval for the treatment effect
#' @field prior_ATE Prior distribution for ATE
#' @field prior_tau_prob_zero Prior distribution for tau_prob_zero
#' @field predict_list List of predictions

hurdleLogNormal <- R6::R6Class(
  classname = "hurdleLogNormal",
  private = list(
    ..mcmc_checks = NULL,          # A mcmc_checks object
    ..version = NULL,
    ..stanfit = NULL,
    ..ATE_draws = NULL,            # Posterior draws for ATE
    ..tau_prob_zero_draws = NULL,  # Posterior draws for change in zero probability
    ..stan_data = NULL,
    ..credible_interval = NULL,
    ..prior_ATE = NULL,
    ..prior_tau_prob_zero = NULL,
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
    #' @description Get the posterior draws for ATE
    ATE_draws = function() {
      return(private$..ATE_draws)
    },
    #' @description Get the posterior draws for change in zero probability
    tau_prob_zero_draws = function() {
      return(private$..tau_prob_zero_draws)
    },
    #' @description Get the MCMC diagnostics
    mcmChecks = function() {
      return(private$..mcmc_checks)
    },
    #' @description Get the credible interval
    credible_interval = function() {
      return(private$..credible_interval)
    },
    #' @description Get the prior for ATE
    prior_ATE = function() {
      return(private$..prior_ATE)
    },
    #' @description Get the prior for change in zero probability
    prior_tau_prob_zero = function() {
      return(private$..prior_tau_prob_zero)
    },
    #' @description Get the list of predictions
    predict_list = function() {
      return(private$..predict_list)
    }
  ),
  public = list(
    #' @description
    #' Create a new Bayesian Hurdle Log-Normal Model object.
    #'
    #' @param data Data frame to be used
    #' @param y Name of the outcome variable in the data frame
    #' @param x Vector of names of all covariates in the data frame
    #' @param treatment Name of the treatment indicator variable in the data frame
    #' @param mean_alpha_logit Prior mean for alpha in the logit part (default: -3)
    #' @param sd_alpha_logit Prior standard deviation for alpha in the logit part (default: 2)
    #' @param mean_beta_logit Prior mean for beta in the logit part (default: 0 vector)
    #' @param sd_beta_logit Prior standard deviation for beta in the logit part (default: 0.5 vector)
    #' @param tau_mean_logit Prior mean for the treatment effect in the logit part (default: 0)
    #' @param tau_sd_logit Prior standard deviation for the treatment effect in the logit part (default: 0.5)
    #' @param mean_tau Prior mean for the treatment effect in the log-normal part (default: 0)
    #' @param sigma_tau Prior standard deviation for the treatment effect in the log-normal part (default: 0.035)
    #' @param seed Seed for Stan fitting
    #' @param fit Flag for fitting the data to the model or not
    #' @param ... Additional arguments for Stan
    #' @return invisible
    initialize = function(data, y, x, treatment,
                          mean_alpha_logit = -3, sd_alpha_logit = 2,
                          mean_beta_logit = NULL, sd_beta_logit = NULL,
                          tau_mean_logit = 0, tau_sd_logit = 0.5,
                          mean_tau = 0, sigma_tau = 0.035,
                          seed = 1982, fit = TRUE, ...) {

      private$..version <- packageVersion("imt")
      private$..var_cols <- x
      private$..treatment <- treatment
      cleaned_data <- cleanData(
        data = data,
        y = y,
        treatment = treatment,
        x = x,
        binary = FALSE
      )

      # Set default values for mean_beta_logit and sd_beta_logit if NULL
      if (is.null(mean_beta_logit)) {
        mean_beta_logit <- rep(0, cleaned_data$K)
      }
      if (is.null(sd_beta_logit)) {
        sd_beta_logit <- rep(0.5, cleaned_data$K)
      }

      stan_data <- list(
        N = cleaned_data$N,
        y = cleaned_data$Y,
        K = cleaned_data$K,
        X = cleaned_data$X,
        mean_alpha_logit = mean_alpha_logit,
        sd_alpha_logit = sd_alpha_logit,
        mean_beta_logit = mean_beta_logit,
        sd_beta_logit = sd_beta_logit,
        treatment = cleaned_data$treat_vec,
        tau_mean_logit = tau_mean_logit,
        tau_sd_logit = tau_sd_logit,
        mean_tau = mean_tau,
        sigma_tau = sigma_tau,
        run_estimation = 0
      )
      private$..stan_data <- stan_data

      # Draw from the prior
      sim_out <- rstan::sampling(stanmodels$hurdlelognormal,
                                 data = private$..stan_data,
                                 seed = seed,
                                 ...)
      stan_data$run_estimation <- 1

      private$..prior_ATE <- as.matrix(sim_out, pars = "ATE")[, 1]
      private$..prior_tau_prob_zero <- as.matrix(sim_out, pars = "tau_prob_zero")[, 1]

      # Fit model
      if (fit) {
        message("Fitting model to the data")
        private$..stanfit <- rstan::sampling(stanmodels$hurdlelognormal,
                                             data = stan_data,
                                             seed = seed,
                                             ...)
        private$..mcmc_checks <- mcmcChecks$new(
          fit = private$..stanfit,
          pars = c("ATE", "tau_prob_zero")
        )

        private$..ATE_draws <-
          as.data.frame(private$..stanfit, pars = "ATE") |>
          dplyr::pull(ATE)

        private$..tau_prob_zero_draws <-
          as.data.frame(private$..stanfit, pars = "tau_prob_zero") |>
          dplyr::pull(tau_prob_zero)

        # Extract and store predictions
        private$..predict_list <- rstan::extract(private$..stanfit, pars = c("y_pred", "y0", "y1"))
      }
      return(invisible())
    },

    #' @description
    #' Plot MCMC trace for the ATE and tau_prob_zero parameters
    #' @param ... Additional arguments for Stan
    #' @return A ggplot object
    tracePlot = function(...) {
      return(
        bayesplot::mcmc_trace(private$..stanfit,
                              pars = c("ATE", "tau_prob_zero"), ...
        )
      )
    },

    #' @description
    #' Calculates the posterior probability of an effect being greater than, less than,
    #' or within a range defined by thresholds.
    #'
    #' @param effect_type The type of effect to calculate probability for ("ATE" or "tau_prob_zero")
    #' @param a Optional. Lower bound for the threshold.
    #' @param b Optional. Upper bound for the threshold.
    #' @param prior Logical. If TRUE, calculates probabilities based on
    #' the prior distribution. If FALSE (default), uses the posterior distribution.
    #'
    #' @return A character string summarizing the estimated probability
    #'
    calcProb = function(effect_type, a = 0, b = NULL, prior = FALSE) {
      # Input validation
      if (!(effect_type %in% c("ATE", "tau_prob_zero"))) {
        stop("Invalid 'effect_type'. Must be 'ATE' or 'tau_prob_zero'.")
      }
      if (is.null(a) && is.null(b)) {
        stop("Either 'a' or 'b' must be provided.")
      }
      if (!is.null(a) && !is.null(b) && b <= a) {
        stop("'b' must be greater than 'a'.")
      }

      if (effect_type == "ATE") {
        if (prior) {
          x <- private$..prior_ATE
          txt <- "Our prior is "
        } else {
          x <- private$..ATE_draws
          txt <- "Given the data, we estimate "
        }
      } else {  # effect_type == "tau_prob_zero"
        if (prior) {
          x <- private$..prior_tau_prob_zero
          txt <- "Our prior is "
        } else {
          x <- private$..tau_prob_zero_draws
          txt <- "Given the data, we estimate "
        }
      }

      if (!is.null(a) && is.null(b)) {
        p <- scales::percent(mean(x > a))
        statement <- glue::glue(
          "{txt} that the ",
          "probability that the {effect_type} is more than {a}",
          " is {p}."
        )
      } else if (is.null(a) && !is.null(b)) {
        p <- scales::percent(mean(x < b))
        statement <- glue::glue(
          "{txt} that the probability that the",
          " {effect_type} is less than {b} is {p}."
        )
      } else { # both 'a' and 'b' are present
        p <- scales::percent(mean(x > a & x < b))
        statement <- glue::glue(
          "{txt} that the probability that the {effect_type}",
          " is between {a} and {b} is {p}."
        )
      }
      return(statement)
    },

    #' @description
    #' Calculate point estimate of the effect
    #'
    #' @param effect_type The type of effect to calculate the point estimate for ("ATE" or "tau_prob_zero")
    #' @param median Logical value. If TRUE (default), the median of
    #'   the draws is returned. If FALSE, the mean is returned.
    #'
    #' @return A numeric value representing the point estimate.
    #'
    pointEstimate = function(effect_type, median = TRUE) {
      if (!(effect_type %in% c("ATE", "tau_prob_zero"))) {
        stop("Invalid 'effect_type'. Must be 'ATE' or 'tau_prob_zero'.")
      }
      draws <- if (effect_type == "ATE") private$..ATE_draws else private$..tau_prob_zero_draws
      if (median) {
        return(median(draws))
      } else {
        return(mean(draws))
      }
    },

    #' @description
    #' Calculates credible interval for the effect of the intervention
    #'
    #' @param effect_type The type of effect to calculate the credible interval for ("ATE" or "tau_prob_zero")
    #' @param width Numeric value between 0 and 1 representing the desired
    #'   width of the credible interval (e.g., 0.95 for a 95% credible interval).
    #' @param round Integer value indicating the number of decimal places to round
    #'   the lower and upper bounds of the credible interval.
    #'
    #' @return A character string with the following information:
    #'  - The probability associated with the specified width
    #'  - The lower and upper bounds of the credible interval, rounded to the
    #'    specified number of decimal places
    #'
    credibleInterval = function(effect_type, width = 0.95, round = 2) {
      if (!(effect_type %in% c("ATE", "tau_prob_zero"))) {
        stop("Invalid 'effect_type'. Must be 'ATE' or 'tau_prob_zero'.")
      }
      draws <- if (effect_type == "ATE") private$..ATE_draws else private$..tau_prob_zero_draws
      private$..credible_interval <- credibleInterval(
        draws = draws, width
      )
      statement <- glue::glue(
        "Given the data, we estimate that there is a ",
        "{scales::percent(width)} probability that the {effect_type} is between ",
        "{round(private$..credible_interval$lower_bound, round)} and ",
        "{round(private$..credible_interval$upper_bound, round)}."
      )
      return(statement)
    },
    #' @description
    #' Plot prior distributions for ATE and tau_prob_zero
    #' @param bins Number of bins for the histograms
    #' @param xlim_ate Optional. Limits for the x-axis of the ATE histogram
    #' @param xlim_tau Optional. Limits for the x-axis of the tau_prob_zero histogram
    #' @return A list containing two ggplot objects: one for ATE and one for tau_prob_zero
    plotPrior = function(bins = 2000, xlim_ate = NULL, xlim_tau = NULL) {
      # Access prior distributions
      ate_prior_draws <- private$..prior_ATE
      tau_prior_prob_zero_draws <- private$..prior_tau_prob_zero

      priors <- tibble::tibble(
        ate_prior_draws = ate_prior_draws,
        tau_prior_prob_zero_draws = round(tau_prior_prob_zero_draws * 100, 1)
      )

      # Plot prior distributions
      p1 <- ggplot(data = priors, aes(x = ate_prior_draws)) +
        ggplot2::geom_histogram(bins = bins) +
        labs(title = "Prior Distribution of ATE", x = "ATE", y = "Frequency")

      p2 <- ggplot(data = priors, aes(x = tau_prior_prob_zero_draws)) +
        ggplot2::geom_histogram(bins = bins) +
        labs(title = "Prior Distribution of tau_prob_zero (%)",
             x = "tau_prob_zero (%)", y = "Frequency")

      # Apply xlim if provided
      if (!is.null(xlim_ate)) {
        p1 <- p1 + ggplot2::coord_cartesian(xlim = xlim_ate)
      }
      if (!is.null(xlim_tau)) {
        p2 <- p2 + ggplot2::coord_cartesian(xlim = xlim_tau)
      }

      return(list(ate_prior = p1, tau_prior = p2))
    },

    #' @description
    #' Performs Posterior Predictive Checks and generates a density overlay plot
    #'
    #' @param n Number of posterior predictive samples to use in the plot
    #' @param xlim Optional. Limits for the x-axis of the plot
    #' @return A ggplot object representing the density overlay plot
    posteriorPredictiveCheck = function(n = 50, xlim = NULL) {
      # Access the predictions
      y_pred <- private$..predict_list$y_pred
      # Sample 'n' rows from the predictions
      y_pred_sample <- y_pred[sample.int(nrow(y_pred), n), ]
      # Create the plot
      p <- bayesplot::ppc_dens_overlay(
        y = fake_data$watch_time, # Assuming 'fake_data' is accessible within the class
        yrep = as.matrix(y_pred_sample)
      )
      # Apply xlim if provided
      if (!is.null(xlim)) {
        p <- p + ggplot2::coord_cartesian(xlim = xlim)
      }

      return(p)
    }
  )
)

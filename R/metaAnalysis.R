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

#' Create a Meta-Analysis Object Using Data From Previous Studies
#'
#' A meta analysis has raw data and draws from the lift's posterior
#' distribution. This is represented by an R6 Class.
#'
#' @export
metaAnalysis <- R6::R6Class(
  classname = "meta-analysis",
  private = list(
    mean_mu = NULL, # real
    sd_mu = NULL, # real
    data = NULL, # data frame or tibble
    fit = NULL, # stan fit object
    mu_draws = NULL, # vector of real numbers
    ci_width = NULL, # real between 0 and 1
    lift = NULL, # real
    lower_bound = NULL, # real
    upper_bound = NULL, # real
    mcmc_checks = NULL # a mcmc_checks object
  ),
  active = list(
    #' @field PosteriorATE Draws from the posterior distribution of the
    #'   average treatment effect.
    PosteriorATE = function() {
      return(private$mu_draws)
    },
    #' @field checks MCMC diagnostics
    checks = function() {
      return(private$mcmc_checks)
    },
    #' @field CredibleInterval Lower and upper bounds of the credible interval
    CredibleInterval = function() {
      return(c(private$lower_bound, private$upper_bound))
    },
    #' @field PointEstimate Point estimate of the average treatment effect
    PointEstimate = function() {
      return(mean(private$mu_draws))
    },
    #' @field fitted Stan fit object
    fitted = function() {
      return(private$fit)
    }
  ),
  public = list(
    #' @description
    #' Create a new meta analysis object.
    #' @param data Data frame with data point estimates and standard
    #'   errors from studies.
    #' @param point_estimates Name of the variable in the data frame that
    #'   contains the point estimates.
    #' @param standard_errors Name of the variable in the data frame that
    #'   contains the standard errors of the point estimates.
    #' @param id Name of the variable in the data frame that contains the
    #'   id of the studies.
    #' @param mean_mu Prior mean for the true lift in the population.
    #' @param sd_mu Prior mean for the standard deviation of the
    #'   true lift in the population.
    #' @param ci_width Credible interval's width.
    #' @param X Covariates matrix.
    #' @param run_estimation Integer flag to control whether estimation is run (1) or not (0).
    #' @param ... other arguments passed to [rstan::sampling()]
    #' @return A new `meta_analysis` object.
    initialize = function(data, point_estimates, standard_errors, id,
                          mean_mu = 0.0, sd_mu = 0.05, ci_width = 0.75,
                          X = NULL, run_estimation = 1, ...) {
      stopifnot(
        ci_width > 0, ci_width < 1,
        sd_mu > 0
      )

      filtered_data <- tidyr::drop_na(data)
      if (nrow(data) > nrow(filtered_data)) {
        warning(glue::glue(
          "Dropped {nrow(filtered_data) - nrow(data)} observations ",
          "due to missing data"
        ))
      }

      y <- dplyr::pull(filtered_data, {{ point_estimates }})

      s <- dplyr::pull(filtered_data, {{ standard_errors }})

      id <- dplyr::pull(filtered_data, {{ id }})

      private$data <- tibble::tibble(id = id, y = y, s = s)

      private$mean_mu <- mean_mu
      private$sd_mu <- sd_mu
      private$ci_width <- ci_width

      if(is.null(X)) {
        stan_data <- list(
          J = length(y),
          y = y,
          sigma = s,
          mean_mu = mean_mu,
          sd_mu = sd_mu,
          run_estimation = run_estimation)

        private$fit <- rstan::sampling(stanmodels$metaanalysisnox,
                                       data = stan_data,
                                       ...
        )

      }else{
        stop("I don't have that model yet")
      }
      private$mcmc_checks <- mcmcChecks$new(private$fit, pars = "mu")
      private$mu_draws <- rstan::extract(private$fit, "mu")$mu
      private$lift <- mean(private$mu_draws)
      private$lower_bound <- quantile(
        private$mu_draws, (1 - private$ci_width) / 2
      )
      private$upper_bound <- quantile(
        private$mu_draws, 1 - (1 - private$ci_width) / 2
      )
      return(invisible())
    },

    #' @description
    #' Plots the raw data.
    #' @return A plot with point estimates and 95% confidence intervals.
    PlotRawData = function() {
      studies <- dplyr::mutate(private$data,
                               LB = y - 1.96 * s,
                               UB = y + 1.96 * s,
                               tooltip = glue::glue(
                                 "lift is {scales::percent(y)} with a 95% confidence interval ",
                                 "between {scales::percent(LB, accuracy=1)} and ",
                                 "{scales::percent(UB, accuracy=1)}"
                               )
      )

      return(ggplot2::ggplot(data = studies, ggplot2::aes(x = id, y = y)) +
               ggplot2::geom_pointrange(
                 ggplot2::aes(ymin = LB, ymax = UB)
               ) +
               ggplot2::scale_y_continuous(
                 labels = scales::percent_format(accuracy = 1)
               ) +
               ggplot2::ylab("Lift") +
               ggplot2::xlab("Study"))
    },

    #' @description
    #' Plots lift's prior and posterior distributions.
    #'
    #' For more details see [vizdraws::vizdraws()].
    #' @param ... other arguments passed to vizdraws.
    #' @return An interactive plot of the prior and posterior distributions.
    PlotLift = function(...) {
      p <- vizdraws::vizdraws(
        prior = glue::glue("N({private$mean_mu}, {private$sd_mu})"),
        posterior = private$mu_draws,
        ...
      )
      return(p)
    },

    #' @description
    #' Update the width of the credible interval.
    #' @param ci_width New width for the credible interval. This number in the
    #'   (0,1) interval.
    UpdateCI = function(ci_width) {
      stopifnot((ci_width > 0 & ci_width < 1))
      # TODO(martinezig): Write a helper function.
      private$ci_width <- ci_width
      private$lower_bound <- quantile(
        private$mu_draws, (1 - private$ci_width) / 2
      )
      private$upper_bound <- quantile(
        private$mu_draws, 1 - (1 - private$ci_width) / 2
      )
      message(glue::glue(
        "The credible interval has been updated to a ",
        "{scales::percent(ci_width, accuracy = 1)} interval."
      ))
    },

    #' @description
    #' Calculates that probability that lift is between a and b.
    #' @param a Lower bound. By default -Inf.
    #' @param b Upper bound. By default Inf.
    #' @param percent A logical that indicates that a and b should
    #'   be converted to percentage.
    #' @return A string with the probability.
    probability = function(a = -Inf, b = Inf, percent = TRUE) {
      stopifnot(b > a)
      prob <- mean((private$mu_draws > a & private$mu_draws < b))
      a_value <- ifelse(percent, scales::percent(a), a)
      b_value <- ifelse(percent, scales::percent(b), b)
      msg <- dplyr::case_when(
        (is.infinite(a) & is.infinite(b)) ~ glue::glue(
          "Set a or b to something"
        ),
        is.infinite(a) ~ glue::glue(
          "The probability that lift is less than {b_value} is ",
          "{scales::percent(prob, accuracy = 1)}."
        ),
        is.infinite(b) ~ glue::glue(
          "The probability that lift is more than {a_value} is ",
          "{scales::percent(prob, accuracy = 1)}."
        ),
        TRUE ~ glue::glue(
          "The probability that lift is between {a_value} and ",
          "{b_value} is {scales::percent(prob, accuracy = 1)}."
        )
      )
      return(msg)
    },
    #' @description
    #' Calculates the point estimate a credible interval for the meta analysis.
    #' @param percent A logical that indicates that the point estimate should be
    #'  converted to percent.
    #' @return A string with the findings
    findings = function(percent = TRUE) {
      lift <- ifelse(percent, scales::percent(mean(private$mu_draws),
                                              accuracy = 1
      ),
      mean(private$mu_draws)
      )
      lower_bound <- ifelse(percent, scales::percent(private$lower_bound,
                                                     accuracy = 1
      ),
      private$lower_bound
      )
      upper_bound <- ifelse(percent, scales::percent(private$upper_bound,
                                                     accuracy = 1
      ),
      private$upper_bound
      )

      msg <- glue::glue(
        "According to this meta analysis, lift is ",
        "{lift}, with a ",
        "{scales::percent(private$ci_width)} probability ",
        "that is as low as {lower_bound} ",
        "and as high as {upper_bound}."
      )
      return(msg)
    }
  )
)

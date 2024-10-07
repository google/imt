#' Computes the BASIE (BAyeSian Interpretation of Estimates) posterior
#' distribution
#'
#' Implementation of [BASIE](https://www.mathematica.org/download-media?MediaItemId=%7B68A84423-AED0-4349-8524-51E91336507C%7D)
#'
#' @param priorMean Prior distribution mean.
#' @param priorSD Prior distribution standard deviation.
#' @param likMean Likelihood mean (point estimate).
#' @param likSD Likelihood standard deviation (standard error of the point
#'   estimate).
#'
#' @return A list containing the mean and standard deviation of the BASIE
#'   posterior distribution.
#' @export
#'
#' @examples
#' priorMean <- 0
#' priorSD <- 0.03
#' likMean <- 0.071
#' likSD <- 0.074
#'
#' basie_estimate <- basie$new(
#'   priorMean = priorMean,
#'   priorSD = priorSD,
#'   likMean = likMean,
#'   likSD = likSD
#' )
#'
#' basie_estimate$vizdraws(
#'   percentage = TRUE,
#'   breaks = c(0, 0.05),
#'   break_names = c("Worse", "Equivalent", "Better"),
#'   height = 450,
#'   width = 750,
#'   display_mode_name = TRUE
#' )
#'
#' txt <-
#'   glue::glue(
#'     "This example assumes that the prior distribution mean is ",
#'     "{scales::percent(priorMean)}, and its standard deviation ",
#'     "is {scales::percent(priorSD)}. ",
#'     "Furthermore, we imagine a study that found ",
#'     "a point estimate for the effect of {scales::percent(likMean)} ",
#'     "with a standard error of {scales::percent(likSD)}. ",
#'     "Finally, we assume that different decisions would be made ",
#'     "if lift is negative, positive but less than 5%, or greater than 5%."
#'   )
#'
#' cat(stringr::str_wrap(txt), "\n")
#'
#' basie_estimate$plotProbabilities()
#' 

basie <- R6::R6Class(
  classname = "basie",
  private = list(
    prior_mean = NULL,
    prior_variance = NULL,
    likelihood_mean = NULL,
    likelihood_variance = NULL
  ),
  active = list(
    #' @field mean Mean of the BASIE posterior distribution.
    mean = function() {
      (
        1 / private$prior_variance * private$prior_mean +
          1 / private$likelihood_variance * private$likelihood_mean
      ) /
        (1 / private$prior_variance + 1 / private$likelihood_variance)
    },
    #' @field sd Standard deviation of the BASIE posterior distribution.
    sd = function() {
      1 / sqrt(1 / private$prior_variance + 1 / private$likelihood_variance)
    },
    #' @field prior Prior distribution.
    prior = function() {
      glue::glue("N({private$prior_mean}, {sqrt(private$prior_variance)})")
    }
  ),
  public = list(
    #' @description
    #' Estimate the BASIE posterior distribution.
    #'
    #' @param priorMean Prior distribution mean.
    #' @param priorSD Prior distribution standard deviation.
    #' @param likMean Likelihood mean (point estimate).
    #' @param likSD Likelihood standard deviation (standard error of the point
    #'   estimate).
    #' @return An object of class `Basie`.
    initialize = function(priorMean, priorSD, likMean, likSD) {
      private$prior_mean <- priorMean
      private$prior_variance <- priorSD^2
      private$likelihood_mean <- likMean
      private$likelihood_variance <- likSD^2
      invisible(self)
    },
    
    #' @description
    #' Plot the prior and BASIE posterior distributions.
    #' See [vizdraws::vizdraws()] for more details.
    #'
    #' @param draws Number of draws for the posterior.
    #' @param ... Other arguments passed to `vizdraws::vizdraws()`.
    vizdraws = function(draws = 5e4L, ...) {
      vizdraws::vizdraws(
        prior = self$prior,
        posterior = rnorm(
          n = draws,
          mean = self$mean,
          sd = self$sd
        ),
        ...
      )
    },
    
    #' @description
    #' Estimates the probability that the effect of the intervention is
    #'   greater than `x`.
    #'
    #' @param x Threshold of interest.
    probability = function(x) {
      pnorm(q = x, mean = self$mean, sd = self$sd, lower.tail = FALSE)
    },
    
    #' @description
    #' Plot the probability that the effect of the intervention is greater
    #' than a given threshold.
    #'
    #' @param from The starting value for the x-axis.
    #' @param to The ending value for the x-axis.
    #' @param length.out The number of points to generate for the x-axis.
    #' @param ... Other arguments passed to `ggplot2::ggplot()`.
    plotProbabilities = function(from = -0.08, to = 0.08,
                                 length.out = 100, ...) {
      df_plot <- tibble::tibble(x = seq(from = from, to = to,
                                        length.out = length.out)) |>
        dplyr::mutate(pr = self$probability(x))
      
      ggplot2::ggplot(data = df_plot, ggplot2::aes(x = x, y = pr)) +
        ggplot2::geom_line() +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::scale_x_continuous(labels = scales::percent) +
        ggplot2::ylab("Probability") +
        ggplot2::xlab("Impact") +
        ggplot2::ggtitle("Probability that the impact is greater than x") +
        ggplot2::theme_bw(base_size = 16)
    }
  )
)
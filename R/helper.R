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

#' Calculate logit link and sample from binomial distribution
#'
#' This function takes prior / posterior draws of the Bayesian logit model,
#' applies theinverse logit (logistic) transformation to obtain probabilities,
#' and then generates random samples from binomial distributions.
#'
#' @param alpha Numeric. A draw of alpha param from the logit model
#' @param tau Numeric. A draw of tau param from the logit model
#' @param beta Vector. A draw of beta params from the logit model
#' @param treat A 0 / 1 vector of treatment indicator
#' @param X Data to be predicted
#' @param N Numeric. Size of the sample, should depend on size of the data.

#' @return A vector of N samples of preditive draws
#'
logitRng <- function(alpha, tau, beta, treat, X, N) {
  link <- alpha + tau * treat + as.matrix(X) %*% as.matrix(beta)
  return(stats::rbinom(N, size = 1, prob = 1 / (1 + exp(-link))))
}

#' Validate a Logical Subgroup Vector
#'
#' This function checks if a given vector is a logical vector (`TRUE`/`FALSE`)
#' and whether its length matches the number of rows in a specified matrix.
#' It is designed to validate subgroup vectors used for subsettin data.
#'
#' @param subgroup A logical vector representing the subgroup to be validated.
#' @param N Length the subgroup should have.
#' @param name (Optional) A string indicating the name of group.
#'
#' @return The original `subgroup` vector if it passes all validation checks.
#'
#' @details This function performs two key validations:
#'    1. Checks if the `subgroup` vector is logical.
#'    2. Checks if the length of the `subgroup` vector matches the N.
#'
validate_logical_vector <- function(subgroup, N, name = NULL) {
  # Check if 'vec' is a logical vector
  if (!is.logical(subgroup)) {
    stop("Input 'subgroup' must be a logical vector.")
  }

  # Check if length of 'vec' matches input N
  if (length(subgroup) != N) {
    if (is.null(name)) {
      message(
        "Length of 'subgroup' must match the number of rows ", N, "."
      )
    } else {
      message(
        "Length of 'subgroup' must match the number of rows in",
        name, " ."
      )
    }
    stop()
  }

  # If both checks pass, return vec
  return(invisible())
}

#' Calculate Point Estimate (Median or Mean) as Percentage
#'
#' This function computes a point estimate from a numeric vector,
#' returning either the median or the mean as a percentage.
#'
#' @param x A numeric vector containing the data from which to calculate the point estimate.
#' @param median A logical value indicating whether to use the median (default: `TRUE`)
#'               or the mean (`FALSE`) as the point estimate.
#'
#' @return A numeric value representing the chosen point estimate (median or mean)
#'         of the input vector `x`, multiplied by 100 to express it as a percentage.
#'
#' @details This function provides a simple way to obtain either the median or mean
#'          of a numeric vector as a percentage. The choice between these two measures
#'          of central tendency can be controlled by the `median` argument.
#'
pointEstimate <- function(x, median = TRUE) {
  if (median) {
    return(median(x) * 100)
  } else {
    return(mean(x) * 100)
  }
}

#' Calculate Probability of Posterior Draws Falling Within a Range
#'
#' This function estimates the probability that a vector of posterior draws,
#' represented by the parameter `eta`, falls within a specified range.
#' It provides flexibility to use either prior distributions or posterior draws,
#' and to specify one-sided or two-sided probability calculations.
#'
#' @param x A numeric vector containing either posterior draws (default)
#'          or prior samples of the treatment effect parameter (`eta`).
#' @param a (Optional) The lower bound of the range (as a proportion, not percentage).
#' @param b (Optional) The upper bound of the range (as a proportion, not percentage).
#' @param prior A logical value indicating whether to use prior samples (`TRUE`)
#'              or posterior draws (`FALSE`, default) for calculation.
#' @param group_name A string describing the group for which the probability
#'                   is being calculated (default: "group average").
#'
#' @return A formatted string stating the calculated probability and the
#'         specified range. The probability is the proportion of samples
#'         (either prior or posterior) that fall within the defined range.
#'
#' @details This function checks the following cases:
#'  * If both `a` and `b` are `NULL`, it returns an empty string.
#'  * If `b` is less than or equal to `a`, it throws an error.
#'
#' The calculated probability and range are presented in a human-readable string
#' using the `glue` package for formatting.

calcProb <- function(
    x,
    a = NULL,
    b = NULL,
    prior = FALSE,
    group_name = "group average") {
  # Input validation (same as before)
  if (is.null(a) && is.null(b)) {
    statement <- ""
  }
  if (!is.null(a) && !is.null(b) && b <= a) {
    stop("'b' must be greater than 'a'.")
  }
  if (prior) {
    txt <- "Our prior is "
  } else {
    txt <- "Given the data, we estimate "
  }
  if (!is.null(a) && is.null(b)) {
    p <- scales::percent(mean(x > a))
    statement <- glue::glue(
      "{txt} that the ",
      "probability that the {group_name} is more than {a*100}",
      " percentage points is {p}."
    )
  } else if (is.null(a) && !is.null(b)) {
    p <- mean(x < b)
    statement <- glue::glue(
      "{txt} that the probability that the",
      " {group_name} is less than {b*100} percentage points is {p}."
    )
  } else { # both 'a' and 'b' are present
    p <- mean(x > a & x < b)
    statement <- glue::glue(
      "{txt} that the probability that the {group_name}",
      " is between {a*100} and {b*100} percentage points is {p}."
    )
  }
  return(statement)
}

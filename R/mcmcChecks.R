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

#' @title MCMC Checks
#' @description Checks convergence, mixing, effective sample size, and divergent transitions
#' @docType class
#' @export
#' @field everything_looks_fine logical indicating whether all MCMC tests passed.
#' @field diagnostics list of the outcome of each MCMC test
#' @field warnings list of the warning messages from failed MCMC tests
#' @section Methods:
#' \describe{
#'   \item{\code{$new(fit, pars)}}{Runs diagnostics on the supplied \code{stanfit}
#'     object, restricted to parameters identified by the character vector
#'     \code{pars}.
#'     \cr Tests include:
#'     \cr Share of specified parameters with an Rhat less than 1.1. If any
#'     have an Rhat > 1.1, \code{everything_looks_fine} is set to \code{FALSE}.
#'     \cr Share of specified parameters with an n_eff at least 0.1% of the
#'     total number of posterior draws. If any have n_eff < 0.001 * N,
#'     \code{everything_looks_fine} is set to \code{FALSE}.
#'     \cr Share of specified parameters with an n_eff of at least 100. If any
#'     have n_eff < 100, \code{everything_looks_fine} is set to \code{FALSE}.
#'     \cr Number of divergent transitions during posterior sampling. If there
#'     are any whatsoever, \code{everything_looks_fine} is set to
#'     \code{FALSE}.
#'     \cr Share of posterior iterations where the sampler reached the
#'     maximum treedepth. If more than 25\% of iterations maxed out,
#'     \code{everything_looks_fine} is set to \code{FALSE}.}
#' }
mcmcChecks <- R6::R6Class(
  classname = "mcmcChecks",
  private = list(
    ..diagnostics = NULL,
    ..warnings = NULL,
    ..everything_looks_fine = TRUE
  ),
  active = list(
    diagnostics = function() {
      return(private$..diagnostics)
    },
    warnings = function() {
      return(private$..warnings)
    },
    everything_looks_fine = function() {
      return(private$..everything_looks_fine)
    }
  ),
  public = list(
    #' @description
    #' Initialize a new mcmcChecks object and run diagnostics
    #' @param fit A stanfit object to check
    #' @param pars A character vector of parameter names to check
    initialize = function(fit, pars) {
      warnings <- NULL

      # Check that the sampler converged: Rhat < 1.1.
      convergence <- try({
        rhats <- rstan::summary(fit, pars = pars)$summary[, "Rhat"]
        scales::percent(sum(rhats < 1.1) / length(rhats))
      })
      if (convergence != "100%") private$..everything_looks_fine <- FALSE

      # Check that mixing is good enough to enable calculation of the effective
      # sample size: n_eff / N > 0.001
      mixing <- try({
        ratio <- rstan::summary(fit, pars = pars)$summary[, "n_eff"] /
          length(rstan::extract(fit, pars = "lp__")[[1]])
        scales::percent(sum(ratio > 0.001) / length(ratio))
      })
      if (mixing != "100%") private$..everything_looks_fine <- FALSE

      # Check that effective sample sizes are large enough to allow accurate
      # calculation of posterior means, standard deviations, and CIs: n_eff > 100.
      eff <- try({
        n_eff <- rstan::summary(fit, pars = pars)$summary[, "n_eff"]
        scales::percent(sum(n_eff > 100) / length(n_eff))
      })
      if (eff != "100%") private$..everything_looks_fine <- FALSE

      # Number of divergent transitions & treedepth saturations
      n_chain <- fit@sim$chains
      max_td <- fit@stan_args[[1]]$control$max_treedepth
      if (is.null(max_td)) max_td <- 10

      divergent <- 0
      treedepth <- 0
      for (i in 1:n_chain) {
        chain_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)[[i]]
        chain_divergent_count <- sum(chain_params[, "divergent__"])
        divergent <- divergent + chain_divergent_count
        # Identify saturated iterations
        saturated_iterations <- chain_params[, "treedepth__"] == max_td
        # Count saturated iterations and add to total
        chain_saturation_count <- sum(saturated_iterations)
        treedepth <- treedepth + chain_saturation_count
      }
      if (divergent != 0) private$..everything_looks_fine <- FALSE
      # Calculate the total number of iterations excluding warmup
      total_iterations <- n_chain * (fit@sim$iter - fit@sim$warmup)
      # Convert treedepth to percentage and round
      treedepth_percent <- round(100 * treedepth / total_iterations)
      if (private$..everything_looks_fine) {
        message("All diagnostics look fine")
      } else {
        if (divergent != 0) {
          warnings <- c(
            warnings,
            glue::glue(
              "There were {divergent}",
              " divergent transitions after warmup."
            )
          )
        }
        if (convergence != "100%") {
          warnings <- c(
            warnings,
            glue::glue("Only {convergence} of the Rhats are less than 1.1")
          )
        }
        if (mixing != "100%") {
          warnings <- c(
            warnings,
            glue::glue("Only {mixing} of n_eff/N are greater than 0.001")
          )
        }
        if (eff != "100%") {
          warnings <- c(
            warnings,
            glue::glue("Only {eff} of n_eff are greater than 100")
          )
        }
        warning(warnings)
        private$..warnings <- warnings
      }

      private$..diagnostics <- tibble::tibble(
        divergent = divergent,
        convergence = convergence,
        mixing = mixing,
        eff = eff,
        treedepth = scales::percent(treedepth / 100)
      )
      return(invisible())
    }
  )
)

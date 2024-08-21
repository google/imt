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

#' Check Baseline Equivalency.
#'
#' @param data dataframe with the pre-intervention variables,
#'   and the treatment indicator.
#' @param variables vector of with the names of the pre-intervention variables.
#' @param treatment name of the treatment indicator.
#'
#' @return tibble with the standardized difference of the pre-intervention
#'    variables. The tibble includes
#'    variables: the variable name,
#'    std_diff: the standardized difference for that variable as a number
#'    balance: the standardized difference for that variable as a factor
#'             variable.
#' For more details about this methodology check
#' \url{https://ies.ed.gov/ncee/wwc/Docs/OnlineTraining/wwc_training_m3.pdf}.
#' @importFrom stats as.formula lm predict sd
#' @export
#'
checkBaseline <- function(data,
                          variables,
                          treatment) {
  effect_sizes <- calculateEffectSizes(data, treatment_column = treatment,
                                       to_check = variables) |>
    dplyr::mutate(
      balance = dplyr::case_when(
        abs(std_diff) >= 0.25 ~ "Very Concerned",
        abs(std_diff) <= 0.05 ~ "Not Concerned",
        TRUE ~ "Concerned" #
      ),
      balance = factor(
        balance,
        levels = c("Not Concerned", "Concerned", "Very Concerned")
      )
    )
  return(effect_sizes)
}

#' Hedges' g Effect Size with Pooled Standard Deviation
#'
#' Calculates Hedges' g, a standardized effect size for comparing means.
#' This version includes a small-sample correction factor (omega)
#' and uses the pooled standard deviation.
#'
#' @param n_t Numeric value representing the sample size of the treatment group.
#' @param n_c Numeric value representing the sample size of the control group.
#' @param y_t Numeric value representing the mean of the treatment group.
#' @param y_c Numeric value representing the mean of the control group.
#' @param s_t Numeric value representing the standard deviation of the treatment group.
#' @param s_c Numeric value representing the standard deviation of the control group.
#'
#' @return The calculated Hedges' g effect size.
#'
#' @details
#' Hedges' g is a variation of Cohen's d that adjusts for small-sample bias.
#' It is calculated as the difference in means divided by the pooled standard 
#' deviation, then multiplied by a correction factor.
#'
#' @references
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators. 
#' Journal of Educational Statistics, 6(2), 107-128.
#'
#' @export
hedgesG <- function(n_t, n_c, y_t, y_c, s_t, s_c) {

  # Input validation
  if (!is.numeric(c(n_t, n_c, y_t, y_c, s_t, s_c))) {
    stop("All inputs must be numeric values.", call. = FALSE)
  }
  if (any(c(n_t, n_c) <= 0)) {
    stop("Sample sizes must be positive.", call. = FALSE)
  }

  # Hedges' g correction factor
  omega <- 1 - 3 / (4 * (n_t + n_c) - 9)

  # Pooled standard deviation
  s_pooled <- sqrt(((n_t - 1) * s_t^2 + (n_c - 1) * s_c^2) / (n_t + n_c - 2))

  # Hedges' g calculation
  g <- omega * (y_t - y_c) / s_pooled

  return(g)
}

#' Cox's Proportional Hazards Index (Cox's C)
#'
#' Calculates Cox's C, a standardized effect size measure for comparing hazard rates
#' between two groups in survival analysis.
#'
#' @param p_t Numeric value representing the proportion of events (e.g., failures, deaths) in the treatment group.
#' @param p_c Numeric value representing the proportion of events in the control group.
#' @param n_t Numeric value representing the sample size of the treatment group.
#' @param n_c Numeric value representing the sample size of the control group.
#'
#' @return The calculated Cox's C effect size.
#'
#' @details
#' Cox's C is a useful effect size for survival analysis when hazard ratios are not constant over time.
#' It's calculated based on the log odds ratio of events and includes a small sample size correction.
#' The value 1.65 is used to approximate a conversion to a Cohen's d-like scale.
#'
#' @references
#' Cox, D. R. (1972). Regression models and life-tables. Journal of the Royal Statistical Society: Series B (Methodological), 34(2), 187-202.
#'
coxsIndex <- function(p_t, p_c, n_t, n_c) {

  # Input validation (optional but recommended)
  if (!is.numeric(c(p_t, p_c, n_t, n_c))) {
    stop("All inputs must be numeric values.")
  }
  if (p_t <= 0 || p_t >= 1 || p_c <= 0 || p_c >= 1) {
    stop("Proportions of events must be between 0 and 1 (exclusive).")
  }
  if (n_t <= 0 || n_c <= 0) {
    stop("Sample sizes must be positive.")
  }

  # Calculate Cox's C
  omega <- 1 - 3 / (4 * (n_t + n_c) - 9)  # Small sample correction factor

  # Log odds ratio calculation with error handling for extreme proportions (0 or 1)
  log_odds_ratio <- tryCatch(
    log(p_t / (1 - p_t)) - log(p_c / (1 - p_c)),
    error = function(e) {
      if (p_t == 0 || p_t == 1 || p_c == 0 || p_c == 1) {
        warning("Event proportions of 0 or 1 detected. Cox's C may not be appropriate.")
        return(NA) # Or handle differently based on your needs
      } else {
        stop(e)
      }
    }
  )

  # Return Cox's C
  return(omega * log_odds_ratio / 1.65)

}

#' Calculate Effect Sizes for Treatment vs. Control
#'
#' This function calculates effect sizes for each variable in a data frame, 
#' comparing treatment and control groups. It handles continuous, binary, 
#' and categorical variables using appropriate effect size measures.
#'
#' @param data A data frame containing the variables and treatment/control indicator.
#' @param treatment_column The name of the column (as a string) in the data frame that 
#'  indicates whether an observation is in the treatment or control group. This column 
#'  must be a factor with exactly two levels.
#' @param to_check (optional) A character vector specifying the names of the
#' variables for which effect sizes should be calculated. If NULL (default),
#' all variables (except the treatment column) are processed.
#'
#' @return A data frame with two columns:
#'  * `Variable`: The name of each variable in the original data frame.
#'  * `EffectSize`: The calculated effect size for each variable.
#'
#' @details
#'  - For continuous variables, Hedges' g effect size is calculated.
#'  - For binary variables, Cox's Proportional Hazards Index (Cox's C) is calculated.
#'  - For categorical variables, the variable is converted into multiple indicator 
#'   (dummy) variables, and the average Cox's C across these indicators is reported.
#' 
#' @export
calculateEffectSizes <- function(data, treatment_column, to_check = NULL) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }
  if (!treatment_column %in% names(data)) {
    stop("Treatment column not found in data.")
  }
  
  # Check and process treatment column
  treatment_var <- data[[treatment_column]]
  if (is.factor(treatment_var)) {
    if (nlevels(treatment_var) != 2) {
      stop("Treatment column as factor must have exactly two levels.")
    }
  } else if (is.numeric(treatment_var)) {
    unique_values <- unique(treatment_var)
    if (length(unique_values) != 2) {
      stop("Treatment column as numeric must have exactly two unique values.")
    }
    treatment_var <- factor(treatment_var)
  } else if (is.logical(treatment_var)) {
    treatment_var <- factor(treatment_var)
  } else {
    stop("Treatment column must be a factor with two levels, a numeric with two unique values, or a logical.")
  }
  
  # Replace the original treatment column with the processed version
  data[[treatment_column]] <- treatment_var
  
  if (!is.null(to_check) && !all(to_check %in% names(data))) {
    stop("Some variables in 'to_check' are not found in the data.")
  }
  
  # Determine which columns to process (either all or just specified ones)
  cols_to_process <- if (is.null(to_check)) {
    setdiff(names(data), treatment_column)
  } else {
    to_check
  }
  
  # Get unique levels of the treatment column
  levels_treatment <- levels(treatment_var)
  
  # Function to calculate effect size for a single variable
  calculate_effect_size <- function(var, treatment, var_name) {
    if (is.numeric(var)) {
      # Continuous variable: use Hedges' G
      var_by_group <- split(var, treatment)
      effect_size <- hedgesG(
        n_t = length(var_by_group[[2]]), n_c = length(var_by_group[[1]]),
        y_t = mean(var_by_group[[2]]), y_c = mean(var_by_group[[1]]),
        s_t = stats::sd(var_by_group[[2]]), s_c = stats::sd(var_by_group[[1]])
      )
      return(tibble::tibble(variables = var_name, std_diff = effect_size))
    } else if (is.logical(var) || is.factor(var) || is.character(var)) {
      # Binary (logical), factor, or character variable: use Cox's Index
      if (is.character(var)) {
        var <- factor(var)
      }
      if (is.factor(var) && nlevels(var) > 2) {
        # Multi-level factor: create dummy variables
        dummy_data <- stats::model.matrix(~ var - 1)
        # Calculate Cox's Index for each dummy variable
        effect_sizes_dummy <- apply(dummy_data, 2, function(dummy_col) {
          coxsIndex(
            p_t = mean(dummy_col[treatment == levels_treatment[2]]),
            p_c = mean(dummy_col[treatment == levels_treatment[1]]),
            n_t = sum(treatment == levels_treatment[2]),
            n_c = sum(treatment == levels_treatment[1])
          )
        })
        # Return effect sizes for each dummy
        return(tibble::tibble(
          variables = paste0(var_name, "_", colnames(dummy_data)),
          std_diff = effect_sizes_dummy
        ))
      } else {
        # Binary variable (logical or two-level factor): use Cox's Index directly
        effect_size <- coxsIndex(
          p_t = mean(as.numeric(var[treatment == levels_treatment[2]])),
          p_c = mean(as.numeric(var[treatment == levels_treatment[1]])),
          n_t = sum(treatment == levels_treatment[2]),
          n_c = sum(treatment == levels_treatment[1])
        )
        return(tibble::tibble(variables = var_name, std_diff = effect_size))
      }
    } else {
      return(tibble::tibble(variables = var_name, std_diff = NA_real_))
    }
  }
  
  # Calculate effect sizes for all specified variables
  effect_sizes <- purrr::map_dfr(cols_to_process, ~calculate_effect_size(data[[.x]], data[[treatment_column]], .x))
  
  # Return the result tibble
  return(effect_sizes)
}

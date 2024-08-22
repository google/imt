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

#' Combine and Unite Columns
#'
#' This function takes a data frame, identifies columns that are not specified
#' in the exclusion list, and combines them into a new column called 'group'.
#' The original columns used for the combination are then removed.
#' Finally, it returns a data frame with only the 'group',
#' 'variables', 'std_diff', and 'balance' columns.
#'
#' @param df A data frame containing the columns to be combined.
#'
#' @return A data frame with the combined 'group' column and specified columns.
#'
.combineColumns <- function(df) {
  # Identify columns that are NOT 'variables', 'std_diff', or 'balance'
  columns_to_unite <-
    names(df)[!names(df) %in% c("variables", "std_diff", "balance")]

  # Create the 'group' column using unite
  df <- df %>%
    tidyr::unite(
      col = "group",
      dplyr::all_of(columns_to_unite),
      sep = "_",
      remove = TRUE
    )
  # Return the desired columns
  df <- df %>%
    dplyr::select(group, variables, std_diff, balance)

  return(df)
}

#' Add a Random Treatment Indicator Column to a Data Frame
#'
#' This function takes a data frame and adds a new column named "treated" with
#' randomly assigned TRUE/FALSE values. Randomization can be done either on
#' the entire data frame or stratified by specified columns. The probability of
#' being assigned to the treatment group can be specified, with a
#' default of 0.5.
#'
#' @param data The input data frame.
#' @param group_by (Optional) A character vector of column names to stratify the
#'   randomization. If provided, the randomization will be done within
#'   d each groupefined by the specified columns.
#' @param seed (Optional) An integer to set the random seed for reproducibility.
#' @param pr_treated (Optional) The probability of a row being assigned to the
#'   treatment group (TRUE). Default is 0.5.
#'
#' @return A new data frame with the added "treated" column.
#'
#' @importFrom dplyr %>% group_by mutate ungroup
#' @importFrom rlang syms
#'
.randomize_internal <-
  function(data,
           group_by = NULL,
           seed = NULL,
           pr_treated = 0.5) {
    # Set random seed if provided
    if (!is.null(seed)) {
      set.seed(as.integer(seed))
    }
    if (is.null(group_by)) {
      # No grouping, simple randomization
      data <- data %>%
        dplyr::mutate(treated = sample(
          c(TRUE, FALSE),
          nrow(data),
          replace = TRUE,
          prob = c(pr_treated, 1 - pr_treated)
        ))
    } else {
      # Stratified randomization with proper handling of multiple groups
      data <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) %>%
        dplyr::mutate(treated = sample(
          c(TRUE, FALSE),
          dplyr::n(),
          replace = TRUE,
          prob = c(pr_treated, 1 - pr_treated)
        )) %>%
        dplyr::ungroup()
    }
    return(data)
  }

#' Randomly Assign Treatment While Controlling for Baseline Equivalency
#'
#' This function repeatedly randomizes treatment
#' assignment (using `.randomizer`) until
#' baseline equivalency is achieved across specified
#' variables, as measured by the `checkBaseline` function from
#' the `im` package. It can optionally stratify the randomization by
#' specified groups.
#'
#' @param data The input data frame containing pre-intervention variables.
#' @param seed (Optional) An integer to set the random
#' seed for reproducibility of the
#' initial randomization attempt. Subsequent attempts will use new random seeds.
#' @param max_attempts The maximum number of randomization
#' attempts to make before stopping and returning an error.
#' @param variables A vector of the names of the pre-intervention
#' variables to check
#' for baseline equivalency.
#' @param standard The desired level of baseline equivalence.
#' Must be one of "Not Concerned",
#' "Concerned", or "Very Concerned". Default is "Not Concerned".
#' ("Not Concerned", "Concerned", or "Very Concerned").
#'  Must be one of "Not Concerned", "Concerned", or "Very Concerned".
#' @param pr_treated (Optional) The probability of a row being assigned to the
#'  treatment group (TRUE). Default is 0.5.
#' @param group_by (Optional) A character vector of column names to stratify the
#'  randomization. If provided, the randomization will be done within each group
#'  defined by the specified columns.
#'
#' @return A new data frame with the added "treated"
#' column, if baseline equivalency is achieved within the specified number of
#' attempts. Otherwise, an error is thrown.
#'
#' @export
#'
randomize <-
  function(data,
           variables,
           standard = "Not Concerned",
           seed = NULL,
           max_attempts = 100,
           pr_treated = 0.5,
           group_by = NULL) {
    # Valid standard levels with their numeric equivalents
    valid_standards <- c(
      "Not Concerned" = 1,
      Concerned = 2,
      "Very Concerned" = 3
    )
    # Check if input standard is valid
    if (!(standard %in% names(valid_standards))) {
      stop(
        "Invalid 'standard' value. Please choose from:",
        " 'Not Concerned', 'Concerned', or 'Very Concerned'."
      )
    }
    standard_levels <-
      standard_levels <-
      c(
        "Not Concerned" = 1,
        Concerned = 2,
        "Very Concerned" = 3
      )
    numeric_standard <- standard_levels[standard]

    # Initial randomization (no seed set for subsequent
    # attempts for true randomness)
    new_data <-
      .randomize_internal(data,
        seed = seed,
        pr_treated = pr_treated,
        group_by = group_by
      )

    attempts <- 1

    # Use repeat loop for better control and readability
    repeat {
      # Check baseline equivalence based on grouping (or lack thereof)
      if (is.null(group_by)) {
        balance_summary <- imt::checkBaseline(new_data, variables, "treated")
        worst_balance <-
          max(as.numeric(balance_summary$balance))
      } else {
        # Group by all columns specified in group_by
        balance_summaries <- new_data %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) %>%
          dplyr::group_modify(~ imt::checkBaseline(., variables, "treated"))
        worst_balance <-
          max(as.numeric(balance_summaries$balance))
        balance_summaries <- .combineColumns(balance_summaries)
      }

      # Success condition (check if ALL worst balances meet the standard)
      if (all(worst_balance <= numeric_standard)) {
        return(
          list(
            data = new_data,
            seed = seed,
            balance = if (is.null(group_by)) {
              balance_summary
            } else {
              balance_summaries
            }
          )
        )
      }

      # Break if max attempts reached
      if (attempts >= max_attempts) {
        break
      }

      # Re-randomize for next attempt
      new_data <-
        .randomize_internal(data, pr_treated = pr_treated, group_by = group_by)
      attempts <- attempts + 1
    }

    stop(
      "Could not achieve baseline equivalence after",
      max_attempts,
      "attempts."
    )
  }

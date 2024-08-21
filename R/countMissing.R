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

#' Count missing values (NA) in a dataframe
#'
#' This function takes a dataframe as input and returns a tibble
#' summarizing the number of missing values (NA) in each column and
#' the number of rows with at least one missing value.
#'
#' @param df A dataframe to analyze.
#'
#' @return A tibble with the following columns:
#'   - NA_<column_name>: Number of NAs in each original column
#'   - missing_rows: Number of rows with at least one NA across all columns
#'
#' @export
#'
countMissing <- function(df) {
  missing <- df |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.))),
      any_missing = sum(rowSums(is.na(.)) > 0)
    ) |>
    tidyr::pivot_longer(
      cols = tidyselect::everything(), names_to = "Variable",
      values_to = "Count", names_prefix = "_"
    )
  if (missing$any_missing == 0) {
    return("There is no missing data")
  } else {
    missing <- missing |> dplyr::filter(Count > 0)
    return(missing)
  }
}

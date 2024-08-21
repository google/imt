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

#' Cleans and prepares data for analysis
#'
#' This function performs a series of data cleaning and preprocessing steps
#' to ensure the data is suitable for analysis. This includes:
#' * Missing data handling
#' * Variable type checks
#' * Collinearity and zero-variance feature removal
#'
#' @param data A data.frame containing the data to be cleaned.
#' @param y Name of the dependent variable (character).
#' @param treatment Name of the treatment variable (character, should be logical).
#' @param x Names of the covariates to include in the model (character vector, optional).
#' @param binary Should the dependent variable be treated as binary? Default is FALSE
#'
#' @return A list containing the cleaned dataset and relevant metadata:
#'  * \code{N}: The number of observations after cleaning.
#'  * \code{K} The number of covariates after cleaning.
#'  * \code{X} The cleaned covariate matrix.
#'  * \code{treat_vec}: Treatment vector as integers (1 for TRUE, 0 for FALSE).
#'  * \code{Y}: The dependent variable vector.
#' @export

cleanData <- function(data,
                      y,
                      treatment,
                      x = NULL,
                      binary = FALSE) {
  stopifnot(is.data.frame(data))
  # Ensure required arguments are provided
  if (missing(y) || missing(treatment) || missing(data)) {
    stop("Please provide values for 'y', 'treatment', and 'data'")
  }
  # Check number of observations in the original data
  n_orig <- nrow(data)
  # Select columns specified in y, x, and treatment
  selected_columns <- c(y, treatment)
  if (!is.null(x)) selected_columns <- c(selected_columns, x)
  # Check if selected columns exist in the data
  invalid_columns <- setdiff(selected_columns, colnames(data))
  if (length(invalid_columns) > 0) {
    stop(
      "Columns not found in the data: ",
      toString(invalid_columns),
      "."
    )
  }
  # Check if the treatment is of type logical
  if (!is.logical(data[[treatment]])) {
    stop("Variable '", treatment, "' should be logical (TRUE or FALSE).")
  }

  # Check if y is logical when binary is TRUE
  if (binary) {
    if (!is.logical(data[[y]])) {
      stop(
        "If 'binary' is TRUE, the dependent variable '",
        y, "' should be logical (TRUE or FALSE)."
      )
    } # Convert logical to integer (1 for TRUE, 0 for FALSE)
    analysis_df <- data |>
      dplyr::mutate(dplyr::across({{ y }}, ~ as.integer(.)))
  }

  # Convert logical to integer (1 for TRUE, 0 for FALSE)
  analysis_df <- data |>
    dplyr::mutate(dplyr::across({{ treatment }}, ~ as.integer(.))) |>
    dplyr::select(tidyselect::all_of(selected_columns)) |>
    tidyr::drop_na()
  # Warning message
  n_dropped <- n_orig - nrow(analysis_df)
  if (n_dropped > 0) {
    warning(glue::glue("Dropped {n_dropped} records due to missing data"))
  }
  # Combine covariates
  if (!is.null(x)) {
    xs <- paste(x, collapse = "+")
    # Build the model formula
    fml <- paste0(y, " ~ 1 + ", treatment, " + ", xs)
    # Extract model matrix and remove intercept
    covars <- stats::model.matrix(stats::as.formula(fml),
      data = analysis_df
    )[, -1]
    # Identify and remove zero-variance columns
    zero_var_cols <- which(apply(covars, 2, stats::sd) == 0)
    if (length(zero_var_cols) > 0) {
      zero_var_names <- colnames(covars)[zero_var_cols]
      warning("Dropping due to zero variance:", zero_var_names)
      covars <- covars[, -zero_var_cols]
    }
    # Identify and remove collinear columns
    collinear_cols <- caret::findLinearCombos(covars)$remove
    if (length(collinear_cols) > 0) {
      collinear_names <- colnames(covars)[collinear_cols]
      warning("Dropping due to collinearity:", collinear_names)
      covars <- covars[, -collinear_cols]
    }
    # Drop treatment if not already dropped before
    if (!(treatment %in% c(names(zero_var_cols), names(collinear_cols)))) {
      covars <- covars[, -1] # removes treatment from covaraites list
    }
  } else {
    fml <- paste0(y, " ~ 1 + ", treatment)
    covars <- NULL
  }
  treat_vec <- analysis_df |> dplyr::pull(treatment)
  Y <- analysis_df |> dplyr::pull(y)
  if (length(x) == 1) {
    K <- 1
    covars <- matrix(covars, ncol = 1)
  } else if (is.null(x)) {
    K <- 0
  } else {
    K <- ncol(covars)
  }
  cleaned_data <- list(
    N = nrow(analysis_df),
    K = K,
    X = covars,
    treat_vec = treat_vec,
    Y = Y
  )
  return(cleaned_data)
}

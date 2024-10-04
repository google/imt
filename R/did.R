#' Calculate Difference-in-Differences Effect
#'
#' This function calculates the difference-in-differences (DID) estimate of the
#' average treatment effect, along with its standard error. It assumes a
#' simple DID design with two groups (control and treatment) and two time
#' periods (pre and post).
#'
#' @param mean_pre_control Mean outcome for the control group in the pre-treatment period.
#' @param sd_pre_control Standard deviation of the outcome for the control group in the pre-treatment period.
#' @param n_pre_control Sample size for the control group in the pre-treatment period.
#' @param mean_post_control Mean outcome for the control group in the post-treatment period.
#' @param sd_post_control Standard deviation of the outcome for the control group in the post-treatment period.
#' @param n_post_control Sample size for the control group in the post-treatment period.
#' @param mean_pre_treat Mean outcome for the treatment group in the pre-treatment period.
#' @param sd_pre_treat Standard deviation of the outcome for the treatment group in the pre-treatment period.
#' @param n_pre_treat Sample size for the treatment group in the pre-treatment period.
#' @param mean_post_treat Mean outcome for the treatment group in the post-treatment period.
#' @param sd_post_treat Standard deviation of the outcome for the treatment group in the post-treatment period.
#' @param n_post_treat Sample size for the treatment group in the post-treatment period.
#'
#' @return A list containing the following components:
#'   \item{did_estimate}{The DID estimate of the average treatment effect.}
#'   \item{se_did}{The standard error of the DID estimate.}
#'
#' @examples
#' \donttest{
#' CalculateDIDEffect(
#'   mean_pre_control = 10, sd_pre_control = 2, n_pre_control = 50,
#'   mean_post_control = 12, sd_post_control = 2.5, n_post_control = 50,
#'   mean_pre_treat = 11, sd_pre_treat = 2.1, n_pre_treat = 50,
#'   mean_post_treat = 15, sd_post_treat = 2.6, n_post_treat = 50
#' )
#' }
#'
#' @export
CalculateDIDEffect <- function(mean_pre_control, sd_pre_control, n_pre_control,
                               mean_post_control, sd_post_control, n_post_control,
                               mean_pre_treat, sd_pre_treat, n_pre_treat,
                               mean_post_treat, sd_post_treat, n_post_treat) {
  # Calculate the change in means for the control group
  delta_control <- mean_post_control - mean_pre_control
  
  # Calculate the change in means for the treatment group
  delta_treat <- mean_post_treat - mean_pre_treat
  
  # Calculate the difference-in-differences estimate (average treatment effect)
  did_estimate <- delta_treat - delta_control
  
  # Calculate the pooled variance for each period
  pooled_var_pre <- ((n_pre_control - 1) * sd_pre_control^2 + (n_pre_treat - 1) * sd_pre_treat^2) /
    (n_pre_control + n_pre_treat - 2)
  pooled_var_post <- ((n_post_control - 1) * sd_post_control^2 + (n_post_treat - 1) * sd_post_treat^2) /
    (n_post_control + n_post_treat - 2)
  
  # Calculate the standard error of the DID estimate
  se_did <- sqrt(pooled_var_pre / n_pre_control + pooled_var_pre / n_pre_treat +
                   pooled_var_post / n_post_control + pooled_var_post / n_post_treat)
  
  # Return a list with the estimate and standard error
  return(list(
    did_estimate = did_estimate,
    se_did = se_did
  ))
}

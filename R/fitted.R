#' Extract Fitted Values of a Dynamite Model
#'
#' Note that these are conditional on the observed data in `newdata`,
#' i.e., these are one-step estimates E(y_t|y_t-1,...,y_1). Often
#' [predict.dynamitefit] is what you want.
#'
#' @export
#' @inheritParams predict.dynamitefit
#' @param newdata \[`data.frame`]\cr Data used in predictions.
#'   If `NULL` (default), the data used in model estimation is used for
#'   predictions as well,
#'   There should be no new time points that were not present in the data that
#'   were used to fit the model, and no new group levels can be included.
#' @examples
#' fitted(gaussian_example_fit, n_draws = 2L)
#' \dontrun{
#' set.seed(1)
#' fit <- dynamite(obs(LakeHuron ~ 1, "gaussian") + lags(),
#'   data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
#'   "id", "time", chains = 1, refresh = 0)
#'
#' # One-step ahead samples (fitted values) from the posterior
#' # (first time point is fixed due to lag in the model):
#' library(ggplot2)
#' fitted(fit) |>
#'   dplyr::filter(time > 2) |>
#'   ggplot(aes(time, LakeHuron_fitted, group = .draw)) +
#'   geom_line(alpha = 0.5) +
#'   # observed values
#'   geom_line(aes(y = LakeHuron), colour = "tomato") +
#'   theme_bw()
#'
#' # Posterior predictive distribution given the first time point:
#' predict(fit, type = "mean") |>
#'   dplyr::filter(time > 2) |>
#'   ggplot(aes(time, LakeHuron_mean, group = .draw)) +
#'   geom_line(alpha = 0.5) +
#'   # observed values
#'   geom_line(aes(y = LakeHuron), colour = "tomato") +
#'   theme_bw()
#' }
#' @srrstats {RE4.9} Returns the fitted values.
fitted.dynamitefit <- function(object, newdata = NULL, n_draws = NULL, ...) {
  predict_dynamitefit(
    object,
    newdata,
    type = "mean",
    eval_type = "fitted",
    impute = "none",
    new_levels = "none",
    n_draws
  )
}

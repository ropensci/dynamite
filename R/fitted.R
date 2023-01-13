#' Extract Fitted Values of a Dynamite Model
#'
#' Fitted values for a `dynamitefit` object. Note that these are conditional on
#' the observed data in `newdata`,i.e., these are one-step estimates
#' E(y_t|y_t-1,...,y_1). Often [dynamite::predict.dynamitefit()] is what you
#' want.
#'
#' @export
#' @inheritParams predict.dynamitefit
#' @param newdata \[`data.frame`]\cr Data used in predictions.
#'   If `NULL` (default), the data used in model estimation is used for
#'   predictions as well,
#'   There should be no new time points that were not present in the data that
#'   were used to fit the model, and no new group levels can be included.
#' @return A `data.frame` containing the fitted values.
#' @srrstats {RE4.9} Returns the fitted values.
#' @examples
#' fitted(gaussian_example_fit, n_draws = 2L)
#' \donttest{
#' set.seed(1)
#' fit <- dynamite(
#'   dformula = obs(LakeHuron ~ 1, "gaussian") + lags(),
#'   data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
#'   group = "id",
#'   time = "time",
#'   chains = 1,
#'   refresh = 0
#' )
#'
#' if (requireNamespace("dplyr") &&
#'     requireNamespace("tidyr") &&
#'     base::getRversion() >= "4.1.0") {
#'   # One-step ahead samples (fitted values) from the posterior
#'   # (first time point is fixed due to lag in the model):
#'   fitted(fit) |>
#'     dplyr::filter(time > 2) |>
#'     ggplot2::ggplot(ggplot2::aes(time, LakeHuron_fitted, group = .draw)) +
#'     ggplot2::geom_line(alpha = 0.5) +
#'     # observed values
#'     ggplot2::geom_line(ggplot2::aes(y = LakeHuron), colour = "tomato") +
#'     ggplot2::theme_bw()
#'
#'   # Posterior predictive distribution given the first time point:
#'   predict(fit, type = "mean") |>
#'     dplyr::filter(time > 2) |>
#'     ggplot2::ggplot(ggplot2::aes(time, LakeHuron_mean, group = .draw)) +
#'     ggplot2::geom_line(alpha = 0.5) +
#'     # observed values
#'     ggplot2::geom_line(ggplot2::aes(y = LakeHuron), colour = "tomato") +
#'     ggplot2::theme_bw()
#' }
#' }
#'
fitted.dynamitefit <- function(object, newdata = NULL, n_draws = NULL,
                               expand = TRUE, df = TRUE, ...) {
  stopifnot_(
    !is.null(object$stanfit),
    "No Stan model fit is available."
  )
  stopifnot_(
    checkmate::test_flag(x = expand),
    "Argument {.arg expand} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_flag(x = df),
    "Argument {.arg df} must be a single {.cls logical} value."
  )
  initialize_predict(
    object,
    newdata,
    type = "mean",
    eval_type = "fitted",
    funs = list(),
    impute = "none",
    new_levels = "none",
    global_fixed = FALSE,
    n_draws,
    expand,
    df
  )
}

#' Extract Fitted Values of a \pkg{dynamite} Model
#'
#' Fitted values for a `dynamitefit` object, i.e.,
#' \eqn{E(y_t | newdata, \theta)} where \eqn{\theta} contains all the
#' model parameters. See also [predict.dynamitefit()] for multi-step
#' predictions.
#'
#' @export
#' @family prediction
#' @inheritParams predict.dynamitefit
#' @param newdata \[`data.frame`]\cr Data used in predictions.
#'   If `NULL` (default), the data used in model estimation is used for
#'   predictions as well.
#'   There should be no new time points that were not present in the data that
#'   were used to fit the model, and no new group levels can be included.
#' @return A `data.frame` containing the fitted values.
#' @srrstats {RE4.9} Returns the fitted values.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' fitted(gaussian_example_fit, n_draws = 2L)
#' \donttest{
#' set.seed(1)
#' # Please update your rstan and StanHeaders installation before running
#' # on Windows
#' if (!identical(.Platform$OS.type, "windows")) {
#'   fit <- dynamite(
#'     dformula = obs(LakeHuron ~ 1, "gaussian") + lags(),
#'     data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
#'     time = "time",
#'     group = "id",
#'     chains = 1,
#'     refresh = 0
#'   )
#'
#'   if (requireNamespace("dplyr") && requireNamespace("tidyr")) {
#'
#'     # One-step ahead samples (fitted values) from the posterior
#'     # (first time point is fixed due to lag in the model):
#'     f <- dplyr::filter(fitted(fit), time > 2)
#'     ggplot2::ggplot(f, ggplot2::aes(time, LakeHuron_fitted, group = .draw)) +
#'       ggplot2::geom_line(alpha = 0.5) +
#'       # observed values
#'       ggplot2::geom_line(ggplot2::aes(y = LakeHuron), colour = "tomato") +
#'       ggplot2::theme_bw()
#'
#'     # Posterior predictive distribution given the first time point:
#'     p <- dplyr::filter(predict(fit, type = "mean"), time > 2)
#'     ggplot2::ggplot(p, ggplot2::aes(time, LakeHuron_mean, group = .draw)) +
#'       ggplot2::geom_line(alpha = 0.5) +
#'       # observed values
#'       ggplot2::geom_line(ggplot2::aes(y = LakeHuron), colour = "tomato") +
#'       ggplot2::theme_bw()
#'   }
#' }
#' }
#'
fitted.dynamitefit <- function(object, newdata = NULL, n_draws = NULL, thin = 1,
                               expand = TRUE, df = TRUE, drop = TRUE, ...) {
  stopifnot_(
    !missing(object),
    "Argument {.arg object} is missing."
  )
  stopifnot_(
    is.dynamitefit(object),
    "Argument {.arg object} must be a {.cls dynamitefit} object."
  )
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
  stopifnot_(
    checkmate::test_flag(x = drop),
    "Argument {.arg drop} must be a single {.cls logical} value."
  )
  if (thin > 1L) {
    idx_draws <- seq.int(1L, ndraws(object), by = thin)
  } else {
    n_draws <- check_ndraws(n_draws, ndraws(object))
    idx_draws <- ifelse_(
      identical(n_draws, ndraws(object)),
      seq_len(n_draws),
      object$permutation[seq_len(n_draws)]
    )
  }
  initialize_predict(
    object,
    newdata,
    type = "mean",
    eval_type = "fitted",
    funs = list(),
    impute = "none",
    new_levels = "none",
    global_fixed = FALSE,
    idx_draws,
    expand,
    df,
    drop
  )
}

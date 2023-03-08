#' Summary for a Dynamite Model Fit
#'
#' The `summary` method provides statistics of the posterior samples of the
#' model; this is an alias of [dynamite::as.data.frame.dynamitefit()] with
#' `summary = TRUE`.
#'
#' @export
#' @family output
#' @rdname dynamite
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @return `summary` returns a `data.frame`.
#' @srrstats {BS6.4, RE4.18} Implements `summary` method.
#' @examples
#' summary(gaussian_example_fit,
#'   types = "beta",
#'   probs = c(0.05, 0.1, 0.9, 0.95)
#' )
#'
summary.dynamitefit <- function(object, ...) {
  as.data.frame.dynamitefit(object, summary = TRUE, ...)
}

#' Summary for a Dynamite Model Fit
#'
#' The `summary` method provides statistics of the posterior samples of the
#' model; this is an alias of [dynamite::as.data.frame.dynamitefit()].
#'
#' @export
#' @rdname dynamite
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param ... Further arguments to [dynamite::as.data.frame.dynamitefit()].
#' @return `summary` returns a `data.frame`.
#' @srrstats {BS6.4, RE4.18} Implements `summary` method.
#' @examples
#' summary(gaussian_example_fit, types = "beta")
#'
summary.dynamitefit <- function(object, ...) {
  as.data.frame.dynamitefit(object, ...)
}

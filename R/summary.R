#' Summary for a Dynamite Model Fit
#'
#' Returns the summary statistics of the posterior samples of the model; this
#' is an alias of `as.data.frame(object)`.
#'
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param ... Further arguments to the `as.data.frame.dynamitefit`.
#' @export
#' @examples
#' summary(gaussian_example_fit, types = "beta")
#'
#' @srrstats {BS6.4, RE4.18} Implements `summary` method.
summary.dynamitefit <- function(object, ...) {
  stopifnot_(
    is.dynamitefit(object),
    "Argument {.arg object} must be a {.cls dynamitefit} object."
  )
  if (!is.null(object$stanfit)) {
    as.data.frame(object, ...)
  } else {
    message_("No Stan model fit is available.")
  }
}

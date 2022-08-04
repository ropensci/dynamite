#' Extract the Number of Observations Used to Fit a Dynamite Model
#'
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param ... Not used.
#' @return Total number of observations as an integer.
#'   Missing values are not accounted for as the
#'   number of complete cases may vary across
#'   channels, time points, and groups.
#' @export
#' @examples
#' nobs(gaussian_example_fit)
#'
#' @srrstats {RE4.5} Provides number a crude number of observations.
nobs.dynamitefit <- function(object, ...) {
  stopifnot_(
    is.dynamitefit(object),
    "Argument {.var object} must be a {.cls dynamitefit} object."
  )
  nrow(object$data)
}

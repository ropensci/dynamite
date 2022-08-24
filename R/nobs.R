#' Extract the Number of Observations Used to Fit a Dynamite Model
#'
#' @export
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param ... Not used.
#' @return Total number of observations as an integer.
#'   Missing values are not accounted for as the
#'   number of complete cases may vary across
#'   channels, time points, and groups.
#' @srrstats {RE4.5} Provides number a crude number of observations.
#' @examples
#' nobs(gaussian_example_fit)
#'
nobs.dynamitefit <- function(object, ...) {
  stopifnot_(
    !missing(object),
    "Argument {.arg object} is missing."
  )
  stopifnot_(
    is.dynamitefit(object),
    "Argument {.var object} must be a {.cls dynamitefit} object."
  )
  nrow(object$data)
}

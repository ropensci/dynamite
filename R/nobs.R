#' Extract the Number of Observations Used to Fit a Dynamite Model
#'
#' @export
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param ... Not used.
#' @return Total number of non-missing observations as an `integer`.
#' @srrstats {RE4.5} Provides the number of observations.
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
  n_obs_vars <- grep(
    pattern = "n_obs_.+",
    x = names(object$stan$sampling_vars),
    value = TRUE
  )
  sum(
    vapply(
      n_obs_vars,
      function(x) {
        sum(object$stan$sampling_vars[[x]])
      },
      integer(1L)
    )
  )
}

#' Return the Number of Posterior Draws of a `dynamitefit` Object
#'
#' @export
#' @export ndraws
#' @family output
#' @aliases ndraws
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @return Number of posterior draws as a single `integer` value.
#' @examples
#' ndraws(gaussian_example_fit)
#'
ndraws.dynamitefit <- function(x) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  stopifnot_(
    !is.null(x$stanfit),
    "No Stan model fit is available."
  )
  as.integer(
    (x$stanfit@sim$n_save[1L] - x$stanfit@sim$warmup2[1L]) *
      x$stanfit@sim$chains
  )
}

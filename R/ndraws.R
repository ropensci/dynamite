#' Return the Number of Posterior Draws of a dynamitefit Object
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @export
#' @export ndraws
#' @rdname ndraws-dynamitefit
#' @aliases ndraws
#' @method ndraws dynamitefit
#' @examples
#' ndraws(gaussian_example_fit)
ndraws.dynamitefit <- function(x) {
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  if (!is.null(x$stanfit)) {
    as.integer(
      (x$stanfit@sim$n_save[1] - x$stanfit@sim$warmup2[1]) *
        x$stanfit@sim$chains
    )
  } else {
    message_("No Stan model fit is available.")
  }
}

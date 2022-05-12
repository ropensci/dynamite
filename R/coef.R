#' Extract Regression Coefficients of Dynamite Model
#'
#' @export
#' @param object An object of class \code{dynamitefit}.
#' @param type  \[`character(1)`]\cr Either `beta` (the default) for
#'   time-invariant coefficients or `delta` for time-varying coefficients.
#' @param summary \[`logical(1)`]\cr If `TRUE` (default), returns posterior
#'   mean and lower and upper limits of the 95% posterior intervals for all
#'   parameters. If `FALSE`, returns all the posterior samples instead.
#' @param ... Ignored.
#' @importFrom stats coef
#' @examples
#' data(gaussian_example_fit)
#' betas <- coef(gaussian_example_fit, type = "beta")
#' deltas <- coef(gaussian_example_fit, type = "delta")
#'
coef.dynamitefit <- function(object, type = c("beta", "delta"),
                             summary = TRUE, ...) {
  type <- match.arg(type)
  as.data.frame(object, types = type, summary = summary)
}

#' Extract regression coefficients of dynamitefit
#'
#'
#' @export
#' @param object An object of class \code{dynamitefit}.
#' @param ... Ignored.
#' @importFrom stats coef
coef.dynamitefit <- function(object, ...) {
  # TODO is this best format?
  as.data.frame(object, parameter_types = "beta")
}

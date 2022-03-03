#' Extract regression coefficients of btvcmfit
#'
#'
#' @export
#' @param object An object of class \code{btvcmfit}.
#' @param ... Ignored.
#' @importFrom stats coef
coef.btvcmfit <- function(object, ...) {
    # TODO is this best format?
    as.data.frame(object, parameter_types = "beta")
}

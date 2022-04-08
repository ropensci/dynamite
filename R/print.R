#' Print Summary of btvcm Object
#'
#' Prints the summary information of time-invariant model parameters.
#'
#' @param x An output from \code{\link{btvcm}}.
#' @param pars A character vector of parameter names.
#' The default is to print summary of spline coefficient standard deviations tau.
#' @param ... Additional arguments to \code{\link{print.stanfit}}.
#' @method print btvcmfit
#' @export
print.btvcmfit <- function(x, pars, ...) {
    if (!is.null(x$stanfit)) {
        if (missing(pars)) {
            all_pars <- x$stanfit@sim$pars_oi
            idx <-  c(grep("^tau", all_pars), grep("lambda",all_pars))
            pars <- all_pars[idx]
        }
        print(x$stanfit, pars = pars, ...)
    } else {
        message("No Stan model fit is available.")
        invisible(x)
    }
}

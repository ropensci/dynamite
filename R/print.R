#' Print Summary of dynamite Object
#'
#' Prints the summary information of time-invariant model parameters.
#'
#' @param x An output from \code{\link{dynamite}}.
#' @param pars A character vector of parameter names.
#' The default is to print summary of spline coefficient standard deviations tau.
#' @param ... Additional arguments to \code{\link{print.stanfit}}.
#' @method print dynamitefit
#' @export
print.dynamitefit <- function(x, pars, ...) {
    if (!is.null(x$stanfit)) {
        if (missing(pars)) {
            all_pars <- x$stanfit@sim$pars_oi
            idx <-  c(grep("^tau", all_pars), grep("lambda",all_pars))
            pars <- all_pars[idx]
        }
        if (length(idx) == 0) {
          message("Model does not contain time-invariant parameters. ")  # TODO what about time-invariant betas...
        } else {
          print(x$stanfit, pars = pars, ...)
        }
    } else {
        message("No Stan model fit is available.")
        invisible(x)
    }
}

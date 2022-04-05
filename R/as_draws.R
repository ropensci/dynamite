#' Convert \code{btvcm} Output to \code{draws_df} Format
#'
#' Converts the output from \code{btvcm} call to a
#' \code{draws_df} format of the \code{posterior} package, enabling the use
#' of diagnostics and plotting methods of \code{posterior} and \code{bayesplot}
#' packages.
#'
#' @param x An object of class \code{btvcmfit}.
#' @param times A vector of indices defining which time points to return?
#' Default is all.
#' @param parameter_types What parameters should be returned? Possible choices are
#' \code{"beta"}, \code{"tau"}, ... #TODO
#' @param ... Ignored.
#' @return A \code{draws_df} object.
#' @importFrom posterior as_draws as_draws_df
#' @aliases as_draws as_draws_df
#' @export
#' @export as_draws_df
#' @rdname as_draws-btvcmfit
#' @method as_draws_df btvcmfit
as_draws_df.btvcmfit <- function(x, parameter_types, ...) {
    if (!requireNamespace("posterior", quietly = TRUE)) {
        stop("This function depends on the 'posterior' package. ", call. = FALSE)
    }
    posterior::as_draws(as.data.frame(x, parameter_types = parameter_types))
    # Not tested
}
#' @export
#' @export as_draws
#' @rdname as_draws-btvcmfit
#' @method as_draws btvcmfit
as_draws.btvcmfit <- function(x, ...) {
    as_draws_df.btvcmfit(x, ...)
}

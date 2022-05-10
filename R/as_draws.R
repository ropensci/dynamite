#' Convert \code{dynamite} Output to \code{draws_df} Format
#'
#' Converts the output from \code{dynamite} call to a
#' \code{draws_df} format of the \code{posterior} package, enabling the use
#' of diagnostics and plotting methods of \code{posterior} and \code{bayesplot}
#' packages. Note that this function extracts all variables from the model fit,
#' and the naming of variables coincides with the internal declarations. For
#' more user-friendly object, use \code{as.data.frame} method on the output of
#' the \code{dynamite} function.
#'
#' @param x An object of class \code{dynamitefit}.
#' @param ... Ignored.
#' @return A \code{draws_df} object.
#' @importFrom posterior as_draws as_draws_df
#' @aliases as_draws as_draws_df
#' @export
#' @export as_draws_df
#' @rdname as_draws-dynamitefit
#' @method as_draws_df dynamitefit
as_draws_df.dynamitefit <- function(x, ...) {

  posterior::as_draws_df(rstan::extract(x$stanfit, permuted = FALSE))
}
#' @export
#' @export as_draws
#' @rdname as_draws-dynamitefit
#' @method as_draws dynamitefit
as_draws.dynamitefit <- function(x, ...) {
  as_draws_df.dynamitefit(x, ...)
}

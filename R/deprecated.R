#' Deprecated Functions in the dynamite Package
#'
#' These functions are provided for compatibility with older versions of the
#' package. They will eventually be completely removed.
#'
#' @rdname dynamite-deprecated
#' @name dynamite-deprecated
#' @usage plot_betas(x, ...)
#' plot_deltas(x,  ...)
#' plot_nus(x, ...)
#' plot_lambdas(x,  ...)
#' plot_psis(x, ...)
#' @return A `ggplot` object.
#' @seealso See [dynamite::plot.dynamitefit()] for documentation of the
#'   parameters of these functions
#' @export  plot_betas plot_deltas plot_nus plot_lambdas plot_psis
#' @aliases plot_betas plot_deltas plot_nus plot_lambdas plot_psis
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param ... Not used.
#' @section Details:
#'
#'   * `plot_betas` is now called via `plot(., plot_type = "beta")`
#'   * `plot_deltas` is now called via `plot(., plot_type = "delta")`
#'   * `plot_nus` is now called via `plot(., plot_type = "nu")`
#'   * `plot_lambdas` is now called via `plot(., plot_type = "lambda")`
#'   * `plot_psis` is now called via `plot(., plot_type = "psi")`
#'
plot_betas <- function(x, ...) {
  .Deprecated("plot")
  plot(x, plot_type = "beta", ...)
}
plot_deltas <- function(x, ...) {
  .Deprecated("plot")
  plot(x, plot_type = "delta", ...)
}
plot_nus <- function(x, ...) {
  .Deprecated("plot")
  plot(x, plot_type = "nu", ...)
}
plot_lambdas <- function(x, ...) {
  .Deprecated("plot")
  plot(x, plot_type = "lambda", ...)
}
plot_psis <- function(x, ...) {
  .Deprecated("plot")
  plot(x, plot_type = "psi", ...)
}

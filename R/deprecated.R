#' Deprecated Functions in the dynamite Package
#'
#' These functions are provided for compatibility with older versions of the
#' package. They will eventually be completely removed.
#'
#' @rdname dynamite-deprecated
#' @name dynamite-deprecated
#' @docType package
#' @usage plot_betas(x, ...)
#' plot_deltas(x, ...)
#' plot_nus(x, ...)
#' plot_lambdas(x, ...)
#' plot_psis(x, ...)
#' @export plot_betas plot_deltas plot_nus plot_lambdas plot_psis
#' @aliases plot_betas plot_deltas plot_nus plot_lambdas plot_psis
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param parameters \[`charecter()`]\cr Parameter name(s) for which the plots
#'   should be drawn. Possible options can be found with the function
#'   [dynamite::get_parameter_names()]. The default is all parameters of a
#'   specific type, limited by `n_params`.
#' @param responses \[`character()`]\cr Response(s) for which the plots should
#'   be drawn. Possible options are `unique(x$priors$response)`. Default is
#'   all responses. Ignored if the argument `parameters` is supplied.
#' @param type \[`character(1)`]\cr Type of the parameter for which the plots
#'   should be drawn. Possible options can be found with the function
#'   [dynamite::get_parameter_types()]. Ignored if the argument `parameters`
#'   is supplied or if `plot_type` is not `"default"`.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity level for `geom_ribbon`.
#'   Default is 0.5.
#' @param scales \[`character(1)`]\cr Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), plots also
#'   the time-varying alphas if such parameters exists in the model.
#' @param include_cutpoints \[`logical(1)`]\cr If `TRUE` (default), plots also
#'   the cutpoints if such parameters exists in the model (either time-varying
#'   or time-invariant).
#' @param times \[`integer()`]\cr Time point(s) for which the plots should be
#'   drawn for time-varying parameters. By default, all time points are
#'   included, up to the maximum number of parameters specified by `n_params`
#'   starting from the first non-fixed time point.
#' @param groups \[`character(1)`]\cr Group name(s) for which the plots
#'   should be drawn when `plot_type` is `"nu"`. The default is all groups.
#' @param n_params \[`integer()`]\cr Maximum number of parameters to plot.
#'   The default value is set by `plot_type`: 5 for `"default"`, 50 for
#'   `"beta"`, `"nu"` and `"lambda"`, and 3 for `"delta"` and `"psi"`. This
#'   argument is intended to prevent accidental plots that may be very large
#'   and time consuming to plot. Please use the `parameters`, `times`, and
#'   `groups` arguments to fine-tune which parameters to plot.
#' @param ... Not used..
#' @return A `ggplot` object.
#' @seealso [dynamite::plot.dynamitefit()]
#' @section Details:
#'
#'   * `plot_betas` is now called via `plot(., plot_type = "beta")`
#'   * `plot_deltas` is now called via `plot(., plot_type = "delta")`
#'   * `plot_nus` is now called via `plot(., plot_type = "nu")`
#'   * `plot_lambdas` is now called via `plot(., plot_type = "lambda")`
#'   * `plot_psis` is now called via `plot(., plot_type = "psi")`
#'
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' plot_betas(gaussian_example_fit)
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

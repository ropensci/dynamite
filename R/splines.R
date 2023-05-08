#' Define the B-splines Used for the Time-varying Coefficients of the Model.
#'
#' This function can be used as part of [dynamiteformula()] to define the
#' splines used for the time-varying coefficients \eqn{\delta}.
#'
#' @export
#' @family formulas
#' @param df \[`integer(1)`]\cr Degrees of freedom, i.e., the total number of
#'   spline coefficients. See [splines::bs()]. Note that the knots are always
#'   defined as an equidistant sequence on the interval starting from the first
#'   non-fixed time point to the last time point in the data. See
#'   [dynamiteformula()] for more information on fixed time points.
#'   Should be an (unrestricted) positive integer.
#' @param degree \[`integer(1)`]\cr See [splines::bs()]. Should be an
#'   (unrestricted) positive integer.
#' @param lb_tau \[`numeric()`]\cr Hard constraint(s) on the lower bound of the
#'   standard deviation parameters \eqn{\tau} of the random walk priors. Can be
#'   useful in avoiding divergences in some cases. See also the `noncentered`
#'   argument. Can be a single positive value, or vector defining the
#'   lower bound separately for each channel, even for channels without
#'   varying effects. The ordering is based on the order of channel definitions
#'   in the `dynamiteformula` object.
#' @param noncentered \[`logical()`]\cr If `TRUE`, use a noncentered
#'   parameterization for the spline coefficients. Default is `FALSE`. Try
#'   changing this if you encounter divergences or other problems in sampling
#'   for example when simulating from prior predictive distribution.
#'   Can be a single logical value, or vector of logical values, defining the
#'   parameterization separately for each channel, even for channels without
#'   varying effects.
#' @param override \[`logical(1)`]\cr If `FALSE` (the default), an existing
#'    definition for the splines will not be overridden by another call to
#'    `splines()`. If `TRUE`, any existing definitions will be replaced.
#' @return An object of class `splines`.
#' @examples
#' # Two channel model with varying effects, with explicit lower bounds for the
#' # random walk prior standard deviations, with noncentered parameterization
#' # for the first channel and centered for the second channel.
#' obs(y ~ 1, family = "gaussian") + obs(x ~ 1, family = "gaussian") +
#'   lags(type = "varying") +
#'   splines(
#'     df = 20, degree = 3, lb_tau = c(0, 0.1),
#'     noncentered = c(TRUE, FALSE)
#'   )
#'
splines <- function(df = NULL, degree = 3L, lb_tau = 0,
                    noncentered = FALSE, override = FALSE) {
  shrinkage <- FALSE
  stopifnot_(
    checkmate::test_flag(x = shrinkage),
    "Argument {.arg shrinkage} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_flag(x = override),
    "Argument {.arg override} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_int(x = df, lower = 1L, null.ok = TRUE),
    "Argument {.arg df} must be a single positive {.cls integer}."
  )
  stopifnot_(
    checkmate::test_int(x = degree, lower = 1L),
    "Argument {.arg degree} must be a single positive {.cls integer}."
  )
  stopifnot_(
    checkmate::test_logical(
      x = noncentered,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg noncentered} must be a {.cls logical} vector."
  )
  stopifnot_(
    checkmate::test_numeric(
      x = lb_tau,
      lower = 0.0,
      min.len = 1L,
      finite = TRUE,
      any.missing = FALSE
    ),
    "Argument {.arg lb_tau} must be a
     {.cls numeric} vector of non-negative values."
  )
  structure(
    list(
      shrinkage = shrinkage,
      lb_tau = lb_tau,
      noncentered = noncentered,
      bs_opts = list(
        df = onlyif(!is.null(df), as.integer(df)),
        degree = as.integer(degree),
        intercept = TRUE
      )
    ),
    override = override,
    class = "splines"
  )
}

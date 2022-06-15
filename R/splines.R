#' Define the B-splines Used for the Time-varying Coefficients of the Model.
#'
#' This function can be used as part of `dynamiteformula` to define the splines
#' for the time-varying coeffients \eqn{\delta}.
#'
#' TODO: Think the naming of this option, we have a "local" shrinkage due to
#' the random walk prior of the spline coefficients, and the this "global" shrinkage over
#' splines. Need to avoid confusion with the traditional global-local shrinkage priories (e.g. horseshoe)
#'
#' @param shrinkage \[`logical(1)`]\cr If `TRUE`, a common global shrinkage
#'   parameter \eqn{\lambda} is used for the splines so that the standard
#'   deviation of the random walk prior is of the spline coefficients is
#'   \eqn{\lambda\tau}. Default is `FALSE`.
#'   TODO: Needs more testing whether this is useful.
#' @param override \[`logical(1)`]\cr If `FALSE` (the default), an existing
#'    definition for the splines will not be overridden by another call to
#'    `splines()`. If `TRUE`, any existing definitions will be replaced.
#' @param lb_tau \[`numeric()`]\cr Hard constraint(s) on the lower bound of the
#'   standard deviation parameters \eqn{\tau} of the random walk priors. Can be
#'   useful in avoiding divergences in some cases. See also `noncentered`
#'   argument. Can be a single positive value, or vector defining the
#'   lower bound separately for each channel.
#' @param noncentered  \[`logical()`]\cr If `TRUE`, use noncentered
#'   parameterization for the spline coefficients. Default is `FALSE`. Try
#'   changing this if you encounter divergences or other problems in sampling.
#'   Can be a single logical value, or vector of logical values, defining the
#'   parameterization separately for each channel.
#' @param df \[`numeric()`]\cr Degree of freedom, i.e. the total number of
#'   spline coefficients. See [splines::bs()]. Note that the knots are always
#'   defined as equidistant sequence on the interval starting from the first
#'   non-fixed time point to the last time point in the data.
#' @param degree \[`integer(1)`]\cr See [splines::bs()].
#' @export
#' @srrstats {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstats {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @srrstats {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#'
#' TODO explain the fixed time points above in more detail (here or
#' somewhere else)?
splines <- function(shrinkage = FALSE, override = FALSE,
                    df = NULL, degree = 3, lb_tau = 0, noncentered = FALSE) {
  shrinkage <- try_(shrinkage, type = "logical")[1] # TODO check the length
  override <- try_(override, type = "logical")[1] # TODO
  df <- try_(df, type = "integer")[1]
  degree <- try_(degree, type = "integer")
  # length of these this argument can be > 1 for channel-wise definitions,
  # check the length later
  noncentered <- try_(noncentered, type = "logical")
  lb_tau <- try_(lb_tau, type = "numeric")
  # TODO: better error message
  if (any(lb_tau < 0)) {
    stop_("Lower bound for `tau` should be non-negative.")
  }
  structure(
    list(
      shrinkage = shrinkage,
      lb_tau = lb_tau,
      noncentered = noncentered,
      bs_opts = list(
        df = df,
        degree = degree,
        intercept = TRUE
      )
    ),
    override = override,
    class = "splines"
  )
}

#' Checks if the argument represents a splines definition
#'
#' @param x An \R object
#'
#' @noRd
is.splines <- function(x) {
  inherits(x, "splines")
}

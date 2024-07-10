#' Define a Common Latent Factor for the Dynamite Model.
#'
#' This function can be used as part of [dynamite::dynamiteformula()] to define
#' a common latent factor component. The latent factor is modeled as a spline
#' similarly as a time-varying intercept, but instead of having equal effect on
#' each group, there is an additional loading variable for each group so that
#' in the linear predictor we have a term \eqn{\lambda_i \psi_t} for each
#' group \eqn{i}.
#' @export
#' @family formulas
#' @param responses \[`character()`]\cr Names of the responses that the
#'   factor should affect. Default is all responses defined with
#'   `obs` except categorical responses, which do not (yet) support the factor
#'   component.
#' @param nonzero_lambda \[`logical()`]\cr If `TRUE` (the default), assumes
#'   that the mean of factor loadings is nonzero or not. Should be a logical
#'   vector matching the length of `responses` or a single logical value in
#'   case `responses` is `NULL`. See details.
#' @param correlated \[`logical()`]\cr If `TRUE` (the default), the latent
#'   factors are assumed to be correlated between channels.
#' @param noncentered_psi \[`logical(1)`]\cr If `TRUE`, uses a
#'   noncentered parametrization for spline coefficients of all the factors.
#'   The number of knots is based `splines()` call. Default is `FALSE`.
#' @param flip_sign \[`logical(1)`]\cr If `TRUE` (default), try to avoid
#'   multimodality due to sign-switching by defining the sign of \eqn{\lambda}
#'   and \eqn{\psi} based on the mean of \eqn{\omega_1,\ldots, \omega_D}
#'   coefficients. This only affects channels with `nonzero_lambda = FALSE`.
#'   If the true mean of \eqn{\omega}s is close to zero, this might not help,
#'   in which case it is better to set `flip_sign = FALSE` and post-process the
#'   samples in other ways (or use only one chain and/or suitable initial
#'   values). This argument is common to all factors.
#' @return An object of class `latent_factor`.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' # three channel model with common factor affecting for responses x and y
#' obs(y ~ 1, family = "gaussian") +
#'   obs(x ~ 1, family = "poisson") +
#'   obs(z ~ 1, family = "gaussian") +
#'   lfactor(
#'     responses = c("y", "x"), nonzero_lambda = c(TRUE, FALSE),
#'     correlated = TRUE, noncentered_psi = FALSE
#'   )
#'
lfactor <- function(responses = NULL, nonzero_lambda = TRUE, correlated = TRUE,
                    noncentered_psi = FALSE, flip_sign = TRUE) {
  stopifnot_(
    checkmate::test_character(x = responses, min.len = 1L, null.ok = TRUE),
    "Argument {.arg responses} must be a {.cls character} vector."
  )
  stopifnot_(
    checkmate::test_logical(
      x = noncentered_psi,
      any.missing = FALSE,
      len = 1
    ),
    "Argument {.arg noncentered_psi} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_logical(
      x = flip_sign,
      any.missing = FALSE,
      len = 1
    ),
    "Argument {.arg flip_sign} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_logical(
      x = nonzero_lambda,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg nonzero_lambda} must be a {.cls logical} vector."
  )
  stopifnot_(
    checkmate::test_logical(
      x = correlated,
      any.missing = FALSE,
      len = 1
    ),
    "Argument {.arg correlated} must be a single {.cls logical} value."
  )
  structure(
    list(
      responses = responses,
      nonzero_lambda = nonzero_lambda,
      correlated = correlated,
      noncentered_psi = noncentered_psi,
      flip_sign = flip_sign
    ),
    class = "latent_factor"
  )
}

#' Define a Common Latent Factor for the Dynamite Model.
#'
#' This function can be used as part of [dynamite::dynamiteformula()] to define
#' a common latent factor component. The latent factor is modeled as a spline
#' similarly as a time-varying intercept, but instead of having equal effect on
#' each group, there is an additional loading variable for each group so that
#' in the linear predictor we have a term \eqn{\lambda_i \psi_t} for each
#' group \eqn{i}. In order to keep the full factor loadings \eqn{\lambda},
#' the latent factor \eqn{\psi} and the full model identifiable, some
#' restrictions are added to the model. Details will be available in an
#' upcoming paper. This component should be treated as experimental feature.
#'
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
#'   The number of knots is based `splines()` call.
#' @return An object of class `latent_factor`.
#' @examples
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
  noncentered_psi = FALSE) {
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
      noncentered_psi = noncentered_psi
    ),
    class = "latent_factor"
  )
}

#' Define a Common Latent Factor for the Dynamite Model.
#'
#' This function can be used as part of [dynamite::dynamiteformula()] to define
#' a common latent factor component. The latent factor is modeled as a spline
#' similarly as time-varying intercept, but instead of having equal effect on
#' each group, there is additional loading variable for each group so that in
#' the linear predictor we have a term \eqn{\lambda_i  \psi_t} for each group i.
#' In order to keep the full the factor loadings \eqn{\lambda}, the latent
#' factor \eqn{\psi} and the full model identifiable, some restrictions are
#' added to the model. Details will be available in an upcoming paper.
#'
#' @export
#' @param responses \[`character()`]\cr Names of the responses for which the
#'   factor should affect. Default is all responses defined with
#'   `obs` except categorical response, which does not (yet) support factor
#'   component.
#' @param noncentered_lambda \[`logical()`]\cr If `TRUE` (the default), use a
#'   noncentered parametrization for factor loadings. Should be a logical
#'   vector matching the length of `responses` or a single logical value in case
#'   `responses` is `NULL`. Try
#'   changing this if you encounter divergences or other problems in sampling.
#'   Use `splines()` to define whether the spline coefficients of the the
#'   factors are should be centered or not.
#' @param noncentered_psi \[`logical(1)`]\cr If `TRUE`, uses a
#'   noncentered parametrization for spline coefficients of all the factors.
#'   The number of knots is based `splines()` call.
#' @param nonzero_lambda \[`logical()`]\cr If `TRUE` (the default), assumes
#'   that the mean of factor loadings is nonzero or not. Should be a logical
#'   vector matching the length of `responses` or a single logical value in
#'   case `responses` is `NULL`. See details.
#' @param correlated \[`logical()`]\cr If `TRUE` (the default), the latent
#'   factors are assumed to be correlated between channels.
#' @return An object of class `latent_factor`.
#' @examples
#' # three channel model with common factor affecting for responses x and y
#' obs(y ~ 1, family = "gaussian") +
#'   obs(x ~ 1, family = "poisson") +
#'   obs(z ~ 1, family = "gaussian") +
#'   lfactor(responses = c("y", "x"), noncentered_lambda = c(FALSE, TRUE),
#'     noncentered_psi = FALSE, nonzero_lambda = c(TRUE, FALSE))
#'
lfactor <- function(responses = NULL, noncentered_lambda = TRUE,
    noncentered_psi = FALSE, nonzero_lambda = TRUE, correlated = TRUE) {
  stopifnot_(
    checkmate::test_character(x = responses, min.len = 1L, null.ok = TRUE),
    "Argument {.arg responses} must be a {.cls character} vector."
  )
  stopifnot_(
    checkmate::test_logical(
      x = noncentered_lambda,
      any.missing = FALSE,
      len = ifelse_(is.null(responses), 1L, length(responses))
    ),
    "Argument {.arg noncentered_lambda} must be a {.cls logical} vector."
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
      len = ifelse_(is.null(responses), 1L, length(responses))
    ),
    "Argument {.arg nonzero_lambda} must be a {.cls logical} vector."
  )
  stopifnot_(
    checkmate::test_logical(
      x = correlated,
      any.missing = FALSE,
      len = 1
    ),
    "Argument {.arg noncentered_psi} must be a single {.cls logical} value."
  )
  structure(
    list(
      responses = responses,
      noncentered_lambda = noncentered_lambda,
      noncentered_psi = noncentered_psi,
      nonzero_lambda = nonzero_lambda,
      correlated = correlated
    ),
    class = "latent_factor"
  )
}

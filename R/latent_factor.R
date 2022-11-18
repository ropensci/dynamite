#' Define Common Latent Factor for the Dynamite Model.
#'
#' This function can be used as part of [dynamiteformula()] to define
#' a common latent factor component. #TODO definition and constraint on lambdas
#'
#'
#'
#' @export
#' @param responses \[`character()`]\cr Names of the responses for which the
#'   factor should affect. Default is all responses defined with
#'   `obs` except categorical response, which does not (yet) support factor
#'   component.
#' @param noncentered_lambda \[`logical()`]\cr If `TRUE` (the default), use a
#'   noncentered parameterization for factor loadings. Should be a logical
#'   vector matching the length of `responses` or a single logical value in case
#'   `responses` is `NULL`. Try
#'   changing this if you encounter divergences or other problems in sampling.
#'   Use `splines()` to define whether the spline coefficients of the the
#'   factors are should be centered or not.
#' @param noncentered_psi \[`logical(1)`]\cr If `TRUE`, uses a
#'   noncentered parameterization for spline coefficients of all the factors.
#'   The number of knots is based `splines()` call.
#' @param nonzero_lambda  \[`logical()`]\cr If `TRUE` (the default), assumes
#'   that the mean of factor loadings is nonzero or not. Should be a logical
#'   vector matching the length of `responses` or a single logical value in case
#'   `responses` is `NULL`. See details.
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

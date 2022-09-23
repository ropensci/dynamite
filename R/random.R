#' Define Random Intercepts for the Dynamite Model.
#'
#' This function can be used as part of [dynamiteformula()] to define
#' (correlated) random intercepts for each group.
#'
#' With a large number of time points these intercepts can become challenging
#' sample with default priors. This is because with large group sizes the
#' group-level intercepts tend to be behave similarly to fixed group-factor
#' variable so the model becomes overparameterized given these and the common
#' intercept term. Another potential cause for sampling problems is relatively
#' large variation in the intercepts (large sigma_nu) compared to the sampling
#' variation (sigma) in the Gaussian case.
#'
#' @export
#' @param responses \[`character()`]\cr Names of the responses for which the
#'   random intercepts should be defined. Default is all responses defined with
#'   `obs`.
#' @param correlated \[`logical(1)`]\cr If `TRUE` (the default), correlations of
#'   intercepts within a group (i.e., between responses) are modeled so that
#'   the intercepts follow a multivariate normal distribution.
#' @param noncentered \[`logical(1)`]\cr If `TRUE` (the default), use a
#'   noncentered parameterization for random intercepts. Try changing this if
#'   you encounter divergences or other problems in sampling.
#' @return An object of class `random`.
#' @examples
#' # three channel model with correlated random effects for responses x and y
#' obs(y ~ 1, family = "gaussian") +
#'   obs(x ~ 1, family = "poisson") +
#'   obs(z ~ 1, family = "gaussian") +
#'   random(responses = c("y", "x"), correlated = TRUE)
#'
random <- function(responses = NULL, correlated = TRUE, noncentered = TRUE) {
  stopifnot_(
    checkmate::test_character(x = responses, min.len = 1L, null.ok = TRUE),
    "Argument {.arg responses} must be a {.cls character} vector."
  )
  stopifnot_(
    checkmate::test_flag(x = correlated),
    "Argument {.arg correlated} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_flag(x = noncentered),
    "Argument {.arg noncentered} must be a single {.cls logical} value."
  )
  structure(
    list(
      responses = responses,
      correlated = correlated,
      noncentered = noncentered
    ),
    class = "random"
  )
}

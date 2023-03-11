#' Additional Specifications for the Group-level Random Effects of the DMPM
#'
#' This function can be used as part of [dynamiteformula()] to define
#' whether the group-level random effects should be modeled as correlated or
#' not.
#'
#' With a large number of time points random intercepts can become challenging
#' sample with default priors. This is because with large group sizes the
#' group-level intercepts tend to be behave similarly to fixed group-factor
#' variable so the model becomes overparameterized given these and the common
#' intercept term. Another potential cause for sampling problems is relatively
#' large variation in the intercepts (large sigma_nu) compared to the sampling
#' variation (sigma) in the Gaussian case.
#'
#' @export
#' @family formulas
#' @param correlated \[`logical(1)`]\cr If `TRUE` (the default), correlations
#'   of random effects are modeled as multivariate normal.
#' @param noncentered \[`logical(1)`]\cr If `TRUE` (the default), use a
#'   noncentered parameterization for random effects. Try changing this if
#'   you encounter divergences or other problems in sampling.
#' @return An object of class `random_spec`.
#' @examples
#' # two channel model with correlated random effects for responses x and y
#' obs(y ~ 1 + random(~1), family = "gaussian") +
#'   obs(x ~ 1 + random(~1 + z), family = "poisson") +
#'   random_spec(correlated = TRUE)
#'
random_spec <- function(correlated = TRUE, noncentered = TRUE) {
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
      correlated = correlated,
      noncentered = noncentered
    ),
    class = "random_spec"
  )
}

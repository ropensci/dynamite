#' Define Random Intercepts for the Dynamite Model.
#'
#' This function can be used as part of `dynamiteformula` to define random
#' intercepts for each group.
#'
#' @param responses \[`character()`]\cr Names of the responses for which the
#'   random intercepts should be defined. Default is all responses defined with
#'   `obs`.
#' @param correlated \[`logical(1)`]\cr If `TRUE` (the default), correlations of
#'   intercepts within a group (i.e., between responses) are modeled so that
#'   the intercepts follow a multivariate normal distribution.
#' @return An object of class `random`.
#' @export
#' @examples
#' # three channel model with correlated random effects for responses x and y
#' obs(y ~ 1, family = "gaussian") +
#'   obs(x ~ 1, family = "poisson") +
#'   obs(z ~ 1, family = "gaussian") +
#'   random(responses = c("y", "x"), correlated = TRUE)
#'
random <- function(responses = NULL, correlated = TRUE) {
  stopifnot_(
    checkmate::test_character(x = responses, min.len = 1L, null.ok = TRUE),
    "Argument {.arg responses} must be a {.cls character} vector."
  )
  stopifnot_(
    checkmate::test_flag(x = correlated),
    "Argument {.arg correlated} must be a single {.cls logical} value."
  )
  structure(
    list(
      responses = responses,
      correlated = correlated
    ),
    class = "random"
  )
}

#' Is The Argument a `random` Definition?
#'
#' @param x An \R object.
#' @noRd
is.random <- function(x) {
  inherits(x, "random")
}

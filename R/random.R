#' Define the Random Intercepts for the Dynamite Model.
#'
#' This function can be used as part of `dynamiteformula` to define random
#' intercepts for each group.
#'
#' @param channels \[`character()`]\cr Names of the channels for which the
#'   random intercepts should be defined.
#' @param correlated \[`logical(1)`]\cr If `TRUE` (the default), correlation of
#'   intercepts within a group (i.e. between channels) is modelled
#'   (as multivariate normal).
#' @return An object of class `random`.
#' @export
random <- function(channels, correlated = TRUE) {
  stopifnot_(
    checkmate::test_character(x = channels, min.len = 1L),
    "Argument {.arg channels} must be a {.cls character} vector."
  )
  stopifnot_(
    checkmate::test_flag(x = correlated),
    "Argument {.arg correlated} must be a single {.cls correlated} value."
  )
  structure(
    list(
      channels = channels,
      correlated = correlated),
    class = "random"
  )
}

#' Is The Argument a Random Definition
#'
#' @param x An \R object
#' @noRd
is.random <- function(x) {
  inherits(x, "random")
}

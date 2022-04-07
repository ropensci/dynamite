#' Bayesian Time-Varying Coefficients Model
#'
#' @param formula todo
#' @param data todo
#' @param group todo
#' @param time todo
#' @param ... todo
#'
#' @export
btvcm <- function(formula, data, group, time, ...) {
    x <- btvcmfit(formula, data, group, time, ...)
    x
}

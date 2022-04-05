#' Bayesian Time-Varying Coefficients Model
#'
#' @param formula
#' @param data
#' @param group
#' @param time
#' @param ...
#'
#' @export
btvcm <- function(formula, data, group, time, ...) {
    x <- btvcmfit(formula, data, group, time, ...)
    x
}

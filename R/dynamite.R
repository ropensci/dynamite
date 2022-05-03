#' Bayesian Time-Varying Coefficients Model
#'
#' @param formula todo
#' @param data todo
#' @param group todo
#' @param time todo
#' @param ... todo
#'
#' @export
dynamite <- function(formula, data, group, time, ...) {
    x <- dynamitefit(formula, data, group, time, ...)
    x
}

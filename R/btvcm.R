#' The main function
#' @export
btvcm <- function(formula, data, group, time, ...) {
    x <- btvcmfit(formula, data, group, time, ...)
    x
}

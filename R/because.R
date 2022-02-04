#' The main function
#' @export
because <- function(formula, data, group, time, ...) {
    x <- becausefit(formula, data, group, time, ...)
    x
}

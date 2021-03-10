#' The main function
#' @export
because <- function(formula, data, family, ...) {
    x <- becausefit(formula, data)
    return(x)
}

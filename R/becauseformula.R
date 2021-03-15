#' Model formula for \pkg{because}
#'
#' @export
becauseformula <- function(formula, ...) {
    out <- formula
    class(out) <- c("becauseformula")
    return(out)
}

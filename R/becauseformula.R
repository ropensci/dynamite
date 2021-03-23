#' Model formula for \pkg{because}
#'
#' @export
becauseformula <- function(formula, ...) {
    out <- formula
    class(out) <- c("becauseformula")
    return(out)
}

#' Checks if argument is a \code{becauseformula} object
#'
#' @param x An \R object
#'
#' @export
is.becauseformula <- function(x) {
    inherits(x, "becauseformula")
}


# Checks if argument is a formula
is.formula <- function(x) {
    inherits(x, "formula")
}

# Get the right hand side of a formula
# @param x A formula
rhs <- function(x) {
    x <- as.formula(x)
    if (length(x) == 3L) x[3L] else NULL
}

# Get the left hand side of a formula
# @param x A formula
lhs <- function(x) {
    x <- as.formula(x)
    if (length(x) == 3L) x[2L] else NULL
}

#' @export
`+.becauseformula` <- function(e1, e2) {
    if (is.becauseformula(e1)) {
        out <- plus_becauseformula(e1, e2)
    } else {
        stop_("TODO")
    }
    out
}

# Internal `+.becauseformula` for joining formulas
plus_becauseformula <- function(e1, e2) {
    #TODO
}
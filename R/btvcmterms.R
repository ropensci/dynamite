#' Parse formulas of \pkg{because} models
#'
#' @export
btvcmterms <- function(formula, ...) {
    UseMethod("becauseterms")
}

#'
#'@export
btvcmterms.default <- function(formula, ...) {
    # TODO
}

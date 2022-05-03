#' Parse formulas of \pkg{because} models
#'
#' @export
dynamiteterms <- function(formula, ...) {
    UseMethod("becauseterms")
}

#'
#'@export
dynamiteterms.default <- function(formula, ...) {
    # TODO
}

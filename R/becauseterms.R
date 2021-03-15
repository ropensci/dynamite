#' Parse formulas of \pkg{because} models
#'
#' @export
becauseterms <- function(formula, ...) {
    UseMethod("becauseterms")
}

#'
#'@export
becauseterms.default <- function(formula, ...) {

}

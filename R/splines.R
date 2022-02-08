#' Define the B-splines used for the regression coefficients of the model.
#' @param shrinkage A logical value. If \code{TRUE} (the default) a shrinkage prior is used
#'    for smoothing the splines.
#' @param override A logical value. If \code{FALSE} (the default), an existing
#'    definition for the splines will not be overridden by another call to \code{splines}.
#'    If \code{TRUE}, any existing definitions will be replaced.
#' @export
splines <- function(shrinkage = TRUE, override = FALSE) {
    # TODO more parameters to control the splines?
    shrinkage <- try_(shrinkage, type = "logical")[1]
    override <- try_(override, type = "logical")[1]
    structure(
        list(shrinkage = shrinkage),
        override = override,
        class = "splines"
    )
}

# Checks if the argument represents a splines definition
is.splines <- function(x) {
    inherits(x, "splines")
}

#' Define the B-splines used for the regression coefficients of the model.
#' TODO: Think the naming of this option, we have a "local" shrinkage due to
#' the random walk prior of the spline coefficients, and the this "global" shrinkage over
#' splines. Need to avoid confusion with the traditional global-local shrinkage priories (e.g. horseshoe)
#' @param shrinkage a logical value. If \code{TRUE} (the default) a shrinkage prior is used
#'    for smoothing the splines.
#' @param override a logical value. If \code{FALSE} (the default), an existing
#'    definition for the splines will not be overridden by another call to \code{splines}.
#'    If \code{TRUE}, any existing definitions will be replaced.
#' @param df see [splines::bs()].
#' @param knots see [splines::bs()].
#' @param degree see [splines::bs()].
#' @param intercept see [splines::bs()].
#' @param Boundary.knots see [splines::bs()].
#' @export
splines <- function(shrinkage = TRUE, override = FALSE,
                    df = NULL, knots = NULL, degree = 3, intercept = FALSE, Boundary.knots = NULL) {
    shrinkage <- try_(shrinkage, type = "logical")[1]
    override <- try_(override, type = "logical")[1]
    df <- try_(df, type = "integer")[1]
    knots <- try_(knots, type = "numeric")
    degree <- try_(degree, type = "integer")
    structure(
        list(shrinkage = shrinkage,
             noncentered = FALSE, # TODO: allow user to switch
             bs_opts = list(
                df = df,
                knots = knots,
                degree = degree,
                intercept = intercept,
                Boundary.knots = Boundary.knots
            )
        ),
        override = override,
        class = "splines"
    )
}

# Checks if the argument represents a splines definition
is.splines <- function(x) {
    inherits(x, "splines")
}

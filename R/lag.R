#' Add each lagged response as a predictor to each channel. Overrides all other lag definitions.
#' @param k An integer indicating how many previous observations should be included. Defaults to 1.
#' @export
lags <- function(k = 1L) {
    k <- try_(k, type = "integer")[1]
    structure(
        list(k = k),
        class = "lags"
    )
}

# Checks if the argument represents a lags definition
is.lags <- function(x) {
    inherits(x, "lags")
}

# Checks if a formula term includes a lagged response
is_lagterm <- function(x) {
    # TODO implement
    x
}

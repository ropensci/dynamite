#' Add each lagged response as a predictor to each channel.
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

# Extract variables and shifts of lagged terms of the form lag(var, shift)
extract_lags <- function(x) {
    x_split <- unlist(strsplit(x, "\\*|:"))
    has_lag <- grep("lag\\(.*\\)", x_split, perl = TRUE)
    lag_terms <- rm_ws(x_split[has_lag])
    lag_comp <- regexpr("lag\\((?<var>[^,]+)(?:,(?<k>[0-9]+)){0,1}\\)", lag_terms, perl = TRUE)
    n_terms <- length(lag_terms)
    lag_map <- data.frame(src = lag_terms,
                          var = character(n_terms),
                          k = integer(n_terms))
    for (comp in attr(lag_comp, 'capture.name')) {
        start <- attr(lag_comp, "capture.start")[,comp]
        end <- start + attr(lag_comp, "capture.length")[,comp] - 1
        lag_map[,comp] <- substr(lag_terms, start, end)
    }
    lag_map[!duplicated(lag_map)]
}

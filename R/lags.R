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

# Create a lagged version of a vector
lag_ <- function(x, k = 1) {
    c(x[1], x[-length(x)])
}

# Find lag terms in a vector
find_lags <- function(x, processed = FALSE) {
    if (processed) {
        grepl("I\\(lag_\\(.+\\)\\)", x, perl = TRUE)
    } else {
        grepl("lag\\(.+\\)", x, perl = TRUE)
    }
}

# Extract variables and shifts of lagged terms of the form lag(var, shift)
extract_lags <- function(x) {
    has_lag <- find_lags(x, processed = FALSE)
    lag_terms <- x[has_lag]
    # TODO allow vector k
    lag_comp <- regexpr("(?<src>lag\\(\\s*(?<def>.*?)\\s*(?:,\\s*(?<k>[0-9]+)){0,1}\\s*\\))", lag_terms, perl = TRUE)
    lag_map <- list()
    for (comp in attr(lag_comp, 'capture.name')) {
        start <- attr(lag_comp, "capture.start")[,comp]
        end <- start + attr(lag_comp, "capture.length")[,comp] - 1
        lag_map[[comp]] <- trimws(substr(lag_terms, start, end))
    }
    lag_map <- list2DF(lag_map)
    lag_map$k <- as.integer(lag_map$k)
    lag_map[!duplicated(lag_map),]
}
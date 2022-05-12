#' Add each lagged response as a predictor to each channel.
#' @param k \[`integer(1)`: \sQuote{1}]\cr Indicates how many previous
#'   observations should be included.
#' @param type \[`integer(1)`: \sQuote{"fixed"}]\cr Either
#'   `"fixed"` or `"varying"` which indicates whether the coefficients of the
#'   lag terms should vary in time or not.
#' @export
lags <- function(k = 1L, type = c("fixed", "varying")) {
  type <- match.arg(type)
  k <- try_(k, type = "integer")[1]
  structure(
    list(k = k, type = type),
    class = "lags"
  )
}

#' Check if the argument represents a lags definition
#'
#' @param x An R object
#'
#' @noRd
is.lags <- function(x) {
  inherits(x, "lags")
}

#' Create a lagged version of a vector
#'
#' @param x \[`vector()`]\cr A vector of values.
#' @param k \[`integer(1)`: \sQuote{1}]\cr Number of positions to lag by.
#'
#' @noRd
lag_ <- function(x, k = 1) {
  xlen <- length(x)
  c(rep(NA, k), x[1:(length(x) - k)])
}

#' Find lag terms in a character vector
#'
#' @param x \[`character(1)`]\cr A character vector of length one.
#' @param processed \[`logical()`: \sQuote(FALSE)]\cr If true, assumes that
#'   the character string does not contain 'as is' definitions via `I()`.
#'
#' @noRd
find_lags <- function(x, processed = FALSE) {
  if (processed) {
    grepl("I\\(lag_\\(.+\\)\\)", x, perl = TRUE)
  } else {
    grepl("lag\\(.+\\)", x, perl = TRUE)
  }
}

#' Extract lag definitions
#'
#' Extract variables and shifts of lagged terms of the form lag(var, shift)
#' and return them as a data frame for post processing.
#'
#' @param x \[`character(1)`]\cr A character vector of length one.
#'
#' @noRd
extract_lags <- function(x) {
  has_lag <- find_lags(x, processed = FALSE)
  lag_terms <- x[has_lag]
  # TODO allow vector k
  lag_regex <- gregexec(
    pattern = "(?<src>lag\\(\\s*(?<resp>[^\\+]+?)\\s*(?:,\\s*(?<k>[0-9]+)){0,1}\\s*\\))",
    text = lag_terms,
    perl = TRUE
  )
  lag_map <- list()
  lag_matches <- regmatches(lag_terms, lag_regex)
  if (length(lag_matches) > 0) {
    lag_map <- do.call("cbind", args = lag_matches)
    lag_map <- as.data.frame(t(lag_map)[, -1, drop = FALSE])
    lag_map$k <- as.integer(lag_map$k)
    lag_map$k[is.na(lag_map$k)] <- 1L
    lag_map$present <- TRUE
    lag_map |>
      dplyr::distinct() |>
      dplyr::group_by(.data$resp) |>
      tidyr::complete(k = tidyr::full_seq(c(1, .data$k), 1),
                      fill = list(src = "", present = FALSE)) |>
      dplyr::arrange(.data$resp, .data$k) |>
      dplyr::ungroup()
  } else {
    data.frame()
  }
}

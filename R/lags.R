#' Add each lagged response as a predictor to each channel.
#' @param k \[`integer()`: \sQuote{1}]\cr
#'   If `k` is a single integer, then all lags up to an including `k` of all
#'   channels will be included for each channel. If `k` is an integer vector,
#'   then lags indicated by `k` will be included for each channel.
#' @param type \[`integer(1)`: \sQuote{"fixed"}]\cr Either
#'   `"fixed"` or `"varying"` which indicates whether the coefficients of the
#'   lag terms should vary in time or not.
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @export
lags <- function(k = 1L, type = c("fixed", "varying")) {
  type <- match.arg(type)
  k <- try_(k = k, type = "integer")
  if (length(k) == 1) {
    k <- 1:k
  }
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
lag_ <- function(x, k) {
  lag_idx <- seq_len(length(x) - k)
  out <- x
  out[1:k] <- NA
  out[k + lag_idx] <- x[lag_idx]
  out
}

#' Find lag terms in a character vector
#'
#' @param x \[`character()`]\cr A character vector
#'
#' @noRd
find_lags <- function(x) {
  grepl("lag\\([^\\)]+\\)", x, perl = TRUE)
}

#' Extract non-lag variables
#'
#' @param x \[`character(1)`]\cr A character vector
#'
#' @noRd
extract_nonlags <- function(x) {
  has_lag <- find_lags(x)
  x[!has_lag]
}

#' Extract lag definitions
#'
#' Extract variables and shifts of lagged terms of the form lag(var, k)
#' and return them as a data frame for post processing.
#'
#' @param x \[`character()`]\cr a character vector
#'
#' @noRd
extract_lags <- function(x) {
  has_lag <- find_lags(x)
  lag_terms <- x[has_lag]
  # TODO allow vector k
  lag_regex <- gregexec(
    pattern = paste0(
      "(?<src>lag\\(\\s*(?<var>[^\\+]+?)\\s*",
      "(?:,\\s*(?<k>\\-{0,1}[0-9]+)){0,1}\\s*\\))"
    ),
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
      dplyr::group_by(.data$var) |>
      tidyr::complete(k = tidyr::full_seq(c(1, .data$k), 1),
                      fill = list(src = "", present = FALSE)) |>
      dplyr::arrange(.data$var, .data$k) |>
      dplyr::ungroup()
  } else {
    data.frame(
      src = character(0),
      var = character(0),
      k = integer(0),
      present = logical(0)
    )
  }
}

#' Extract lag shift values of a specific variable from a character string
#'
#' @param x \[`character(1)`]\cr a character vector of length one
#' @param self \[`character(1)`]\cr variable whose lags to look for
#'
#' @noRd
extract_self_lags <- function(x, self) {
  lag_regex <-  gregexec(
    pattern = paste0(
      "lag\\(", self, ",\\s*(?<k>[0-9]+)\\s*\\)"
    ),
    text = x,
    perl = TRUE
  )
  lag_matches <- regmatches(x, lag_regex)[[1]]
  if (length(lag_matches) > 0) {
    as.integer(lag_matches["k", ])
  } else {
    0
  }
}

#' Add Each Lagged Response as a Predictor to Each Channel.
#'
#' @param k \[`integer()`: \sQuote{1}]\cr
#'   Lagged values indicated by `k`  of each observed response variable
#'   will be added as a predictor for each channel.
#' @param type \[`integer(1)`: \sQuote{"fixed"}]\cr Either
#'   `"fixed"` or `"varying"` which indicates whether the coefficients of the
#'   added lag terms should vary in time or not.
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @srrstats {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstats {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @export
lags <- function(k = 1L, type = c("fixed", "varying")) {
  type <- match.arg(type)
  k <- try_type(k, "integer")
  structure(
    list(k = k, type = type),
    class = "lags"
  )
}

#' Is the Argument a `lags` Definition
#'
#' @param x An R object
#' @noRd
is.lags <- function(x) {
  inherits(x, "lags")
}

#' Create a Lagged Version of a Vector
#'
#' @param x \[`vector()`]\cr A vector of values.
#' @param k \[`integer(1)`: \sQuote{1}]\cr Number of positions to lag by.
#' @noRd
lag_ <- function(x, k = 1) {
  lag_idx <- seq_len(length(x) - k)
  out <- x
  out[1:k] <- NA
  out[k + lag_idx] <- x[lag_idx]
  out
}

#' Find Lag Terms in a Character Vector
#'
#' @param x \[`character()`]\cr A character vector.
#' @noRd
find_lags <- function(x) {
  grepl("lag\\([^\\)]+\\)", x, perl = TRUE)
}

#' Extract Non-lag Variables
#'
#' @param x \[`character(1)`]\cr A character vector.
#' @noRd
extract_nonlags <- function(x) {
  has_lag <- find_lags(x)
  x[!has_lag]
}

#' Extract Lag Definitions
#'
#' Extract variables and shifts of lagged terms of the form `lag(var, k)`
#' and return them as a data frame for post processing.
#'
#' @param x \[`character()`]\cr a character vector.
#' @noRd
extract_lags <- function(x) {
  has_lag <- find_lags(x)
  if (any(has_lag)) {
    lag_terms <- paste0(x[has_lag], " ")
  } else {
    lag_terms <- character(0)
  }
  lag_regex <- gregexec(
    pattern = paste0(
      "(?<src>lag\\(\\s*(?<var>[^\\+\\)\\,]+?)\\s*",
      "(?:,\\s*(?<k>.+?)){0,1}\\))\\s"
    ),
    text = lag_terms,
    perl = TRUE
  )
  lag_matches <- regmatches(lag_terms, lag_regex)
  if (length(lag_matches) > 0L) {
    lag_map <- do.call("cbind", args = lag_matches)
    lag_map <- as.data.frame(t(lag_map)[, -1, drop = FALSE])
    n_lag <- nrow(lag_map)
    k_str <- lag_map$k
    lag_map$k <- integer(n_lag)
    lag_map$present <- TRUE
    lag_map$increment <- FALSE
    for (j in seq_len(n_lag)) {
      if (nzchar(k_str[j])) {
        k_expr <- try(str2lang(k_str[j]), silent = TRUE)
        if ("try-error" %in% class(k_expr)) {
          stop_("Invalid shifted value expression {.code {k_str[j]}}.")
        }
        k_coerce <- try(eval(k_expr), silent = TRUE)
        if ("try-error" %in% class(k_coerce)) {
          stop_("Invalid lagged value definition {.code {lag_map$src[j]}}.")
        }
        failed <- FALSE
        k_coerce <- tryCatch(
          expr = as.integer(k_coerce),
          error = function(e) NULL,
          warning = function(w) NULL
        )
        if (is.null(k_coerce) ||
            identical(length(k_coerce), 0L) ||
            any(is.na(k_coerce))) {
          stop_("Unable to coerce lag shift value to {.cls integer} in
                 {.code {lag_map$src[j]}}."
          )
        }
        if (any(k_coerce <= 0L)) {
          stop_(c(
            "Shift values must be positive in {.fun lag}:",
            `x` = "Nonpositive shift value was found in
                   {.code {lag_map$src[j]}}."
          ))
        }
        lag_map$k[j] <- k_coerce[1L]
        if (length(k_coerce) > 1L) {
          new_lags <- data.frame(
            src = lag_map$src[j],
            var = lag_map$var[j],
            k = k_coerce[-1L],
            present = TRUE,
            increment = TRUE
          )
          lag_map <- rbind(lag_map, new_lags)
        }
      } else {
        lag_map$k[j] <- 1L
      }
    }
    lag_map |>
      dplyr::distinct() |>
      dplyr::group_by(.data$var) |>
      tidyr::complete(k = tidyr::full_seq(c(1L, .data$k), 1L),
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

#' Extract Lag Shift Values of a Specific Variable from a Character String
#'
#' @param x \[`character(1)`]\cr a character vector of length one.
#' @param self \[`character(1)`]\cr variable whose lags to look for.
#' @noRd
extract_self_lags <- function(x, self) {
  lag_regex <-  gregexec(
    pattern = paste0(
      "lag\\(", self, ",\\s*(?<k>.+?)\\s*\\)"
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

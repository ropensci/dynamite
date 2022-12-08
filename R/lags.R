#' Add Lagged Responses as Predictors to Each Channel of a Dynamite Model
#'
#' Adds the lagged value of the response of each channel specified via
#' [dynamiteformula()] as a predictor to each channel. The added predictors
#' can be either time-varying or time-invariant.
#'
#' @export
#' @param k \[`integer()`]\cr
#'   Values lagged by `k` units of time of each observed response variable
#'   will be added as a predictor for each channel. Should be a positive
#'   (unrestricted) integer.
#' @param type \[`integer(1)`]\cr Either
#'   `"fixed"` or `"varying"` which indicates whether the coefficients of the
#'   added lag terms should vary in time or not.
#' @return An object of class `lags`.
#' @srrstats {G2.3a} Uses match.arg
#' @examples
#' obs(y ~ -1 + varying(~x), family = "gaussian") +
#'   lags(type = "varying") + splines(df = 20)
#'
#' # A two-channel categorical model with time-invariant predictors
#' # here, lag terms are specified manually
#' obs(x ~ z + lag(x) + lag(y), family = "categorical") +
#'   obs(y ~ z + lag(x) + lag(y), family = "categorical")
#'
#' # The same categorical model as above, but with the lag terms
#' # added using 'lags'
#' obs(x ~ z, family = "categorical") +
#'   obs(y ~ z, family = "categorical") +
#'   lags(type = "fixed")
#'
lags <- function(k = 1L, type = c("fixed", "varying")) {
  type <- onlyif(is.character(type), tolower(type))
  type <- try(match.arg(type, c("fixed", "varying")), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be \"fixed\" or \"varying\"."
  )
  stopifnot_(
    checkmate::test_integerish(
      x = k,
      lower = 1L,
      any.missing = FALSE,
      min.len = 1L,
      unique = TRUE,
    ),
    "Argument {.arg k} must be an {.cls integer} vector with positive values."
  )
  structure(
    list(k = as.integer(k), type = type),
    class = "lags"
  )
}

#' Create a Lagged Version of a Vector
#'
#' @param x \[`vector()`]\cr A vector of values.
#' @param k \[`integer(1)`]\cr Number of positions to lag by.
#' @noRd
lag_ <- function(x, k = 1) {
  lag_idx <- seq_len(length(x) - k)
  out <- x
  out[seq_len(k)] <- NA
  out[k + lag_idx] <- x[lag_idx]
  out
}

#' Adds Default Shift Values to Terms of the Form `lag(y)`
#'
#' @param x A `language` object.
#' @noRd
complete_lags <- function(x) {
  if (identical(length(x), 1L)) {
    return(x)
  }
  if (identical(deparse1(x[[1L]]), "lag")) {
    xlen <- length(x)
    if (identical(xlen, 2L)) {
      x <- str2lang(
        paste0("lag(", deparse1(x[[2L]]), ", ", "1)")
      )
    } else if (identical(xlen, 3L)) {
      k <- verify_lag(x[[3L]], deparse1(x))
      x <- str2lang(
        paste0("lag(", deparse1(x[[2L]]), ", ", k, ")")
      )
    } else {
      stop_(c(
        "Invalid lag definition {.code {deparse1(x)}}:",
        `x` = "Too many arguments supplied to {.fun lag}."
      ))
    }
  } else {
    for (i in seq_along(x)) {
      x[[i]] <- complete_lags(x[[i]])
    }
  }
  x
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
  x[!find_lags(x)]
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
  lag_terms <- ifelse_(
    any(has_lag),
    paste0(x[has_lag], " "),
    character(0)
  )
  lag_regex <- gregexec(
    pattern = paste0(
      "(?<src>lag\\((?<var>[^\\+\\)\\,]+?)",
      "(?:,\\s*(?<k>[0-9]+)){0,1}\\))"
    ),
    text = lag_terms,
    perl = TRUE
  )
  lag_matches <- regmatches(lag_terms, lag_regex)
  if (length(lag_matches) > 0L) {
    lag_map <- do.call("cbind", args = lag_matches)
    lag_map <- as.data.frame(t(lag_map)[, -1L, drop = FALSE])
    lag_map$k <- as.integer(lag_map$k)
    lag_map$k[is.na(lag_map$k)] <- 1L
    lag_map$present <- TRUE
    lag_map <- unique(lag_map)
    lag_var <- sort(unique(lag_map$var))
    expanded <- vector(mode = "list", length = length(lag_var))
    for (i in seq_along(lag_var)) {
      v <- lag_var[i]
      tmp <- lag_map[lag_map$var == v, ]
      tmp <- tmp[order(tmp$k), ]
      full <- data.frame(
        src = "",
        var = v,
        k = seq.int(1L, max(tmp$k)),
        present  = FALSE
      )
      expanded[[i]] <- full[full$k %in% tmp$k, ] <- tmp
    }
    lag_map <- do.call("rbind", args = expanded)
  } else {
    data.frame(
      src = character(0L),
      var = character(0L),
      k = integer(0L),
      present = logical(0L)
    )
  }
}

#' Verify that `k` in `lag(y, k)` Represents a Valid Shift Value Expression
#'
#' @param k \[`language`]\cr The shift value definition.
#' @param lag_str \[`character(1)`]\cr The full lag term definition.
#' @noRd
verify_lag <- function(k, lag_str) {
  k_str <- deparse1(k)
  k_coerce <- try(eval(k), silent = TRUE)
  stopifnot_(
    !inherits(k_coerce, "try-error"),
    "Invalid shift value expression {.code {k_str}}."
  )
  k_coerce <- tryCatch(
    expr = as.integer(k_coerce),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  k_len <- length(k_coerce)
  stopifnot_(
    !is.null(k_coerce) && !identical(k_len, 0L) && !any(is.na(k_coerce)),
    "Unable to coerce shift value to {.cls integer} in {.code {lag_str}}."
  )
  stopifnot_(
    identical(k_len, 1L),
    c(
      "Shift value must be a single {.cls integer} in {.fun lag}:",
      `x` = "Multiple shift values were found in {.code {lag_str}}."
    )
  )
  stopifnot_(
    k_coerce > 0L,
    c(
      "Shift value must be positive in {.fun lag}:",
      `x` = "Nonpositive shift value was found in {.code {lag_str}}."
    )
  )
  k_coerce
}

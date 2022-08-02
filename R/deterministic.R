#' Initialize Deterministic Channels
#'
#' Creates new columns in the data for each deterministic channel and ensures
#' that the column types are correct (e.g., a lagged value has the same type
#' as its "parent").
#'
#' @param data \[`data.table`]\cr Data to assign initial values into.
#' @param dd  \[`dynamiteformula`]\cr Formula of deterministic channels.
#' @param dlp \[`dynamiteformula`]\cr Formula of lagged predictors.
#' @param dld \[`dynamiteformula`]\cr Formula of lagged deterministic channels.
#' @param dls \[`dynamiteformula`]\cr Formula of lagged stochastic channels.
#' @noRd
initialize_deterministic <- function(data, dd, dlp, dld, dls) {
  resp_pred <- attr(dlp, "original_response")
  for (i in seq_along(dlp)) {
    data[ ,(dlp[[i]]$response) := data[[resp_pred[i]]]]
    data[ ,(dlp[[i]]$response) := NA]
  }
  rhs_stoch <- get_predictors(dls)
  for (i in seq_along(dls)) {
    stopifnot_(
      rhs_stoch[i] %in% names(data),
      "Can't find variable{?s} {.var {rhs_stoch[i]}} in {.arg data}."
    )
    data[ ,(dls[[i]]$response) := data[[rhs_stoch[i]]]]
    data[ ,(dls[[i]]$response) := NA]
  }
  for (i in seq_along(dd)) {
    as_fun <- paste0("as.", dd[[i]]$specials$resp_type)
    past <- do.call(as_fun, args = list(0))
    data[, dd[[i]]$response := past]
    data[, dd[[i]]$response := NA]
  }
  rhs_det <- get_predictors(dld)
  ro_det <- attr(dld, "rank_order")
  init <- has_past(dld)
  for (k in ro_det) {
    if (init[k]) {
      as_fun <- paste0("as.", dld[[k]]$specials$resp_type)
      past <- do.call(as_fun, args = list(dld[[k]]$specials$past))
      data[ ,(dld[[k]]$response) := past]
      data[ ,(dld[[k]]$response) := NA]
    } else {
      data[ ,(dld[[k]]$response) := data[[rhs_det[k]]]]
      data[ ,(dld[[k]]$response) := NA]
    }
  }
  if (length(dd) > 0L) {
    cl <- get_quoted(dd)
    res <- try(assign_deterministic(data, cl, 1L), silent = TRUE)
    stopifnot_(
      !"try-error" %in% class(res),
      c(
        "Unable to evaluate definitions of deterministic channels:",
        `i` = "Some variables are possibly missing or incorrect."
      )
    )
  }
}

#' Assign Initial Values of Deterministic Channels
#'
#' @inheritParams initialize_deterministic
#' @param idx \[`integer()`] Vector of indices.
#' @param fixed \[`integer(1)`] Number of fixed time points.
#' @param group_var \[`character(1)`] Grouping variable name.
#' @noRd
assign_initial_values <- function(data, dd, dlp, dld, dls,
                                  idx, fixed, group_var) {
  n_det <- length(dd)
  n_lag_pred <- length(dlp)
  n_lag_det <- length(dld)
  n_lag_stoch <- length(dls)
  ro_pred <- ifelse_(
    is.null(attr(dlp, "rank_order")),
    integer(0L),
    attr(dlp, "rank_order")
  )
  ro_det <- ifelse_(
    is.null(attr(dld, "rank_order")),
    integer(0L),
    attr(dld, "rank_order")
  )
  ro_stoch <- seq_len(n_lag_stoch)
  resp_pred <- attr(dlp, "original_response")
  lhs_det <- get_responses(dld)
  rhs_det <- get_predictors(dld)
  lhs_stoch <- get_responses(dls)
  rhs_stoch <- get_predictors(dls)
  cl <- get_quoted(dd)
  for (k in ro_pred) {
    data[, (dlp[[k]]$response) := lapply(.SD, lag_, k = 1L),
         .SDcols = resp_pred[k], by = group_var]
  }
  skip <- TRUE # skip lags at first index to avoid overriding initial values
  for (i in seq_len(fixed + 1L)) {
    idx <- idx + 1L
    assign_lags(data, ro_det, idx, lhs_det, rhs_det, skip)
    assign_lags(data, ro_stoch, idx, lhs_stoch, rhs_stoch, skip)
    assign_deterministic(data, cl, idx)
    skip <- FALSE
  }
  init <- has_past(dld)
  for (i in seq_along(dld)) {
    if (init[i]) {
      as_fun <- paste0("as.", dld[[i]]$specials$resp_type)
      past <- do.call(as_fun, args = list(dld[[i]]$specials$past))
      data[idx, (lhs_det[i]) := past]
    }
  }
  assign_deterministic(data, cl, idx)
}

#' Evaluate Definitions and Assign Values of Deterministic Channels
#'
#' @param data \[`data.table`]\cr Data table to assign the values into.
#' @param cl \[`language`]\cr A quoted expression defining the channels.
#' @param idx \[`integer()`]\cr A vector of indices to assign values into.
#' @noRd
assign_deterministic <- function(data, cl, idx) {
  # This should work in the development version of data.table
  data[idx, cl, env = list(cl = cl)]
}

#' Assign Values of Lagged Channels
#'
#' @param data \[`data.table`]\cr Data to assign the values into.
#' @param ro \[`integer`]\cr Rank order of the channels (the evaluation order).
#' @param idx \[`integer()`]\cr A vector of indices.
#' @param lhs \[`character()`]\cr The lagged outcome variable names.
#' @param rhs \[`character()`]\cr The names of the variables being lagged.
#' @param skip \[`logical(1)`]\cr Skip evaluation on this iteration.
#' @param offset \[`integer(1)`: \sQuote{1L}] The distance between consequent
#'   observations in `data`.
#' @noRd
assign_lags <- function(data, ro, idx, lhs, rhs, skip = FALSE, offset = 1L) {
  if (!skip) {
    for (k in ro) {
      set(data, i = idx, j = lhs[k], value = data[[rhs[k]]][idx - offset])
    }
  }
}

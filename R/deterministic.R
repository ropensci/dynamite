#' Initialize Deterministic Channels
#'
#' Creates new columns in the data for each deterministic channel and ensures
#' that the column types are correct (e.g., a lagged value has the same type
#' as its "parent").
#'
#' @param data \[`data.table()`]\cr Data to assign initial values into.
#' @param dd  \[`dynamiteformula`]\cr Formula of deterministic channels.
#' @param dlp \[`dynamiteformula`]\cr Formula of lagged predictors.
#' @param dld \[`dynamiteformula`]\cr Formula of lagged deterministic channels.
#' @param dls \[`dynamiteformula`]\cr Formula of lagged stochastic channels.
#' @noRd
initialize_deterministic <- function(data, dd, dlp, dld, dls) {
  resp_pred <- attr(dlp, "original_response")
  for (i in seq_along(dlp)) {
    data[, (dlp[[i]]$response) := data[[resp_pred[i]]]]
    data[, (dlp[[i]]$response) := NA]
  }
  rhs_ls <- get_predictors(dls)
  for (i in seq_along(dls)) {
    stopifnot_(
      rhs_ls[i] %in% names(data),
      "Can't find variable{?s} {.var {rhs_ls[i]}} in {.arg data}."
    )
    data[, (dls[[i]]$response) := data[[rhs_ls[i]]]]
    data[, (dls[[i]]$response) := NA]
  }
  for (i in seq_along(dd)) {
    as_fun <- paste0("as.", dd[[i]]$specials$resp_type)
    past <- do.call(as_fun, args = list(0))
    data[, dd[[i]]$response := past]
    data[, dd[[i]]$response := NA]
  }
  rhs_ld <- get_predictors(dld)
  ro_ld <- attr(dld, "rank_order")
  init <- has_past(dld)
  for (k in ro_ld) {
    if (init[k]) {
      as_fun <- paste0("as.", dld[[k]]$specials$resp_type)
      past <- do.call(as_fun, args = list(dld[[k]]$specials$past[1L]))
      data[, (dld[[k]]$response) := past]
      data[, (dld[[k]]$response) := NA]
    } else {
      data[, (dld[[k]]$response) := data[[rhs_ld[k]]]]
      data[, (dld[[k]]$response) := NA]
    }
  }
  if (length(dd) > 0L) {
    cl <- get_quoted(dd)
    res <- try(assign_deterministic(data, 1L, cl), silent = TRUE)
    stopifnot_(
      !inherits(res, "try-error"),
      c(
        "Unable to evaluate definitions of deterministic channels:",
        `x` = attr(res, "condition")$message
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
assign_initial_values <- function(data, idx, dd, dlp, dld, dls,
                                  fixed, group_var) {
  ro_lp <- ifelse_(
    is.null(attr(dlp, "rank_order")),
    integer(0L),
    attr(dlp, "rank_order")
  )
  ro_ld <- ifelse_(
    is.null(attr(dld, "rank_order")),
    integer(0L),
    attr(dld, "rank_order")
  )
  ro_ls <- seq_along(dls)
  resp_lp <- attr(dlp, "original_response")
  k_lp <- attr(dlp, "original_shift")
  lhs_ld <- get_responses(dld)
  rhs_ld <- get_predictors(dld)
  lhs_ls <- get_responses(dls)
  rhs_ls <- get_predictors(dls)
  cl <- get_quoted(dd)
  # avoid NSE notes in R CMD check
  ..i <- NULL
  ..k_lp <- NULL
  for (i in ro_lp) {
    data[, (dlp[[i]]$response) := lapply(.SD, lag_, ..k_lp[..i]),
      .SDcols = resp_lp[i], by = group_var
    ]
  }
  init <- which(has_past(dld))
  for (i in init) {
    as_fun <- paste0("as.", dld[[i]]$specials$resp_type)
    past <- do.call(as_fun, args = list(dld[[i]]$specials$past))
    data[, (lhs_ld[i]) := past]
  }
  idx <- idx + 1L
  assign_deterministic(data, idx, cl)
  for (i in seq_len(fixed)) {
    idx <- idx + 1L
    assign_lags_init(data, idx, ro_ld, lhs_ld, rhs_ld)
    assign_lags_init(data, idx, ro_ls, lhs_ls, rhs_ls)
    assign_deterministic(data, idx, cl)
  }
  assign_deterministic(data, idx, cl)
}

#' Evaluate Definitions and Assign Values of Deterministic Channels
#'
#' @param data \[`data.table`]\cr Data table to assign the values into.
#' @param idx \[`integer()`]\cr A vector of indices to assign values into.
#' @param cl \[`language`]\cr A quoted expression defining the channels.
#' @noRd
assign_deterministic <- function(data, idx, cl) {
  # Requires data.table version 1.14.3 or higher
  data[idx, cl, env = list(cl = cl)]
}

#' Assign Values of Lagged Channels
#'
#' @param data \[`data.table`]\cr Data to assign the values into.
#' @param idx \[`integer()`]\cr A vector of indices.
#' @param ro \[`integer`]\cr Rank order of the channels (the evaluation order).
#' @param lhs \[`character()`]\cr The lagged outcome variable names.
#' @param rhs \[`character()`]\cr The names of the variables being lagged.
#' @param skip \[`logical(1)`]\cr Skip evaluation on this iteration.
#' @param offset \[`integer(1)`]\cr The distance between consequent
#'   observations in `data`.
#' @noRd
assign_lags <- function(data, idx, ro, lhs, rhs, skip = FALSE, offset = 1L) {
  if (!skip) {
    for (k in ro) {
      data.table::set(
        data, i = idx, j = lhs[k], value = data[[rhs[k]]][idx - offset]
      )
    }
  }
}

#' Assign Values of Lagged Channels Without Overriding Initial Values
#'
#' @inheritParams assign_lags
#' @noRd
assign_lags_init <- function(data, idx, ro, lhs, rhs, offset = 1L) {
  for (k in ro) {
    val <- data[[rhs[k]]][idx - offset]
    na_val <- is.na(val)
    val[na_val] <- data[[lhs[k]]][idx][na_val]
    data.table::set(data, i = idx, j = lhs[k], value = val)
  }
}

#' Evaluate Definitions of Deterministic Channels
#'
#' @inheritParams parse_data
#' @param dformulas \[list()]\cr The return object of [dynamite::parse_lags()].
#' @noRd
evaluate_deterministic <- function(dformulas, data, group_var, time_var) {
  fixed <- as.integer(attr(dformulas$all, "max_lag"))
  n_time <- n_unique(data[[time_var]])
  dd <- dformulas$det
  dlp <- dformulas$lag_pred
  dld <- dformulas$lag_det
  dls <- dformulas$lag_stoch
  cl <- get_quoted(dd)
  initialize_deterministic(data, dd, dlp, dld, dls)
  idx <- seq.int(1L, nrow(data), by = n_time) - 1L
  assign_initial_values(data, idx, dd, dlp, dld, dls, fixed, group_var)
  if (n_time > fixed + 1L) {
    ro_ld <- ifelse_(
      is.null(attr(dld, "rank_order")),
      integer(0L),
      attr(dld, "rank_order")
    )
    ro_ls <- seq_along(dls)
    lhs_ld <- get_responses(dld)
    rhs_ld <- get_predictors(dld)
    lhs_ls <- get_responses(dls)
    rhs_ls <- get_predictors(dls)
    idx <- idx + fixed + 1L
    for (i in seq.int(fixed + 2L, n_time)) {
      idx <- idx + 1L
      assign_lags(data, idx, ro_ld, lhs_ld, rhs_ld)
      assign_lags(data, idx, ro_ls, lhs_ls, rhs_ls)
      assign_deterministic(data, idx, cl)
    }
  }
}

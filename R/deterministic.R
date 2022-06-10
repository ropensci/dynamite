# TODO documentation
initialize_deterministic <- function(data, dd, dlp, dld, dls) {
  resp_pred <- attr(dlp, "original_response")
  for (i in seq_along(dlp)) {
    data[ ,(dlp[[i]]$response) := data[[resp_pred[i]]]]
    data[ ,(dlp[[i]]$response) := NA]
  }
  rhs_stoch <- get_predictors(dls)
  for (i in seq_along(dls)) {
    data[ ,(dls[[i]]$response) := data[[rhs_stoch[i]]]]
    data[ ,(dls[[i]]$response) := NA]
  }
  n_det <- length(dd)
  if (n_det > 0) {
    for (i in seq_len(n_det)) {
      as_fun <- paste0("as.", dd[[i]]$specials$resp_type)
      past <- do.call(as_fun, args = list(0))
      data[, dd[[i]]$response := past]
      data[, dd[[i]]$response := NA]
    }
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
  if (n_det > 0) {
    cl <- get_quoted(dd)
    res <- try(assign_deterministic(data, cl, 1L), silent = TRUE)
    if ("try-error" %in% class(res)) {
      stop_(c(
        "Unable to evaluate definitions of deterministic channels:",
        `i` = "Some variables are possibly missing or incorrect."
      ))
    }
  }
}

assign_initial_values <- function(data, dd, dlp, dld, dls,
                                  idx, fixed, group_var) {
  n_det <- length(dd)
  n_lag_pred <- length(dlp)
  n_lag_det <- length(dld)
  n_lag_stoch <- length(dls)
  ro_pred <- attr(dlp, "rank_order")
  ro_det <- attr(dld, "rank_order")
  ro_stoch <- 1L:n_lag_stoch
  resp_pred <- attr(dlp, "original_response")
  lhs_det <- get_responses(dld)
  rhs_det <- get_predictors(dld)
  lhs_stoch <- get_responses(dls)
  rhs_stoch <- get_predictors(dls)
  cl <- get_quoted(dd)
  for (k in ro_pred) {
    if (is.null(group_var)) {
      data[, (dlp[[k]]$response) := lapply(.SD, lag_, k = 1L),
           .SDcols = resp_pred[k]]
    } else {
      data[, (dlp[[k]]$response) := lapply(.SD, lag_, k = 1L),
           .SDcols = resp_pred[k], by = group_var]
    }
  }
  for (i in seq_len(fixed + 1L)) {
    idx <- idx + 1L
    if (n_lag_det > 0 && i > 1) {
      assign_lags(data, ro_det, idx, lhs_det, rhs_det)
    }
    if (n_lag_stoch > 0 && i > 1) {
      assign_lags(data, ro_stoch, idx, lhs_stoch, rhs_stoch)
    }
    if (n_det > 0) {
      assign_deterministic(data, cl, idx)
    }
  }
  if (length(dld) > 0) {
    init <- has_past(dld)
    for (i in seq_along(dld)) {
      if (init[i]) {
        as_fun <- paste0("as.", dld[[i]]$specials$resp_type)
        past <- do.call(as_fun, args = list(dld[[i]]$specials$past))
        data[idx, (lhs_det[i]) := past]
      }
    }
  }
}

assign_deterministic <- function(data, cl, idx) {
  # This should work in the next version of data.table
  data[idx, cl, env = list(cl = cl)]
  invisible(NULL)
}

assign_lags <- function(data, ro, idx, lhs, rhs) {
  for (k in ro) {
    set(data, i = idx, j = lhs[k], value = data[idx - 1][[rhs[k]]])
  }
}

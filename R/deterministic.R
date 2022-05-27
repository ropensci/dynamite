# TODO documentation
assign_initial_values <- function(data, dd, dl, idx) {
  resp_det <- get_responses(dd)
  init <- has_past(dd)
  if (any(init)) {
    init_idx <- which(init)
    for (i in init_idx) {
      data[idx, (resp_det[i]) := dd[[i]]$specials$past]
    }
  }
  if (length(dl) > 0) {
    ro <- attr(dl, "rank_order")
    init <- has_past(dl)
    lag_lhs <- get_responses(dl)
    lag_rhs <- get_predictors(dl)
    for (i in seq_along(dl)) {
      if (init[i]) {
        data[idx, (lag_lhs[i]) := dl[[i]]$specials$past]
      } else {
        data[idx, (lag_lhs[i]) := data[idx,][[lag_rhs[i]]]]
        data[idx, (lag_lhs[i]) := NA]
      }
    }
  }
}

assign_deterministic <- function(data, cl, idx) {
  # This should work in the next version of data.table
  data[idx, cl, env = list(cl = cl)]
}

assign_lags <- function(data, ro, idx, lag_lhs, lag_rhs) {
  for (k in ro) {
    set(data, i = idx, j = lag_lhs[k],
        value = data[idx - 1][[lag_rhs[k]]])
  }
}

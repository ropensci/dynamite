# TODO documentation
assign_initial_values <- function(data, dd, dl, idx) {
  resp_det <- get_responses(dd)
  for (i in seq_along(dd)) {
    data[idx, (resp_det[i]) := NA_real_]
  }
  if (length(dl) > 0) {
    init <- has_past(dl)
    lag_lhs <- get_responses(dl)
    lag_rhs <- get_predictors(dl)
    for (i in seq_along(dl)) {
      if (init[i]) {
        fixed <- dl[[i]]$specials$past_offset
        idx_fixed <- idx + fixed
        data[idx_fixed, (lag_lhs[i]) := dl[[i]]$specials$past]
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

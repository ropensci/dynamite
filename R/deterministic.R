# TODO documentation
assign_initial_values <- function(data, dd, dld, dls, idx) {
  resp_det <- get_responses(dd)
  for (i in seq_along(dd)) {
    # TODO force type definition in aux
    as_fun <- paste0("as.", dd[[i]]$specials$resp_type)
    data[idx, (resp_det[i]) := NA]
    data[, (resp_det[i]) := lapply(.SD, as_fun), .SDcols = resp_det[i]]
  }
  if (length(dls) > 0) {
    lhs <- get_responses(dls)
    rhs <- get_predictors(dls)
    for (i in seq_along(dls)) {
      data[idx, (lhs[i]) := data[idx,][[rhs[i]]]]
      data[idx, (lhs[i]) := NA]
    }
  }
  if (length(dld) > 0) {
    init <- has_past(dld)
    lhs <- get_responses(dld)
    rhs <- get_predictors(dld)
    for (i in seq_along(dld)) {
      if (init[i]) {
        as_fun <- paste0("as.", dld[[i]]$specials$resp_type)
        idx_fixed <- idx + dld[[i]]$specials$past_offset
        past <- do.call(as_fun, args = list(dld[[i]]$specials$past))
        data[idx_fixed, (lhs[i]) := past]
      }# else {
      #  data[idx, (lhs[i]) := data[idx,][[rhs[i]]]]
      #  data[idx, (lhs[i]) := NA]
      #}
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

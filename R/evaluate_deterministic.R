# TODO documentation
evaluate_deterministic <- function(data, dformula, group_var, time_var) {

  n_time <- length(unique(data[[time_var]]))
  n_id <- length(unique(data[[group_var]]))
  idx_base <- seq(0, n_time * (n_id - 1), by = n_time)
  assign_initial_values(data, dformula, group_var)
  if (length(dformula) > 0 && n_time > 1) {
    exprs <- lapply(get_formulas(dformula), function(x) {
      formula2expression(formula_rhs(x))
    })
    resp <- get_responses(dformula)
    ro <- attr(dformula, "rank_order")
    for (i in 2:n_time) {
      idx <- rep(idx_base, each = 2) + (i-1):i
      idx_i <- idx_base + i
      assign_deterministic(data, resp, exprs, ro, idx, idx_i, groups = group_var)
    }
  }
}

assign_initial_values <- function(data, dformula, group_var) {
  resp <- get_responses(dformula)
  init <- has_past(dformula)
  first_obs <- data[, .I[1], by = group_var]$V1
  if (any(init)) {
    init_idx <- which(init)
    for (i in init_idx) {
      data[first_obs, (dformula[[i]]$response) := dformula[[i]]$specials$past]
    }
  }
  #if (any(!init)) {
  #  det_from_stoch <- sapply(dformula, function(y){
  #    isTRUE(attr(y, "stoch_origin"))
  #  })
  #  stoch_pre <- !init & det_from_stoch
  #  det_pre <- !init & !det_from_stoch
  #  if (any(stoch_pre)) {
  #    formulas <- get_formulas(dformula[stoch_pre])
  #    for (i in seq_along(formulas)) {
  #      e <- formula2expression(formula_rhs(formulas[[i]]))
  #      data[first_obs, (resp[stoch_pre[i]]) := eval(e), by = group_var]
  #    }
  #  }
  #  if (any(det_pre)) {
  #    formulas <- get_formulas(dformula[det_pre])
  #    for (i in seq_along(formulas)) {
  #      e <- formula2expression(formula_rhs(formulas[[i]]))
  #      data[first_obs, (resp[det_pre[i]]) := eval(e), by = group_var]
  #    }
  #  }
  #}
}

assign_deterministic <- function(data, resp, exprs, ro, idx, idx_i, groups) {
  for (j in seq_along(ro)) {
    data[idx_i, (resp[ro[j]]) :=
           data[idx, eval(exprs[[ro[j]]])[2], by = groups]$V1]
  }
}

# TODO documentation
evaluate_deterministic <- function(e, dformula, group_var, time_var) {

  n_time <- length(unique(e$data[[time_var]]))
  n_id <- length(unique(e$data[[group_var]]))
  idx_base <- seq(0, n_time * (n_id - 1), by = n_time)
  assign_initial_values(e, dformula, idx_base + 1, n_time, n_id)
  if (length(dformula) > 0 && n_time > 1) {
    for (i in 2:n_time) {
      idx <- rep(idx_base, each = 2) + (i-1):i
      idx_i <- idx_base + i
      assign_deterministic(e, dformula, idx, idx_i, n_id)
    }
  }
}

assign_initial_values <- function(e, dformula, idx, n_time, n_id) {
  resp <- get_responses(dformula)
  init <- has_past(dformula)
  if (any(init)) {
    init_idx <- which(init)
    for (i in init_idx) {
      e$data[idx, dformula[[i]]$response] <- dformula[[i]]$specials$past
    }
  }
  if (any(!init)) {
    det_from_stoch <- sapply(dformula, function(y){
      isTRUE(attr(y, "stoch_origin"))
    })
    det_pre <- !init & det_from_stoch
    stoch_pre <- !init & !det_from_stoch
    if (any(det_pre)) {
      data_eval <- e$data[rep(idx, each = 2),]
      eval_idx <- 1 + seq(1, 2 * n_id, by = 2)
      data_eval[eval_idx - 1,] <- NA
      e$data[idx, resp[det_pre]] <-
        full_model.matrix_pseudo(get_formulas(dformula[det_pre]),
                                 data_eval)[eval_idx, ]
    }
    if (any(stoch_pre)) {
      e$data[idx, resp[stoch_pre]] <-
        full_model.matrix_pseudo(get_formulas(dformula[stoch_pre]),
                                 e$data[idx, ])
    }
  }
}

assign_deterministic <- function(e, dformula, idx, idx_i, k) {
  rl <- attr(dformula, "rank_list")
  n_rank <- length(rl)
  data_eval <- e$data[rep(1, k * 3),]
  data_eval[,] <- NA
  eval_idx <- rep(1:2, k)  + rep(seq(1, 3 * k, by = 3), each = 2)
  model_idx <- seq(3, 3 * k, by = 3)
  for (j in seq_len(n_rank)) {
    det_idx <- rl[[j]]
    data_eval[eval_idx,] <- e$data[idx,]
    e$data[idx_i, get_responses(dformula[det_idx])] <-
      full_model.matrix_pseudo(get_formulas(dformula[det_idx]),
                               data_eval)[model_idx,]
  }
}

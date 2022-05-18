# TODO documentation
evaluate_deterministic <- function(e, dformula, group_var, time_var) {
  resp <- get_responses(dformula)
  rank <- get_ranks(dformula)
  n_time <- length(unique(e$data[[time_var]]))
  n_id <- length(unique(e$data[[group_var]]))
  id_offset <- seq(0, n_time * (n_id - 1), by = n_time)
  init <- has_past(dformula)
  if (any(init)) {
    idx <- which(init)
    for (i in idx) {
      e$data[1 + id_offset, dformula[[i]]$response] <-
        dformula[[i]]$specials$past
    }
  }
  if (any(!init)) {
    det_from_stoch <- sapply(dformula, function(y){
      isTRUE(attr(y, "stoch_origin"))
    })
    det_pre <- !init & det_from_stoch
    stoch_pre <- !init & !det_from_stoch
    if (any(det_pre)) {
      data_eval <- e$data[rep(1 + id_offset, each = 2),]
      eval_idx <- seq(1, 2 * n_id, by = 2)
      data_eval[eval_idx,] <- NA
      e$data[1 + id_offset, resp[det_pre]] <-
        full_model.matrix_pseudo(get_formulas(dformula[det_pre]),
                                 data_eval)[eval_idx + 1, ]
    }
    if (any(stoch_pre)) {
      e$data[1 + id_offset, resp[stoch_pre]] <-
        full_model.matrix_pseudo(get_formulas(dformula[stoch_pre]),
                                 e$data[1 + id_offset, ])
    }
  }
  if (length(dformula) > 0 && n_time > 1) {
    id_offset_vec <- rep(id_offset, each = 2)
    model_idx <- seq(3, 3 * n_id, by = 3)
    # TODO check if separate evaluation data frame is needed
    data_eval <- e$data[rep(1, n_id * 3),]
    data_eval[,] <- NA
    eval_idx <- rep(1:2, n_id)  + rep(seq(1, 3 * n_id, by = 3), each = 2)
    u_rank <- sort(unique(rank))
    for (i in 2:n_time) {
      past_idx <- (i - 1):i + id_offset_vec
      data_idx <- i + id_offset
      for (j in u_rank) {
        data_eval[eval_idx,] <- e$data[past_idx,]
        idx <- which(rank == j)
        e$data[data_idx, get_responses(dformula[idx])] <-
          full_model.matrix_pseudo(get_formulas(dformula[idx]),
                                   data_eval)[model_idx,]
      }
    }
  }
}

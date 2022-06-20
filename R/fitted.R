#' Extract fitted values of dynamitefit
#'
#' Note that these are conditional on the observed data i.e., we don't
#' simulate new lagged values for covariates, so we underestimate the
#' uncertainty. It is typically better to use predict with type = "mean".
#' These fitted value are mostly useful only for studying one-step ahead
#' estimates.
#' @export
#' @inheritParams predict.dynamitefit
#' @srrstats {RE4.9} *Modelled values of response variables.*
fitted.dynamitefit <- function(object, newdata = NULL,
                               n_draws = NULL,  ...) {
  if (is.null(n_draws)) {
    n_draws <- ndraws(object)
  }
  n_draws <- try(as.integer(n_draws), silent = TRUE)
  if ("try-error" %in% class(n_draws)) {
    stop_("Unable to coerce {.var n_draws} to {.cls integer}.")
  }
  fixed <- as.integer(attr(object$dformulas$all, "max_lag"))
  if (is.null(newdata)) {
    newdata <- data.table::copy(object$data)
  } else {
    if (data.table::is.data.table(newdata)) {
      newdata <- data.table::copy(newdata)
    } else {
      newdata <- data.table::as.data.table(newdata)
    }
  }
  group_var <- object$group_var
  time_var <- object$time_var
  formulas_stoch <- get_formulas(object$dformulas$stoch)
  families_stoch <- get_families(object$dformulas$stoch)
  categories <- lapply(attr(object$stan$responses, "resp_class"),
                       "attr", "levels")
  resp_stoch <- get_responses(object$dformulas$stoch)
  check_newdata(newdata, object$data, type = "response",
                families_stoch, resp_stoch, categories,
                group_var, time_var)
  group <- NULL
  n_id <- 1L
  if (!is.null(group_var)) {
    group <- unique(newdata[[group_var]])
    n_id <- length(group)
  }
  time <- unique(newdata[[time_var]])
  n_time <- length(time)
  n_new <- nrow(newdata)
  n_time <- length(time)
  model_matrix <- full_model.matrix_fast(formulas_stoch, newdata,
                                         object$stan$u_names)
  model_matrix <- model_matrix[rep(seq_len(n_new), n_draws), ]
  newdata <- data.table::as.data.table(newdata)
  newdata <- newdata[rep(seq_len(n_new), each = n_draws), ]
  newdata[, ("draw") := rep(1:n_draws, n_new)]
  eval_envs <- prepare_eval_envs(object, newdata, eval_type = "fitted",
                                 predict_type = "response",
                                 resp_stoch, n_id, n_draws, group_var)
  specials <- evaluate_specials(object$dformulas$stoch, newdata)
  idx <- as.integer(newdata[ ,.I[newdata[[time_var]] == fixed]])
  for (i in seq.int(fixed + 1L, n_time)) {
    idx <- idx + n_draws
    model_matrix_sub <- model_matrix[idx, ]
    for (j in seq_along(resp_stoch)) {
      e <- eval_envs[[j]]
      e$idx <- idx
      e$time <- i - fixed
      e$model_matrix <- model_matrix_sub
      e$offset <- specials[[j]]$offset[idx]
      e$trials <- specials[[j]]$trials[idx]
      e$a_time <- ifelse_(NCOL(e$alpha) == 1, 1, i - fixed)
      eval(e$call, envir = e)
    }
  }
  data.table::setkeyv(newdata, cols = c("draw", group_var, time_var))
  data.table::setDF(newdata)
}

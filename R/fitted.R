#' Extract fitted values of dynamitefit
#'
#' TODO Note that these are conditional on the observed data i.e., we don't
#' simulate new lagged values for covariates, so we underestimate the uncertainty.
#' Better to use predict with type = "mean" TODO do it
#' @export
#' @param object An object of class \code{dynamitefit}.
#' @param newdata TODO
#' @param n_draws TODO
#' @param n_fixed TODO
#' @param ... Ignored.
fitted.dynamitefit <- function(object, newdata = NULL,
                               n_draws = NULL, n_fixed = NULL, ...) {
  if (is.null(n_draws)) {
    n_draws <- ndraws(object)
  }
  fixed <- as.integer(attr(object$dformulas$all, "max_lag"))
  if (is.null(n_fixed)) {
    n_fixed <- fixed
  } else if (n_fixed < fixed) {
    stop_("The model implies at least ", fixed, " fixed time points, ",
          "but only ", n_fixed, " were specified")
  }
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
  n_id <- length(group)
  model_matrix <- full_model.matrix_fast(formulas_stoch, newdata,
                                         object$stan$u_names)
  model_matrix <- model_matrix[rep(seq_len(n_new), n_draws), ]
  newdata <- data.table::as.data.table(newdata)
  newdata <- newdata[rep(seq_len(n_new), n_draws), ]
  newdata[, ("draw") := rep(1:n_draws, each = n_new)]
  data.table::setkeyv(newdata, c("draw", group_var, time_var))
  idx <- seq.int(fixed, newdata[,.N], by = n_time)
  eval_envs <- prepare_eval_envs(object, newdata, eval_type = "fitted",
                                 predict_type = "response",
                                 resp_stoch, n_id, n_draws)
  for (i in seq.int(fixed + 1L, n_time)) {
    idx <- idx + 1L
    model_matrix_sub <- model_matrix[idx, ]
    for (j in seq_along(resp_stoch)) {
      e <- eval_envs[[j]]
      e$idx <- idx
      e$time <- i
      e$model_matrix <- model_matrix_sub
      e$a_time <- ifelse_(NCOL(e$alpha) == 1, 1, i)
      eval(e$call, envir = e)
    }
  }
  newdata
}

#' Extract fitted values of dynamitefit
#'
#' TODO Note that these are conditional on the observed data i.e., we don't
#' simulate new lagged values for covariates, so we underestimate the uncertainty.
#' Better to use predict with type = "mean" TODO do it
#' @export
#' @param object An object of class \code{dynamitefit}.
#' @param newdata TODO
#' @param n_draws TODO
#' @param ... Ignored.
fitted.dynamitefit <- function(object, newdata = NULL, n_draws = NULL, ...) {
  if (is.null(n_draws)) {
    n_draws <- ndraws(object)
  }
  fixed <- attr(object$dformulas$lag, "max_lag")
  if (is.null(n_fixed)) {
    n_fixed <- fixed
  } else if (n_fixed < fixed) {
    stop_("The model implies at least ", fixed, " fixed time points, ",
          "but only ", n_fixed, " were specified")
  }
  if (is.null(newdata)) {
    newdata <- object$data
  } else {
    data.table::setDT(newdata)
  }
  group_var <- object$group_var
  time_var <- object$time_var
  formulas_stoch <- get_formulas(object$dformulas$stoch)
  families_stoch <- get_families(object$dformulas$stoch)
  resp_stoch <- get_responses(object$dformulas$stoch)
  check_newdata(newdata, object$data, type = "response",
                families_stoch, resp_stoch,
                group_var, time_var)

  group <- unique(newdata[[group_var]])
  time <- unique(newdata[[time_var]])
  n_time <- length(time)
  n_id <- length(group)
  n_new <- nrow(newdata)
  n_time <- length(time)
  n_id <- length(group)
  k <- n_id * n_draws
  data.table::setDT(newdata)
  newdata <- newdata[rep(seq_len(n_new), n_draws), ]
  newdata[, ("draw") := rep(1:n_draws, each = n_new)]
  data.table::setkeyv(newdata, c("draw", group_var, time_var))
  n <- newdata[,.N]
  idx <- seq.int(1L, n, by = n_time)
  idx_par <- rep(1L:n_draws, each = n_id)
  model_matrix <- full_model.matrix_fast(formulas_stoch, newdata,
                                         object$stan$u_names)
  eval_envs <- prepare_eval_envs(object, type = "fitted",
                                 newdata, resp_stoch, model_matrix)
  for (i in (fixed + 1):n_time) {
    idx <- idx + 1L
    for (j in seq_along(resp_stoch)) {
      e <- eval_envs[[j]]
      e$idx <- idx
      e$time <- i - 1
      e$a_time <- ifelse_(ncol(e$alpha) == 1, 1, i - 1)
      if (any(idx_na)) {
        eval(e$call, envir = e)
      }
    }
  }
  newdata
}

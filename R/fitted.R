#' Extract Fitted Values of a Dynamite Model
#'
#' Note that these are conditional on the observed data i.e., we don't
#' simulate new lagged values for covariates, so we underestimate the
#' uncertainty. It is typically better to use predict with `type = "mean"`.
#' These fitted value are mostly useful only for studying one-step ahead
#' estimates.
#'
#' @export
#' @inheritParams predict.dynamitefit
#' @examples
#' fitted(gaussian_example_fit, n_draws = 2L)
#'
#' @srrstats {RE4.9} Returns the fitted values.
fitted.dynamitefit <- function(object, n_draws = NULL, ...) {
  n_draws <- check_ndraws(n_draws, ndraws(object))
  fixed <- as.integer(attr(object$dformulas$all, "max_lag"))
  group_var <- object$group_var
  time_var <- object$time_var
  formulas_stoch <- get_formulas(object$dformulas$stoch)
  families_stoch <- get_families(object$dformulas$stoch)
  categories <- lapply(
    attr(object$stan$responses, "resp_class"),
    "attr",
    "levels"
  )
  resp_stoch <- get_responses(object$dformulas$stoch)
  newdata <- data.table::setDF(data.table::copy(object$data))
  newdata <- parse_newdata(
    newdata,
    object$data,
    type = "response",
    families_stoch,
    resp_stoch,
    categories,
    new_levels = "none",
    group_var,
    time_var
  )
  groups <- !is.null(group_var)
  group <- onlyif(groups, unique(newdata[[group_var]]))
  n_id <- ifelse_(groups, length(group), 1L)
  time <- unique(newdata[[time_var]])
  n_time <- length(time)
  n_new <- nrow(newdata)
  n_time <- length(time)
  model_matrix <- full_model.matrix_fast(
    formulas_stoch,
    newdata,
    object$stan$u_names
  )
  model_matrix <- model_matrix[rep(seq_len(n_new), n_draws), ]
  newdata <- data.table::as.data.table(newdata)
  newdata <- newdata[rep(seq_len(n_new), each = n_draws), ]
  newdata[, ("draw") := rep(seq_len(n_draws), n_new)]
  eval_envs <- prepare_eval_envs(
    object,
    newdata,
    eval_type = "fitted",
    predict_type = "response",
    resp_stoch,
    n_id,
    n_draws,
    new_levels = "none",
    group_var
  )
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
      e$a_time <- ifelse_(NCOL(e$alpha) == 1L, 1L, i - fixed)
      eval(e$call, envir = e)
    }
  }
  # remove extra columns
  for (i in seq_along(resp_stoch)) {
    resp <- resp_stoch[i]
    store <- glue::glue("{resp}_store")
    newdata[, c(store) := NULL]
  }
  lhs_det <- get_responses(object$dformulas$lag_det)
  lhs_stoch <- get_responses(object$dformulas$lag_stoch)
  newdata[, c(lhs_det, lhs_stoch) := NULL]
  data.table::setkeyv(newdata, cols = c("draw", group_var, time_var))
  data.table::setDF(newdata)
  newdata
}

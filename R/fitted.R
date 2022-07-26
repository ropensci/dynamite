#' Extract Fitted Values of a Dynamite Model
#'
#' Note that these are conditional on the observed data in `newdata`,
#' i.e. these are one-step estimates E(y_t|y_t-1,...,y_1). Often the
#' [predict.dynamitefit] is what you want.
#'
#' @export
#' @inheritParams predict.dynamitefit
#' @param newdata \[`data.frame`]\cr Data used in predictions.
#'   If `NULL` (default), the data used in model estimation is used for
#'   predictions as well,
#'   There should be no new time points that were not present in the data that
#'   were used to fit the model, and now new group levels can be included. Note
#'   that as `newdata` is expanded with predictions, it can be beneficial
#'   in terms of memory usage to first remove redundant columns from `newdata`.
#' @seealso
#' @examples
#' fitted(gaussian_example_fit, n_draws = 2L)
#' \dontrun{
#' set.seed(1)
#' fit <- dynamite(obs(LakeHuron ~ 1, "gaussian") + lags(),
#'   data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
#'   "id", "time", chains = 1, refresh = 0)
#'
#' # One-step ahead samples (fitted values) from the posterior
#' # (first time point is fixed due to lag in the model):
#' library(ggplot2)
#' fitted(fit) |>
#'   dplyr::filter(time > 2) |>
#'   ggplot(aes(time, LakeHuron_fitted, group = .draw)) +
#'   geom_line(alpha = 0.5) +
#'   # observed values
#'   geom_line(aes(y = LakeHuron), colour = "tomato") +
#'   theme_bw()
#'
#' # Posterior predictive distribution given the first time point:
#' predict(fit, type = "mean") |>
#'   dplyr::filter(time > 2) |>
#'   ggplot(aes(time, LakeHuron_mean, group = .draw)) +
#'   geom_line(alpha = 0.5) +
#'   # observed values
#'   geom_line(aes(y = LakeHuron), colour = "tomato") +
#'   theme_bw()
#' }
#' @srrstats {RE4.9} Returns the fitted values.
fitted.dynamitefit <- function(object, newdata = NULL, n_draws = NULL, ...) {
  n_draws <- check_ndraws(n_draws, ndraws(object))
  newdata_null <- is.null(newdata)
  if (newdata_null) {
    newdata <- data.table::setDF(data.table::copy(object$data))
  } else if (data.table::is.data.table(newdata) || is.data.frame(newdata)) {
    newdata <- data.table::setDF(data.table::copy(newdata))
  } else if (!is.data.frame(newdata)) {
    stop_("Argument {.arg newdata} must be a {.cls data.frame} object.")
  }

  fixed <- as.integer(attr(object$dformulas$all, "max_lag"))
  group_var <- object$group_var
  time_var <- object$time_var
  dd <- object$dformulas$det
  ds <- object$dformulas$stoch
  dlp <- object$dformulas$lag_pred
  dld <- object$dformulas$lag_det
  dls <- object$dformulas$lag_stoch
  formulas_stoch <- get_formulas(ds)
  families_stoch <- get_families(ds)
  categories <- lapply(
    attr(object$stan$responses, "resp_class"),
    "attr", "levels"
  )
  resp_stoch <- get_responses(ds)
  resp_det <- get_responses(dd)
  lhs_det <- get_responses(dld)
  rhs_det <- get_predictors(dld)
  lhs_stoch <- get_responses(dls)
  rhs_stoch <- get_predictors(dls)
  predictors <- setdiff(colnames(newdata),
    c(resp_stoch, lhs_stoch, group_var, time_var))
  newdata <- parse_newdata(
    newdata,
    object$data,
    type = "fitted",
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
  cl <- get_quoted(object$dformulas$det)
  n_time <- length(time)
  n_new <- nrow(newdata)
  n_det <- length(resp_det)
  n_lag_det <- length(lhs_det)
  n_lag_stoch <- length(lhs_stoch)
  ro_det <- onlyif(
    n_lag_det > 0L,
    attr(object$dformulas$lag_det, "rank_order")
  )
  ro_stoch <- seq_len(n_lag_stoch)
  initialize_deterministic(newdata, dd, dlp, dld, dls)
  idx <- seq.int(1L, n_new, by = n_time) - 1L
  assign_initial_values(newdata, dd, dlp, dld, dls, idx, fixed, group_var)
  newdata <- newdata[rep(seq_len(n_new), each = n_draws), ]
  newdata[, (".draw") := rep(seq.int(1L, n_draws), n_new)]
  eval_envs <- prepare_eval_envs(
    object,
    newdata,
    eval_type = "fitted",
    predict_type = "mean",
    resp_stoch,
    n_id,
    n_draws,
    new_levels = "none",
    group_var
  )
  specials <- evaluate_specials(object$dformulas$stoch, newdata)
  idx <- as.integer(newdata[ ,.I[newdata[[time_var]] == ..time[1L]]]) +
    (fixed - 1L) * n_draws
  skip <- TRUE
  original_times <- unique(object$data[[time_var]])
  for (i in seq.int(fixed + 1L, n_time)) {
    time_i <- which(original_times == time[i]) - fixed
    idx <- idx + n_draws
    assign_lags(newdata, ro_det, idx, lhs_det, rhs_det, skip, n_draws)
    assign_lags(newdata, ro_stoch, idx, lhs_stoch, rhs_stoch, skip, n_draws)
    model_matrix <- full_model.matrix_predict(
      formulas_stoch,
      newdata,
      idx,
      object$stan$u_names
    )
    for (j in seq_along(resp_stoch)) {
      e <- eval_envs[[j]]
      e$idx <- idx
      e$time <- time_i
      e$model_matrix <- model_matrix
      e$offset <- specials[[j]]$offset[idx]
      e$trials <- specials[[j]]$trials[idx]
      e$a_time <- ifelse_(identical(NCOL(e$alpha), 1L), 1L, time_i)
      eval(e$call, envir = e)
    }
    assign_deterministic(newdata, cl, idx)
    skip <- FALSE
  }
  lhs_det <- get_responses(object$dformulas$lag_det)
  lhs_stoch <- get_responses(object$dformulas$lag_stoch)
  newdata[, c(lhs_det, lhs_stoch) := NULL]
  data.table::setkeyv(newdata, cols = c(".draw", group_var, time_var))
  data.table::setDF(newdata)
  newdata
}

#' Predict Method for a Dynamite Model
#'
#' Obtain counterfactual predictions for a `dynamitefit` object. Note that
#' forecasting (i.e., predictions for time indices beyond the last time index
#' in the original data) is not supported by the \pkg{dynamite} package.
#'
#' @export
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param newdata \[`data.frame`]\cr Data used in predictions. Predictions are
#'   computed for missing (`NA`) values in the response variable columns, and
#'   non-missing values are assumed fixed.
#'   If `NULL` (default), the data used in model estimation is used for
#'   predictions as well, after all values in the response variable columns
#'   after the first `fixed` time points are converted to `NA` values.
#'   Missing values in predictor columns can be imputed (argument `impute`).
#'   There should be no new time points that were not present in the data that
#'   were used to fit the model. New group levels can be included, but if the
#'   model contains random intercepts, an options for the random
#'   effects for the new levels must be chosen (argument `new_levels`). Note
#'   that as `newdata` is expanded with predictions, it can be beneficial
#'   in terms of memory usage to first remove redundant columns from `newdata`.
#' @param type \[`character(1)`]\cr Type of prediction,
#'   `"response"` (default), `"mean"`, or `"link"`.
#' @param funs \[`list()`]\cr A list containing a summary function for each
#'   channel of the model that should be applied to the predicted responses
#'   for each combination of the posterior draws and the time points.
#'   If empty, the full individual level predictions are returned. Otherwise,
#'   this argument will be recycled if fewer functions are provided than
#'   the number of channels.
#' @param impute \[`character(1)`]\cr Which imputation scheme to use for
#'   missing predictor values. Currently supported options are
#'   no imputation: `"none"` (default), and
#'   last observation carried forward: `"locf"`.
#' @param new_levels \[`character(1)`]\cr
#'   Defines if and how to sample the random intercepts for observations whose
#'   group level was not present in the original data. The options are
#'     * `"none"` (the default) which will signal an error if new levels
#'       are encountered.
#'     * `"bootstrap"` which will randomly draw from the posterior samples of
#'       the random intercepts across all original levels.
#'     * `"gaussian"` which will randomly draw from a gaussian
#'       distribution using the posterior samples of the random intercept
#'       standard deviation.
#'     * `"original"` which will randomly match each new level to one of
#'       the original levels. The posterior samples of the random intercept of
#'       the matched levels will then be used for the new levels.
#'   This argument is ignored if model does not contain random intercepts.
#' @param global_fixed \[`logical(1)`]\cr If `FALSE` (the default),
#'   the first non-fixed time point is counted from the the first non-NA
#'   observation for each group member separately. Otherwise, the first
#'   non-fixed time point is counted from the first time point globally.
#'   If there are no groups, then the options are equivalent.
#' @param n_draws \[`integer(1)`]\cr Number of posterior samples to use,
#'   default is `NULL` which uses all samples.
#' @param ... Ignored.
#' @return A `data.frame` containing the predicted values.
#' @srrstats {G2.3a} Uses match.arg.
#' @srrstats {RE2.2} Missing responses can be predicted.
#' @srrstats {G2.13, G2.14, G2.14a, G2.14b, G2.14c, G2.15}
#'   Missing data is supported in `predict`, with options to
#'   impute missing predictors using the `impute` parameter.
#' @srrstats {BS3.0} NAs are supported, but non-finite values are not.
#' @srrstats {RE4.16} New group levels are supported via
#'   `new_levels` parameter.
#' @examples
#' predict(gaussian_example_fit, type = "response", n_draws = 2L)
#'
predict.dynamitefit <- function(object, newdata = NULL,
                                type = c("response", "mean", "link"),
                                funs = list(),
                                impute = c("none", "locf"),
                                new_levels = c(
                                  "none", "bootstrap", "gaussian", "original"
                                ),
                                global_fixed = FALSE,
                                n_draws = NULL, ...) {
  stopifnot_(
    !is.null(object$stanfit),
    "No Stan model fit is available."
  )
  type <- onlyif(is.character(type), tolower(type))
  type <- try(match.arg(type, c("response", "mean", "link")), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either \"response\", \"mean\", or \"link\"."
  )
  stopifnot_(
    identical(length(funs), 0L) || all(vapply(funs, is.function, logical(1L))),
    "Argument {.arg funs} must be a list of functions."
  )
  impute <- onlyif(is.character(impute), tolower(impute))
  impute <- try(match.arg(impute, c("none", "locf")), silent = TRUE)
  stopifnot_(
    !inherits(impute, "try-error"),
    "Argument {.arg type} must be either \"none\", or \"locf\"."
  )
  new_levels <- onlyif(is.character(new_levels), tolower(new_levels))
  new_levels <- try(
    match.arg(
      new_levels,
      c("none", "bootstrap", "gaussian", "original")
    ),
    silent = TRUE
  )
  stopifnot_(
    !inherits(new_levels, "try-error"),
    "Argument {.arg type} must be either \"none\", \"bootstrap\",
     \"gaussian\" or \"original\"."
  )
  stopifnot_(
    checkmate::test_flag(x = global_fixed),
    "Argument {.arg global_fixed} must be a single {.cls logical} value."
  )
  initialize_predict(
    object,
    newdata,
    type,
    eval_type = "predicted",
    funs,
    impute,
    new_levels,
    global_fixed,
    n_draws
  )
}

#' Internal Function for Both Predict and Fitted Methods
#'
#' @inheritParams predict.dynamitefit
#' @param eval_type \[`character(1)`]\cr Either `"predicted"` or `"fitted"`.
#' @noRd
initialize_predict <- function(object, newdata, type, eval_type, funs,
                               impute, new_levels, global_fixed, n_draws) {
  n_draws <- check_ndraws(n_draws, ndraws(object))
  newdata_null <- is.null(newdata)
  newdata <- check_newdata(object, newdata)
  fixed <- as.integer(attr(object$dformulas$all, "max_lag"))
  group_var <- object$group_var
  time_var <- object$time_var
  resp_stoch <- get_responses(object$dformulas$stoch)
  new_levels <- ifelse_(
    identical(length(attr(object$dformulas$stoch, "random")$responses), 0L),
    "ignore",
    new_levels
  )
  newdata <- parse_newdata(
    dformula = object$dformulas$all,
    newdata = newdata,
    data = object$data,
    type = type,
    eval_type = eval_type,
    families_stoch = get_families(object$dformulas$stoch),
    resp_stoch = resp_stoch,
    categories = lapply(
      attr(object$stan$responses, "resp_class"),
      "attr", "levels"
    ),
    clear_names = c(
      get_responses(object$dformulas$det),
      get_responses(object$dformulas$lag_det),
      get_responses(object$dformulas$lag_stoch)
    ),
    new_levels = new_levels,
    group_var = group_var,
    time_var = time_var
  )
  impute_newdata(
    newdata = newdata,
    impute = impute,
    predictors = setdiff(
      colnames(newdata),
      c(
        resp_stoch,
        get_responses(object$dformulas$lag_stoch),
        group_var,
        time_var
      )
    ),
    group_var = group_var
  )
  clear_nonfixed(
    newdata = newdata,
    newdata_null = newdata_null,
    resp_stoch = resp_stoch,
    eval_type = eval_type,
    group_var = group_var,
    time_var = time_var,
    fixed = fixed,
    global_fixed = global_fixed
  )
  initialize_deterministic(
    data = newdata,
    dd = object$dformulas$det,
    dlp = object$dformulas$lag_pred,
    dld = object$dformulas$lag_det,
    dls = object$dformulas$lag_stoch
  )
  assign_initial_values(
    data = newdata,
    idx = seq.int(1L, nrow(newdata), by = n_unique(newdata[[time_var]])) - 1L,
    dd = object$dformulas$det,
    dlp = object$dformulas$lag_pred,
    dld = object$dformulas$lag_det,
    dls = object$dformulas$lag_stoch,
    fixed = fixed,
    group_var = group_var
  )
  newdata_names <- names(newdata)
  resp_draw <- c(
    grep("_mean", newdata_names, value = TRUE),
    grep("_link", newdata_names, value = TRUE),
    grep("_fitted", newdata_names, value = TRUE),
    setdiff(
      c(
        resp_stoch,
        get_responses(object$dformulas$det),
        get_responses(object$dformulas$lag_det),
        get_responses(object$dformulas$lag_stoch)
      ),
      get_responses(object$dformulas$lag_pred)
    )
  )
  draw_dep <- newdata[, .SD, .SDcols = resp_draw]
  draw_indep <- newdata[, .SD, .SDcols = setdiff(names(newdata), resp_draw)]
  if (length(funs) > 0L) {
    predict_summary(
      object = object,
      storage = draw_dep,
      observed = draw_indep,
      funs = funs,
      new_levels = new_levels,
      n_draws = n_draws,
      fixed = fixed,
      group_var = group_var,
      time_var = time_var
    )
  } else {
    predict_full(
      object = object,
      simulated = draw_dep,
      observed = draw_indep,
      type = type,
      eval_type = eval_type,
      new_levels = new_levels,
      n_draws = n_draws,
      fixed = fixed,
      group_var = group_var,
      time_var = time_var
    )
  }
}

#' Obtain Individual Level Predictions Or Fitted Values
#'
#' @inheritParams initialize_predict
#' @noRd
predict_full <- function(object, simulated, observed,
                         type, eval_type, new_levels, n_draws, fixed,
                         group_var, time_var) {
  formulas_stoch <- get_formulas(object$dformulas$stoch)
  resp_stoch <- get_responses(object$dformulas$stoch)
  resp_det <- get_responses(object$dformulas$det)
  lhs_ld <- get_responses(object$dformulas$lag_det)
  rhs_ld <- get_predictors(object$dformulas$lag_det)
  lhs_ls <- get_responses(object$dformulas$lag_stoch)
  rhs_ls <- get_predictors(object$dformulas$lag_stoch)
  ro_ld <- onlyif(
    length(lhs_ld) > 0L,
    attr(object$dformulas$lag_det, "rank_order")
  )
  ro_ls <- seq_along(lhs_ls)
  n_new <- nrow(simulated)
  simulated <- simulated[rep(seq_len(n_new), each = n_draws), ]
  simulated[, (".draw") := rep(seq.int(1L, n_draws), n_new)]
  eval_envs <- prepare_eval_envs(
    object,
    simulated,
    observed,
    type,
    eval_type,
    resp_stoch,
    n_draws,
    new_levels,
    group_var
  )
  cl <- get_quoted(object$dformulas$det)
  specials <- evaluate_specials(object$dformulas$stoch, observed)
  time <- observed[[time_var]]
  u_time <- unique(time)
  n_time <- length(u_time)
  time_offset <- which(unique(object$data[[time_var]]) == u_time[1L]) - 1L
  draw_time <- rep(time, each = n_draws)
  idx <- which(draw_time == u_time[1L]) + (fixed - 1L) * n_draws
  idx_obs <- rep(
    which(time == u_time[1L]) + (fixed - 1L),
    each = n_draws
  )
  skip <- TRUE
  for (i in seq.int(fixed + 1L, n_time)) {
    time_i <- time_offset + i - fixed
    idx <- idx + n_draws
    idx_obs <- idx_obs + 1L
    assign_lags(simulated, idx, ro_ld, lhs_ld, rhs_ld, skip, n_draws)
    assign_lags(simulated, idx, ro_ls, lhs_ls, rhs_ls, skip, n_draws)
    model_matrix <- full_model.matrix_predict(
      formulas_stoch,
      simulated,
      observed,
      idx,
      idx_obs,
      object$stan$u_names
    )
    for (j in seq_along(resp_stoch)) {
      e <- eval_envs[[j]]
      e$idx <- idx
      e$time <- time_i
      e$model_matrix <- model_matrix
      e$offset <- specials[[j]]$offset[idx_obs]
      e$trials <- specials[[j]]$trials[idx_obs]
      e$a_time <- ifelse_(identical(NCOL(e$alpha), 1L), 1L, time_i)
      if (identical(eval_type, "predicted")) {
        idx_na <- is.na(simulated[idx, .SD, .SDcols = resp_stoch[j]])
        e$idx_out <- which(idx_na)
        e$idx_data <- idx[e$idx_out]
        if (any(idx_na)) {
          eval(e$call, envir = e)
        }
      } else {
        eval(e$call, envir = e)
      }
    }
    assign_deterministic(simulated, idx, cl)
    skip <- FALSE
  }
  finalize_predict(type, eval_type, resp_stoch, simulated, observed)
  simulated[, c(lhs_ld, lhs_ls) := NULL]
  data.table::setDF(simulated)
  data.table::setDF(observed)
  list(simulated = simulated, observed = observed)
}

#' Obtain Summarized Predictions
#'
#' @inheritParams initialize_predict
#' @noRd
predict_summary <- function(object, storage, observed, funs, new_levels,
                            n_draws, fixed, group_var, time_var) {
  formulas_stoch <- get_formulas(object$dformulas$stoch)
  resp_stoch <- get_responses(object$dformulas$stoch)
  resp_det <- get_responses(object$dformulas$det)
  lhs_ld <- get_responses(object$dformulas$lag_det)
  rhs_ld <- get_predictors(object$dformulas$lag_det)
  lhs_ls <- get_responses(object$dformulas$lag_stoch)
  rhs_ls <- get_predictors(object$dformulas$lag_stoch)
  ro_ld <- onlyif(
    length(lhs_ld) > 0L,
    attr(object$dformulas$lag_det, "rank_order")
  )
  ro_ls <- seq_along(lhs_ls)
  n_new <- nrow(storage)
  n_group <- ifelse_(is.null(group_var), 1L, n_unique(observed[[group_var]]))
  time <- observed[[time_var]]
  u_time <- unique(time)
  n_time <- length(u_time)
  idx_obs <- rep(
    which(time == u_time[1L]) + (fixed - 1L),
    each = n_draws
  )
  n_sim <- n_group * n_draws
  idx_prev <- seq.int(1L, n_sim)
  idx <- seq.int(n_sim + 1L, 2L * n_sim)
  simulated <- storage[1L, ]
  simulated[, (names(simulated)) := .SD[NA]]
  simulated <- simulated[rep(1L, 2L * n_sim), ]
  summaries <- NULL
  assign_from_storage(
    storage,
    simulated,
    idx_obs,
    idx
  )
  eval_envs <- prepare_eval_envs(
    object,
    simulated,
    observed,
    type = "response",
    eval_type = "predicted",
    resp_stoch,
    n_draws,
    new_levels,
    group_var
  )
  cl <- get_quoted(object$dformulas$det)
  specials <- evaluate_specials(object$dformulas$stoch, observed)
  #idx_summ <-
  time_offset <- which(unique(object$data[[time_var]]) == u_time[1L]) - 1L
  skip <- TRUE
  for (i in seq.int(fixed + 1L, n_time)) {
    time_i <- time_offset + i - fixed
    idx_obs <- idx_obs + 1L
    shift_simulated_values(simulated, idx_prev, idx)
    assign_from_storage(storage, simulated, idx_obs, idx)
    assign_lags(simulated, idx, ro_ld, lhs_ld, rhs_ld, skip, n_sim)
    assign_lags(simulated, idx, ro_ls, lhs_ls, rhs_ls, skip, n_sim)
    model_matrix <- full_model.matrix_predict(
      formulas_stoch,
      simulated,
      observed,
      idx,
      idx_obs,
      object$stan$u_names
    )
    for (j in seq_along(resp_stoch)) {
      e <- eval_envs[[j]]
      e$idx <- idx
      e$time <- time_i
      e$model_matrix <- model_matrix
      e$offset <- specials[[j]]$offset[idx_obs]
      e$trials <- specials[[j]]$trials[idx_obs]
      e$a_time <- ifelse_(identical(NCOL(e$alpha), 1L), 1L, time_i)
      idx_na <- is.na(simulated[idx, .SD, .SDcols = resp_stoch[j]])
      e$idx_out <- which(idx_na)
      e$idx_data <- idx[e$idx_out]
      if (any(idx_na)) {
        eval(e$call, envir = e)
      }
    }
    assign_deterministic(simulated, idx, cl)
    assign_summaries(summaries, simulated, funs)
    skip <- FALSE
  }
  finalize_predict(type, eval_type, resp_stoch, simulated, observed)
  data.table::setDF(observed)
  list(simulated = summaries, observed = observed)
}

assign_from_storage <- function(storage, simulated, idx_obs, idx) {
  for (n in names(storage)) {
    data.table::set(
      simulated, i = idx, j = n, value = storage[[n]][idx_obs]
    )
  }
}

shift_simulated_values <- function(simulated, idx_prev, idx) {
  for (n in names(simulated)) {
    data.table::set(
      simulated, i = idx_prev, j = n, value = simulated[[n]][idx]
    )
  }
}

assign_summaries <- function(summaries, simulated, funs) {
  invisible(NULL)
}

#' Clean Up Prediction Data Tables
#'
#' @inheritParams predict_full
#' @param resp_stoch \[character()]\cr Vector of response variable names.
#' @noRd
finalize_predict <- function(type, eval_type, resp_stoch, simulated, observed) {
  for (resp in resp_stoch) {
    store <- glue::glue("{resp}_store")
    if (identical(type, "response")) {
      simulated[, glue::glue("{resp}_new") := simulated[[resp]]]
    }
    observed[, c(resp) := observed[[store]]]
    observed[, c(store) := NULL]
    simulated[, c(resp) := NULL]
  }
}

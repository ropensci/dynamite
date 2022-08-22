#' Predict Method for a Dynamite Model
#'
#' Obtain counterfactual predictions for a `dynamitefit` object. Note that
#' forecasting (i.e., predictions for time indices beyond the last time index
#' in the original data) is not supported by the \pkg{dynamite} package.
#'
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
#' @param impute \[`character(1)`]\cr Which imputation scheme to use for
#'   missing predictor values. Currently supported options are
#'   no imputation: `"none"` (default), and
#'   last observation carried forward: `"locf"`.
#' @param new_levels \[`character(1)`: \sQuote{NULL}]\cr
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
#' @export
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
predict.dynamitefit <- function(object, newdata = NULL,
                                type = c("response", "mean", "link"),
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
  predict_dynamitefit(
    object,
    newdata,
    type,
    eval_type = "predict",
    impute,
    new_levels,
    global_fixed,
    n_draws
  )
}

#' Internal Function for Both Predict and Fitted Methods
#'
#' @inheritParams predict.dynamitefit
#' @param eval_type \[`character(1)`]\cr Either `"predict"` or `"fitted"`.
#' @noRd
predict_dynamitefit <- function(object, newdata, type, eval_type,
                                impute, new_levels, global_fixed, n_draws) {
  n_draws <- check_ndraws(n_draws, ndraws(object))
  newdata_null <- is.null(newdata)
  newdata <- check_newdata(object, newdata)
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
  predictors <- setdiff(
    colnames(newdata),
    c(resp_stoch, lhs_stoch, group_var, time_var)
  )
  new_levels <- ifelse_(
    identical(length(attr(object$dformulas$stoch, "random")$responses), 0L),
    "ignore",
    new_levels
  )
  newdata <- parse_newdata(
    object$dformulas$all,
    newdata,
    object$data,
    type,
    eval_type,
    families_stoch,
    resp_stoch,
    categories,
    new_levels,
    group_var,
    time_var
  )
  impute_newdata(newdata, impute, predictors, group_var)
  time <- unique(newdata[[time_var]])
  cl <- get_quoted(object$dformulas$det)
  n_time <- length(time)
  n_new <- nrow(newdata)
  ro_det <- onlyif(
    length(lhs_det) > 0L,
    attr(object$dformulas$lag_det, "rank_order")
  )
  ro_stoch <- seq_len(length(lhs_stoch))
  clear_nonfixed(
    newdata,
    newdata_null,
    resp_stoch,
    eval_type,
    group_var,
    time_var,
    clear_names = c(resp_det, lhs_det, lhs_stoch),
    fixed,
    global_fixed
  )
  initialize_deterministic(newdata, dd, dlp, dld, dls)
  idx <- seq.int(1L, n_new, by = n_time) - 1L
  assign_initial_values(newdata, dd, dlp, dld, dls, idx, fixed, group_var)
  newdata <- newdata[rep(seq_len(n_new), each = n_draws), ]
  newdata[, (".draw") := rep(seq.int(1L, n_draws), n_new)]
  eval_envs <- prepare_eval_envs(
    object,
    newdata,
    type,
    eval_type,
    resp_stoch,
    n_draws,
    new_levels,
    group_var
  )
  specials <- evaluate_specials(object$dformulas$stoch, newdata)
  idx <- which(newdata[[time_var]] == time[1L]) + (fixed - 1L) * n_draws
  time_offset <- which(unique(object$data[[time_var]]) == time[1L]) - 1L
  skip <- TRUE
  for (i in seq.int(fixed + 1L, n_time)) {
    time_i <- time_offset + i - fixed
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
      if (identical(eval_type, "predict")) {
        idx_na <- is.na(newdata[idx, .SD, .SDcols = resp_stoch[j]])
        e$idx_pred <- which(idx_na)
        e$idx_data <- idx[e$idx_pred]
        if (any(idx_na)) {
          eval(e$call, envir = e)
        }
      } else {
        eval(e$call, envir = e)
      }
    }
    assign_deterministic(newdata, cl, idx)
    skip <- FALSE
  }
  if (identical(eval_type, "predict")) {
    for (i in seq_along(resp_stoch)) {
      resp <- resp_stoch[i]
      store <- glue::glue("{resp}_store")
      if (identical(type, "response")) {
        newdata[, glue::glue("{resp}_new") := newdata[[resp]]]
      }
      newdata[, c(resp) := newdata[[store]]]
      newdata[, c(store) := NULL]
    }
  }
  newdata[, c(lhs_det, lhs_stoch) := NULL]
  data.table::setkeyv(newdata, cols = c(".draw", group_var, time_var))
  data.table::setDF(newdata)
  newdata
}

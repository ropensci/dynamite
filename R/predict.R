#' Predict method for a Bayesian Time-Varying Coefficients Model
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
                                n_draws = NULL, ...) {

  type <- onlyif(is.character(type), tolower(type))
  type <- try(match.arg(type, c("response", "mean", "link")), silent = TRUE)
  stopifnot_(
    !"try-error" %in% class(type),
    "Argument {.arg type} must be either \"response\", \"mean\", or \"link\"."
  )
  impute <- onlyif(is.character(impute), tolower(impute))
  impute <- try(match.arg(impute, c("none", "locf")), silent = TRUE)
  stopifnot_(
    !"try-error" %in% class(impute),
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
    !"try-error" %in% class(new_levels),
    "Argument {.arg type} must be either \"none\", \"bootstrap\",
     \"gaussian\" or \"original\"."
  )
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
    type,
    families_stoch,
    resp_stoch,
    categories,
    new_levels,
    group_var,
    time_var
  )
  impute_newdata(newdata, impute, predictors, group_var)
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
  clear_nonfixed(
    newdata,
    newdata_null,
    resp_stoch,
    group_var,
    clear_names = c(resp_det, lhs_det, lhs_stoch),
    fixed,
    n_id,
    n_time
  )
  initialize_deterministic(newdata, dd, dlp, dld, dls)
  idx <- seq.int(1L, n_new, by = n_time) - 1L
  assign_initial_values(newdata, dd, dlp, dld, dls, idx, fixed, group_var)
  newdata <- newdata[rep(seq_len(n_new), each = n_draws), ]
  newdata[, (".draw") := rep(seq.int(1L, n_draws), n_new)]
  n <- newdata[, .N]
  eval_envs <- prepare_eval_envs(
    object,
    newdata,
    eval_type = "predict",
    predict_type = type,
    resp_stoch,
    n_id,
    n_draws,
    new_levels,
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
      idx_na <- is.na(newdata[idx, .SD, .SDcols = resp_stoch[j]])
      e$idx <- idx
      e$time <- time_i
      e$idx_pred <- which(idx_na)
      e$idx_data <- idx[e$idx_pred]
      e$model_matrix <- model_matrix
      e$offset <- specials[[j]]$offset[idx]
      e$trials <- specials[[j]]$trials[idx]
      e$a_time <- ifelse_(identical(NCOL(e$alpha), 1L), 1L, time_i)
      if (any(idx_na)) {
        eval(e$call, envir = e)
      }
    }
    assign_deterministic(newdata, cl, idx)
    skip <- FALSE
  }
  for (i in seq_along(resp_stoch)) {
    resp <- resp_stoch[i]
    store <- glue::glue("{resp}_store")
    if (identical(type, "response")) {
      newdata[, glue::glue("{resp}_new") := newdata[[resp]]]
    }
    newdata[, c(resp) := newdata[[store]]]
    newdata[, c(store) := NULL]
  }
  newdata[, c(lhs_det, lhs_stoch) := NULL]
  data.table::setkeyv(newdata, cols = c(".draw", group_var, time_var))
  data.table::setDF(newdata)
  newdata
}

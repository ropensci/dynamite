#' Predict method for a Bayesian Time-Varying Coefficients Model
#'
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param newdata \[`data.frame`]\cr Data used in predictions.
#'   If `NULL` (default), the data used in model estimation is used for
#'   predictions as well. Predictions are computed for missing values in the
#'   response variable columns, and non-missing values are assumed fixed. If
#'   `newdata` is `NULL`, all values in the response variable columns
#'   after the first `fixed` time points will be converted to `NA` values.
#'   Missing values in predictor columns can be imputed (argument `impute`).
#'   There should be no new time points that were not present in the data that
#'   were used to fit the model. New group levels can be included, but if the
#'   model contains random intercepts, an options for the random
#'   effects for the new levels must be chosen (argument `new_levels`).
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
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @srrstats {RE2.2} *Regression Software should provide different options for processing missing values in predictor and response data. For example, it should be possible to fit a model with no missing predictor data in order to generate values for all associated response points, even where submitted response values may be missing.*
#' @srrstats {G2.8} *Software should provide appropriate conversion or dispatch routines as part of initial pre-processing to ensure that all other sub-functions of a package receive inputs of a single defined class or type.*
#' @srrstats {G2.9} *Software should issue diagnostic messages for type conversion in which information is lost (such as conversion of variables from factor to character; standardisation of variable names; or removal of meta-data such as those associated with [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as insertion of variable or column names where none were provided).*
#' @srrstats {G2.13} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#' @srrstats {G2.14} *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
#' @srrstats {G2.14a} *error on missing data*
#' @srrstats {G2.14b} *ignore missing data with default warnings or messages issued*
#' @srrstats {G2.14c} *replace missing data with appropriately imputed values*
#' @srrstats {G2.15} *Functions should never assume non-missingness, and should never pass data with potential missing values to any base routines with default `na.rm = FALSE`-type parameters (such as [`mean()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/mean.html), [`sd()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/sd.html) or [`cor()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html)).*
#' @srrstats {BS3.0} *Explicitly document assumptions made in regard to missing values; for example that data is assumed to contain no missing (`NA`, `Inf`) values, and that such values, or entire rows including any such values, will be automatically removed from input data.*
#' @srrstats {RE4.14} *Where possible, values should also be provided for extrapolation or forecast *errors*.*
#' @srrstats {RE4.16} *Regression Software which models distinct responses for different categorical groups should include the ability to submit new groups to `predict()` methods.*
predict.dynamitefit <- function(object, newdata = NULL,
                                type = c("response", "mean", "link"),
                                impute = c("none", "locf", "linear"),
                                new_levels = c("none", "bootstrap", "gaussian",
                                               "original"),
                                n_draws = NULL, ...) {
  type <- match.arg(type)
  impute <- match.arg(impute)
  new_levels <- match.arg(new_levels)
  n_draws <- check_ndraws(n_draws, ndraws(object))
  newdata_null <- is.null(newdata)
  if (newdata_null) {
    newdata <- data.table::setDF(data.table::copy(object$data))
  } else if (data.table::is.data.table(newdata)) {
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
  categories <- lapply(attr(object$stan$responses, "resp_class"),
                       "attr", "levels")
  resp_stoch <- get_responses(ds)
  resp_det <- get_responses(dd)
  lhs_det <- get_responses(dld)
  rhs_det <- get_predictors(dld)
  lhs_stoch <- get_responses(dls)
  rhs_stoch <- get_predictors(dls)
  newdata <- parse_newdata(newdata, object$data, type, families_stoch,
                           resp_stoch, categories, new_levels,
                           group_var, time_var)
  if (!identical(impute, "none")) {
    predictors <- setdiff(names(newdata), resp_stoch)
    impute_newdata(newdata, impute, predictors, group_var)
  }
  group <- NULL
  n_id <- 1L
  if (!is.null(group_var)) {
    group <- unique(newdata[[group_var]])
    n_id <- length(group)
  }
  time <- unique(newdata[[time_var]])
  cl <- get_quoted(object$dformulas$det)
  n_time <- length(time)
  n_new <- nrow(newdata)
  n_det <- length(resp_det)
  n_lag_det <- length(lhs_det)
  n_lag_stoch <- length(lhs_stoch)
  if (n_lag_det > 0) {
    ro_det <- attr(object$dformulas$lag_det, "rank_order")
  }
  if (n_lag_stoch > 0) {
    ro_stoch <- 1:n_lag_stoch
  }
  clear_nonfixed(newdata, newdata_null, resp_stoch, group_var,
                 clear_names = c(resp_det, lhs_det, lhs_stoch),
                 fixed, n_id, n_time)
  initialize_deterministic(newdata, dd, dlp, dld, dls)
  idx <- seq.int(1L, n_new, by = n_time) - 1L
  assign_initial_values(newdata, dd, dlp, dld, dls, idx, fixed, group_var)
  newdata <- newdata[rep(seq_len(n_new), each = n_draws), ]
  newdata[, ("draw") := rep(1:n_draws, n_new)]
  n <- newdata[, .N]
  eval_envs <- prepare_eval_envs(object, newdata,
                                 eval_type = "predict", predict_type = type,
                                 resp_stoch, n_id, n_draws,
                                 new_levels, group_var)
  specials <- evaluate_specials(object$dformulas$stoch, newdata)
  idx <- as.integer(newdata[ ,.I[newdata[[time_var]] == fixed]])
  for (i in (fixed + 1L):n_time) {
    idx <- idx + n_draws
    if (n_lag_det > 0 && i > fixed + 1L) {
      assign_lags(newdata, ro_det, idx, lhs_det, rhs_det, n_draws)
    }
    if (n_lag_stoch > 0 && i > fixed + 1L) {
      assign_lags(newdata, ro_stoch, idx, lhs_stoch, rhs_stoch, n_draws)
    }
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
      e$time <- i - fixed
      e$idx_pred <- which(idx_na)
      e$idx_data <- idx[e$idx_pred]
      e$model_matrix <- model_matrix
      e$offset <- specials[[j]]$offset[idx]
      e$trials <- specials[[j]]$trials[idx]
      e$a_time <- ifelse_(NCOL(e$alpha) == 1, 1, i - fixed)
      if (any(idx_na)) {
        eval(e$call, envir = e)
      }
    }
    if (n_det > 0) {
      assign_deterministic(newdata, cl, idx)
    }
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
  if (n_lag_det > 0 || n_lag_stoch > 0) {
    newdata[, c(lhs_det, lhs_stoch) := NULL]
  }
  data.table::setkeyv(newdata, cols = c("draw", group_var, time_var))
  # for consistency with other output types
  data.table::setDF(newdata)
}

#' Predict method for a Bayesian Time-Varying Coefficients Model
#'
#' @param object Object of class `dynamitefit`.
#' @param newdata Data frame used in predictions. If `NULL` (default), the
#'   data used in model estimation is used for predictions as well.
#' @param mode TODO: Think about the names, if we predict missing observations
#'   does the term counterfactual make sense?
#' @param type Type of prediction, `"response"` (default), `"mean"`
#'   or `"link"`.
#' @param n_draws Number of posterior samples to use, default is all.
#' @param n_fixed Number of observations per individual that should be assumed
#'  fixed and will not be predicted (counting from the first non-NA observation)
#' @param ... Ignored.
#' @export
predict.dynamitefit <- function(object, newdata = NULL,
                                mode = c("counterfactual", "forecast"),
                                type = c("response", "mean", "link"),
                                n_draws = NULL, n_fixed = NULL, ...) {
  mode <- match.arg(mode)
  type <- match.arg(type)
  # TODO check args
  do.call(paste0("predict.dynamitefit_", mode),
          list(object = object, newdata = newdata,
               type = type, n_draws = n_draws, n_fixed = n_fixed))
}

predict.dynamitefit_counterfactual <- function(object, newdata, type,
                                               n_draws, n_fixed) {
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
  newdata_null <- is.null(newdata)
  if (newdata_null) {
    newdata <- data.table::copy(object$data)
  } else {
    if (data.table::is.data.table(newdata)) {
      newdata <- data.table::copy(newdata)
    } else {
      newdata <- data.table::as.data.table(newdata)
    }
  }
  # TODO impute predictor values
  group_var <- object$group_var
  time_var <- object$time_var
  formulas_stoch <- get_formulas(object$dformulas$stoch)
  families_stoch <- get_families(object$dformulas$stoch)
  categories <- lapply(attr(object$stan$responses, "resp_class"),
                       "attr", "levels")
  resp_stoch <- get_responses(object$dformulas$stoch)
  resp_det <- get_responses(object$dformulas$det)
  lhs_det <- get_responses(object$dformulas$lag_det)
  rhs_det <- get_predictors(object$dformulas$lag_det)
  lhs_stoch <- get_responses(object$dformulas$lag_stoch)
  rhs_stoch <- get_predictors(object$dformulas$lag_stoch)
  check_newdata(newdata, object$data, type, families_stoch,
                resp_stoch, categories, group_var, time_var)
  group <- unique(newdata[[group_var]])
  time <- unique(newdata[[time_var]])
  cl <- get_quoted(object$dformulas$det)
  n_time <- length(time)
  n_id <- length(group)
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
                 n_fixed, n_id, n_time)
  newdata <- newdata[rep(seq_len(n_new), n_draws), ]
  newdata[, ("draw") := rep(1:n_draws, each = n_new)]
  n <- newdata[,.N]
  idx <- seq.int(1L, n, by = n_time)
  assign_initial_values(newdata, object$dformulas$det,
                        object$dformulas$lag_det, object$dformulas$lag_stoch,
                        idx)
  eval_envs <- prepare_eval_envs(object, newdata,
                                 eval_type = "predict", predict_type = type,
                                 resp_stoch, n_id, n_draws)
  idx <- idx + fixed - 1L
  for (i in (fixed + 1):n_time) {
    idx <- idx + 1L
    if (n_lag_det > 0 && idx > fixed + 1) {
      assign_lags(newdata, ro_det, idx, lhs_det, rhs_det)
    }
    if (n_lag_stoch > 0) {
      assign_lags(newdata, ro_stoch, idx, lhs_stoch, rhs_stoch)
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
      e$time <- i - 1
      e$idx_pred <- idx[which(idx_na)]
      e$model_matrix <- model_matrix
      e$a_time <- ifelse_(NCOL(e$alpha) == 1, 1, i - 1)
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
    if (is_categorical(object$dformulas$stoch[[i]]$family)) {
      #newdata[[resp]] <- NULL
      # TODO categorical case
    }
    if (identical(type, "response")) {
      newdata[[glue::glue("{resp}_new")]] <- newdata[[resp]]
    }
    newdata[[resp]] <- newdata[[store]]
    newdata[[store]] <- NULL
  }
  newdata[,c(lhs_det, lhs_stoch) := NULL]

  # for consistency with other output types
  data.table::setDF(newdata)
}

predict.dynamitefit_forecast <- function(object, newdata, type, n_draws) {
  stop_("Forecasting is not yet supported")
}

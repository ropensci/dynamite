#' Predict Method for a Dynamite Model
#'
#' Obtain counterfactual predictions for a `dynamitefit` object.
#'
#' Note that forecasting (i.e., predictions for time indices beyond the last
#' time index in the original data) is not supported by the \pkg{dynamite}
#' package. However, such predictions can be obtained by augmenting the
#' original data with NA values before model estimation.
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
#'   model contains random effects, an option for the random
#'   effects for the new levels must be chosen (argument `new_levels`).
#'   If the grouping variable of the original data is missing, it is assumed
#'   that all observations in `newdata` belong to the first group in the
#'   original data.
#' @param type \[`character(1)`]\cr Type of prediction,
#'   `"response"` (default), `"mean"`, or `"link"`.
#' @param funs \[`list()`]\cr A named list whose names should correspond to the
#'   response variables of the model. Each element of `funs` should be a
#'   a named `list` of functions that will be applied to the
#'   corresponding predicted `type` of the channel over the individuals
#'   for each combination of the posterior draws and time points.
#'   In other words, the resulting predictions will be averages
#'   over the individuals. The functions should take the corresponding
#'   `type` variable values as their only argument.
#'   If `funs` is empty, the full individual level values are returned
#'   instead. Note that this argument can only be used
#'   if there are multiple individuals (i.e., `group` was not `NULL` in the
#'   `dynamite` call).
#' @param impute \[`character(1)`]\cr Which imputation scheme to use for
#'   missing predictor values. Currently supported options are
#'   no imputation: `"none"` (default), and
#'   last observation carried forward: `"locf"`.
#' @param new_levels \[`character(1)`]\cr
#'   Defines if and how to sample the random effects for observations whose
#'   group level was not present in the original data. The options are
#'     * `"none"` (the default) which will signal an error if new levels
#'       are encountered.
#'     * `"bootstrap"` which will randomly draw from the posterior samples of
#'       the random effects across all original levels.
#'     * `"gaussian"` which will randomly draw from a gaussian
#'       distribution using the posterior samples of the random effects
#'       standard deviation (and correlation matrix if applicable).
#'     * `"original"` which will randomly match each new level to one of
#'       the original levels. The posterior samples of the random effects of
#'       the matched levels will then be used for the new levels.
#'   This argument is ignored if model does not contain random effects.
#' @param global_fixed \[`logical(1)`]\cr If `FALSE` (the default),
#'   the first non-fixed time point is counted from the the first non-NA
#'   observation for each group member separately. Otherwise, the first
#'   non-fixed time point is counted from the first time point globally.
#'   If there are no groups, then the options are equivalent.
#' @param n_draws \[`integer(1)`]\cr Number of posterior samples to use,
#'   default is `NULL` which uses all samples.
#' @param expand \[`logical(1)`]\cr If `TRUE` (the default), the output
#'   is a single `data.frame` containing the original `newdata` and the
#'   predicted values. Otherwise, a `list` is returned with two components,
#'   `simulated` and `observed`, where the first contains only the
#'   predicted values, and the second contains the original `newdata`.
#'   Setting `expand` to `FALSE` can help conserve memory because `newdata`
#'   is not replicated `n_draws` times in the output.
#'   This argument is ignored if `funs` are provided.
#' @param df \[`logical(1)`]\cr If `TRUE` (default) the output
#'   consists of `data.frame` objects, and `data.table` objects otherwise.
#' @param ... Ignored.
#' @return A `data.frame` containing the predicted values or a `list` of two
#'   `data.frames`. See the `expand` argument for details. Note that the
#'   `.draw` column is not the same as `.draw` from `as.data.frame` and
#'   `as_draws` methods as `predict` uses permuted samples. A mapping between
#'   these variables can be done using information in
#'   `object$stanfit@sim$permutation`.
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
#' \donttest{
#' # Simulate from the prior predictive distribution
#'
#' f <- obs(y ~ lag(y) + varying(~ -1 + x), "gaussian") +
#'   splines(df = 10, noncentered = TRUE)
#'
#' # Create data with missing observations
#' # Note that due to the lagged term in the model,
#' # we need to fix the first time point
#' d <- data.frame(y = c(0, rep(NA, 49)), x = rnorm(50), time = 1:50)
#'
#' # suppress warnings due to the lack of data
#' suppressWarnings(
#'   priors <- get_priors(f, data = d, time = "time")
#' )
#'
#' # modify default priors which can produce exploding behavior when used
#' # without data
#' priors$prior <- c(
#'   "normal(0, 1)",
#'   "normal(0.6, 0.1)",
#'   "normal(-0.2, 0.5)",
#'   "normal(0.2, 0.1)",
#'   "normal(0.5, 0.1)"
#' )
#'
#' # samples from the prior conditional on the first time point and x
#' fit <- dynamite(
#'   dformula = f,
#'   data = d,
#'   time = "time",
#'   verbose = FALSE,
#'   priors = priors,
#'   chains = 1
#' )
#'
#' # simulate new data
#' pp <- predict(fit)
#'
#' ggplot2::ggplot(pp, ggplot2::aes(time, y_new, group = .draw)) +
#'   ggplot2::geom_line(alpha = 0.1) +
#'   ggplot2::theme_bw()
#' }
#'
predict.dynamitefit <- function(object, newdata = NULL,
                                type = c("response", "mean", "link"),
                                funs = list(), impute = c("none", "locf"),
                                new_levels = c(
                                  "none", "bootstrap", "gaussian", "original"
                                ),
                                global_fixed = FALSE, n_draws = NULL,
                                expand = TRUE, df = TRUE, ...) {
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
  stopifnot_(
    checkmate::test_flag(x = expand),
    "Argument {.arg expand} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_flag(x = df),
    "Argument {.arg df} must be a single {.cls logical} value."
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
    n_draws,
    expand,
    df
  )
}

#' Internal Function for Both Predict and Fitted Methods
#'
#' @inheritParams predict.dynamitefit
#' @param eval_type \[`character(1)`]\cr Either `"predicted"`, `"fitted"`, or
#'   `"loglik"`.
#' @noRd
initialize_predict <- function(object, newdata, type, eval_type, funs, impute,
                               new_levels, global_fixed, n_draws, expand, df) {
  n_draws <- check_ndraws(n_draws, ndraws(object))
  newdata_null <- is.null(newdata)
  newdata <- check_newdata(object, newdata)
  fixed <- as.integer(attr(object$dformulas$all, "max_lag"))
  group_var <- object$group_var
  time_var <- object$time_var
  resp_stoch <- get_responses(object$dformulas$stoch)
  new_levels <- ifelse_(
    length(which_random(object$dformulas$all)) == 0L,
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
    grep("_loglik", newdata_names, value = TRUE),
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
  draw_dep <- newdata[, .SD, .SDcols = c(resp_draw, group_var, time_var)]
  draw_indep <- newdata[, .SD, .SDcols = setdiff(names(newdata), resp_draw)]
  if (length(funs) > 0L) {
    funs <- parse_funs(object, type, funs)
    data.table::setcolorder(
      x = draw_dep,
      neworder = c(group_var, time_var, resp_draw)
    )
    predict_summary(
      object = object,
      storage = draw_dep,
      observed = draw_indep,
      type = type,
      funs = funs,
      new_levels = new_levels,
      n_draws = n_draws,
      fixed = fixed,
      group_var = group_var,
      time_var = time_var,
      df = df
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
      time_var = time_var,
      expand = expand,
      df = df
    )
  }
}

#' Obtain Full Individual Level Predictions Or Fitted Values
#'
#' @inheritParams initialize_predict
#' @param simulated A `data.table` containing the simulated values.
#' @param observed A `data.table` containing fixed predictors (values
#'   independent of the posterior draw).
#' @noRd
predict_full <- function(object, simulated, observed, type, eval_type,
                         new_levels, n_draws, fixed,
                         group_var, time_var, expand, df) {
  formulas_stoch <- object$dformulas$stoch
  resp <- get_responses(object$dformulas$all)
  resp_stoch <- get_responses(object$dformulas$stoch)
  families <- get_families(object$dformulas$all)
  model_topology <- attr(object$dformulas$all, "model_topology")
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
  simulated <- simulated[rep(seq_len(n_new), each = n_draws)]
  simulated[, (".draw") := rep(seq.int(1L, n_draws), n_new)]
  eval_envs <- prepare_eval_envs(
    object,
    simulated,
    observed,
    type,
    eval_type,
    n_draws,
    new_levels,
    group_var
  )
  specials <- evaluate_specials(object$dformulas$all, observed)
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
    for (j in model_topology) {
      if (is_deterministic(families[[j]])) {
        assign_deterministic_predict(
          simulated,
          idx,
          resp[j],
          formula_rhs(object$dformulas$all[[j]]$formula)
        )
      } else {
        model_matrix <- full_model.matrix_predict(
          formulas_stoch,
          simulated,
          observed,
          idx,
          idx_obs,
          object$stan$u_names
        )
        e <- eval_envs[[j]]
        e$idx <- idx
        e$time <- time_i
        e$model_matrix <- model_matrix
        e$offset <- specials[[j]]$offset[idx_obs]
        e$trials <- specials[[j]]$trials[idx_obs]
        e$y <- observed[[paste0(resp[j], "_store")]][idx_obs]
        e$y <- ifelse_(
          is_categorical(families[[j]]),
          as.integer(e$y),
          e$y
        )
        e$a_time <- ifelse_(identical(NCOL(e$alpha), 1L), 1L, time_i)
        if (identical(eval_type, "predicted")) {
          idx_na <- is.na(simulated[idx, .SD, .SDcols = resp[j]])
          e$idx_out <- which(idx_na)
          e$idx_data <- idx[e$idx_out]
          if (any(idx_na)) {
            eval(e$call, envir = e)
          }
        } else {
          eval(e$call, envir = e)
        }
      }
    }
    skip <- FALSE
  }
  finalize_predict(type, resp_stoch, simulated, observed)
  simulated[, c(lhs_ld, lhs_ls) := NULL]
  data.table::setkeyv(simulated, cols = c(".draw", group_var, time_var))
  if (expand) {
    expand_predict_output(simulated, observed, df)
  } else {
    if (df) {
      data.table::setDF(simulated)
      data.table::setDF(observed)
    }
    list(simulated = simulated, observed = observed)
  }
}

#' Obtain Summarized Predictions
#'
#' @inheritParams initialize_predict
#' @param storage A `data.table` container for storing the temporary
#'   individual level predictions (current and previous time point at each index
#'   only.)
#' @param observed A `data.table` containing fixed predictors (values
#'   independent of the posterior draw)
#' @noRd
predict_summary <- function(object, storage, observed, type, funs, new_levels,
                            n_draws, fixed, group_var, time_var, df) {
  formulas_stoch <- object$dformulas$stoch
  resp <- get_responses(object$dformulas$all)
  resp_stoch <- get_responses(object$dformulas$stoch)
  families <- get_families(object$dformulas$all)
  model_topology <- attr(object$dformulas$all, "model_topology")
  lhs_ld <- get_responses(object$dformulas$lag_det)
  rhs_ld <- get_predictors(object$dformulas$lag_det)
  lhs_ls <- get_responses(object$dformulas$lag_stoch)
  rhs_ls <- get_predictors(object$dformulas$lag_stoch)
  ro_ld <- onlyif(
    length(lhs_ld) > 0L,
    attr(object$dformulas$lag_det, "rank_order")
  )
  ro_ls <- seq_along(lhs_ls)
  n_group <- n_unique(observed[[group_var]])
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
  #simulated <- simulated[rep(1L, 2L * n_sim), , env = list(n_sim = n_sim)]
  simulated <- simulated[rep(1L, 2L * n_sim)]
  data.table::set(
    simulated,
    j = ".draw",
    value = rep(seq.int(1L, n_draws), 2L * n_group)
  )
  assign_from_storage(
    storage,
    simulated,
    idx,
    idx_obs
  )
  summaries <- storage[1L, ]
  for (f in funs) {
    #target <- f$fun(storage[[f$target]][1L])
    #summaries[, (f$name) := target, env = list(target = target)]
    summaries[, (f$name) := f$fun(storage[[f$target]][1L])]
  }
  summaries[, (names(storage)) := NULL]
  summaries[, (names(summaries)) := .SD[NA]]
  #summaries <- summaries[
  #  rep(1L, n_time * n_draws),
  #  env = list(n_time = n_time, n_draws = n_draws)
  #]
  summaries <- summaries[rep(1L, n_time * n_draws)]
  data.table::set(
    summaries,
    j = time_var,
    value =  rep(u_time, n_draws)
  )
  data.table::set(
    summaries,
    j = ".draw",
    value = rep(seq_len(n_draws), each = n_time)
  )
  idx_summ <- which(summaries[[time_var]] == u_time[1L]) + (fixed - 1L)
  eval_envs <- prepare_eval_envs(
    object,
    simulated,
    observed,
    type,
    eval_type = "predicted",
    n_draws,
    new_levels,
    group_var
  )
  specials <- evaluate_specials(object$dformulas$all, observed)
  time_offset <- which(unique(object$data[[time_var]]) == u_time[1L]) - 1L
  skip <- TRUE
  for (i in seq.int(fixed + 1L, n_time)) {
    time_i <- time_offset + i - fixed
    idx_obs <- idx_obs + 1L
    idx_summ <- idx_summ + 1L
    shift_simulated_values(simulated, idx, idx_prev)
    assign_from_storage(storage, simulated, idx, idx_obs)
    assign_lags(simulated, idx, ro_ld, lhs_ld, rhs_ld, skip, n_sim)
    assign_lags(simulated, idx, ro_ls, lhs_ls, rhs_ls, skip, n_sim)
    for (j in model_topology) {
      if (is_deterministic(families[[j]])) {
        assign_deterministic_predict(
          simulated,
          idx,
          resp[j],
          formula_rhs(object$dformulas$all[[j]]$formula)
        )
      } else {
        model_matrix <- full_model.matrix_predict(
          formulas_stoch,
          simulated,
          observed,
          idx,
          idx_obs,
          object$stan$u_names
        )
        e <- eval_envs[[j]]
        e$idx <- idx
        e$time <- time_i
        e$model_matrix <- model_matrix
        e$offset <- specials[[j]]$offset[idx_obs]
        e$trials <- specials[[j]]$trials[idx_obs]
        e$a_time <- ifelse_(identical(NCOL(e$alpha), 1L), 1L, time_i)
        idx_na <- is.na(
          #simulated[idx, .SD, .SDcols = resp_stoch[j], env = list(idx = idx)]
          simulated[idx, .SD, .SDcols = resp[j]]
        )
        e$idx_out <- which(idx_na)
        e$idx_data <- idx[e$idx_out]
        if (any(idx_na)) {
          eval(e$call, envir = e)
        }
      }
    }
    assign_summaries(summaries, simulated, funs, idx, idx_summ)
    skip <- FALSE
  }
  finalize_predict(type = NULL, resp_stoch, simulated, observed)
  if (df) {
    data.table::setDF(summaries)
    data.table::setDF(observed)
  }
  list(simulated = summaries, observed = observed)
}

#' Get And Assign Stored Values as Simulation Template
#'
#' @inheritParams predict_summary
#' @inheritParams predict_full
#' @param idx Indices of the simulated values of the current time point.
#' @param idx_obs Indices of the stored values of which to assign.
#' @noRd
assign_from_storage <- function(storage, simulated, idx, idx_obs) {
  for (n in names(storage)) {
    data.table::set(
      simulated, i = idx, j = n, value = storage[[n]][idx_obs]
    )
  }
}

#' Compute And Assing Summary Predictions
#'
#' @inheritParams predict_summary
#' @inheritParams predict_full
#' @param idx Indices of the simulated values of the current time point.
#' @param idx_summ Indices of the summarized predictions to be assigned.
#' @noRd
assign_summaries <- function(summaries, simulated, funs, idx, idx_summ) {
  for (f in funs) {
    data.table::set(
      summaries,
      i = idx_summ,
      j = f$name,
      #value = simulated[
      #  idx, lapply(.SD, fun),
      #  by = ".draw",
      #  .SDcols = target,
      #  env = list(fun = fun, target = target)
      #][[target]]
      value = simulated[
        idx,
        lapply(.SD, f$fun),
        by = ".draw",
        .SDcols = f$target
      ][[f$target]]
    )
  }
}

#' Shift Current Individuals Level Predictions Back
#'
#' @inheritParams predict_full
#' @param idx_prev Indices of the simulated values of the previous time point.
#' @param idx Indices of the simulated values of the current time point.
#' @noRd
shift_simulated_values <- function(simulated, idx, idx_prev) {
  for (n in names(simulated)) {
    data.table::set(
      simulated, i = idx_prev, j = n, value = simulated[[n]][idx]
    )
  }
}

#' Clean Up Prediction Data Tables
#'
#' @inheritParams predict_full
#' @param resp_stoch \[character()]\cr Vector of response variable names.
#' @noRd
finalize_predict <- function(type, resp_stoch, simulated, observed) {
  for (resp in resp_stoch) {
    store <- glue::glue("{resp}_store")
    if (identical(type, "response")) {
      data.table::set(
        simulated, j = glue::glue("{resp}_new"), value = simulated[[resp]]
      )
    }
    data.table::set(observed, j = resp, value = observed[[store]])
    observed[, c(store) := NULL]
    simulated[, c(resp) := NULL]
  }
}

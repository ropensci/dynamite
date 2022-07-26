#' Check Validity of `n_draws`
#'
#' @param n_draws \[`integer(1)`]\cr
#'   Number of draws to use for `fitted` or `predict`.
#' @param full_draws \[`integer(1)`]\cr
#'   Number of draws avaiable in the `dynamitefit` object.
#' @noRd
check_ndraws <- function(n_draws, full_draws) {
  if (is.null(n_draws)) {
    n_draws <- full_draws
  }
  stopifnot_(
    checkmate::test_int(
      x = n_draws,
      lower = 1L
    ),
    "Argument {.var n_draws} must be a positive {.cls integer}."
  )
  if (n_draws > full_draws) {
    warning_(c(
      "You've supplied {.arg n_draws} = {n_draws} but there are only
      {full_draws} samples available:",
      `i` = "The available samples will be used for prediction."
    ))
    n_draws <- full_draws
  }
  as.integer(n_draws)
}

#' Parse and Prepare New Data for Prediction
#'
#' @param newdata \[`data.frame`]\cr The data to be used for prediction.
#' @param data \[`data.frame`]\cr The original data used to fit the model.
#' @param type \[`character`]\cr Either `"response"`, `"mean"`, `"link"` or
#'   `"fitted"`.
#' @param families_stoch \[`list`]\cr of `dynamitefamily` object.
#' @param resp_stoch \[`character`]\cr A vector of response variables
#' @param categories \[`list`]\cr Response variable categories
#'   for categorical families.
#' @param new_levels \[`character(1)`]\cr Either `"none"`, `"bootstrap"`,
#'  `"gaussian"`, or `"original"`.
#' @param group_var \[`character(1)`,`NULL`]\cr Group variable name or `NULL`
#'   if there is only one group.
#' @param time_var \[`character(1)`]\cr Time index variable name.
#' @noRd
parse_newdata <- function(newdata, data, type, families_stoch, resp_stoch,
                          categories, new_levels, group_var, time_var) {
  if (!is.null(group_var)) {
    stopifnot_(
      group_var %in% names(newdata),
      "Can't find grouping variable {.var {group_var}} in {.var newdata}."
    )
    group <- newdata[[group_var]]
    group <- unique(group)
    extra_levels <- unique(group[!group %in% data[[group_var]]])
    stopifnot_(
      all(group %in% data[[group_var]]) || !identical(new_levels, "none"),
      c(
        "Grouping variable {.var {group_var}} contains unknown levels:",
        `x` = "Level{?s} {.val {as.character(extra_levels)}}
               {?is/are} not present in the original data.",
        `i` = "{.strong Note:} argument {.var new_levels} is {.val none}
               which disallows new levels."
      )
    )
  }
  stopifnot_(
    time_var %in% names(newdata),
    "Can't find time index variable {.var {time_var}} in {.var newdata}."
  )
  if (is.factor(newdata[[time_var]])) {
    newdata[[time_var]] <- as.integer(newdata[[time_var]])
  }
  time <- unique(newdata[[time_var]])
  original_times <- unique(data[[time_var]])
  extra_times <- time[!time %in% original_times]
  stopifnot_(
    all(time %in% original_times),
    c(
      "Time index variable {.var {time_var}} contains unknown time points:",
      `x` = "Time point{?s} {.val {as.character(extra_times)}}
             {?is/are} not present in the original data."
    )
  )
  for (resp in resp_stoch) {
    stopifnot_(
      !is.null(newdata[[resp]]),
      "Can't find response variable {.var {resp}} in {.var newdata}."
    )
  }
  newdata <- fill_time_predict(newdata, group_var, time_var,
    original_times[2] - original_times[1])
  data.table::setDT(newdata, key = c(group_var, time_var))
  # create separate column for each level of categorical variables
  for (i in seq_along(resp_stoch)) {
    resp <- resp_stoch[i]
    if (is_categorical(families_stoch[[i]])) {
      if (type %in% c("mean", "link", "fitted")) {
        resp_levels <- categories[[resp]]
        newdata[, (glue::glue("{resp}_{type}_{resp_levels}")) := NA_real_]
      }
    } else {
      if (type %in% c("mean", "link", "fitted")) {
        newdata[, (glue::glue("{resp}_{type}")) := NA_real_]
      }
    }
    if (!identical(type, "fitted")) {
      newdata[, (glue::glue("{resp}_store")) := newdata[[resp]]]
    }
  }
  newdata
}
#' Adds NA gaps to fill in missing time points in a data frame for predictions
#'
#' Note that if `impute` is `none` and model contains lagged predictors,
#' predictions will eventually fail.
#'
#' @inheritParams dynamite
#' @noRd
fill_time_predict <- function(data, group_var, time_var, time_scale) {

  time <- sort(unique(data[[time_var]]))
  if (length(time) > 1L) {
    original_order <- colnames(data)
    full_time <- seq(time[1L], time[length(time)], by = time_scale)
    groups <- !is.null(group_var)
    if (groups) {
      time_groups <- data |>
        dplyr::group_by(.data[[group_var]]) |>
        dplyr::summarise(
          has_missing = !identical(sort(.data[[time_var]]), full_time)
        )
      if (any(time_groups$has_missing)) {
        warning_(c(
          "Time index variable {.var {time_var}} of {.arg newdata} has gaps:",
          `i` = "Filling the {.arg newdata} to regular time points. This will
          lead to propagation of NA values if the model contains
          exogenous predictors and {.arg impute} is {.val none}."
        ))
        full_data_template <- expand.grid(
          time = full_time,
          group = unique(data[[group_var]])
        )
        names(full_data_template) <- c(time_var, group_var)
        data <- full_data_template |>
          dplyr::left_join(data, by = c(group_var, time_var)) |>
          dplyr::select(original_order)
      }
    } else {
      if (!identical(data[[time_var]], full_time)) {
        warning_(c(
          "Time index variable {.var {time_var}} of {.arg newdata} has gaps:",
          `i` = "Filling the {.arg newdata} to regular time points. This will
          lead to propagation of NA values if the model contains
          exogenous predictors and {.arg impute} is {.val none}."
        ))
        full_data_template <- data.frame(time = full_time)
        names(full_data_template) <- time_var
        data <- full_data_template |>
          dplyr::left_join(data, by = time_var) |>
          dplyr::select(original_order)
      }
    }
  }
  data
}
#' Impute Predictor values in New Data
#'
#' @inheritParams predict
#' @param predictors \[`character()`] A vector of predictor column names
#' @param group_var \[`character(1)`] Grouping variable name.
#' @noRd
impute_newdata <- function(newdata, impute, predictors, group_var) {
  if (identical(impute, "locf")) {
    if (is.null(group_var)) {
      newdata[, (predictors) := lapply(.SD, locf), .SDcols = predictors]
    } else {
      newdata[, (predictors) := lapply(.SD, locf), .SDcols = predictors,
              by = group_var]
    }
  }
}

#' Assign NAs to Time Indices Beyond Fixed Time Points
#'
#' @inheritParams parse_newdata
#' @param newdata_null \[logical(1)]\cr
#'   Was `newdata` `NULL` when `predict` was called?
#' @param resp_stoch \[character()]\cr Vector of response variable names.
#' @param clear_names \[character()]\cr
#'   Vector of column names to remove from `newdata`.
#' @param fixed \[integer(1)]\cr The number of fixed time points.
#' @param n_id \[integer(1)]\cr Number of unique groups.
#' @param n_time \[integer(1)]\cr Number of unique time points.
#' @noRd
clear_nonfixed <- function(newdata, newdata_null, resp_stoch,
                           group_var, clear_names, fixed, n_id, n_time) {

  if (newdata_null) {
    predict_idx <- rep(seq.int(fixed + 1L, n_time), n_id) +
      rep(seq.int(0L, n_id - 1L) * n_time, each = n_time - fixed)
    newdata[predict_idx, c(resp_stoch) := NA]
    newdata_names <- names(newdata)
    for (name in newdata_names) {
      if (name %in% clear_names) {
        newdata[ , (name) := NULL]
      }
    }
  }
}

#' Generate Random Intercepts for New Group Levels
#'
#' @inheritParams predict
#' @param nu Posterior draws of the random intercept parameters.
#' @param sigma_nu Posterior draws of the standard deviations of the random
#'   intercept parameters.
#' @param corr_matrix_nu Posterior draws of the correlation matrix of
#'   intercepts (within-group).
#' @param orig_ids Group levels of the original data.
#' @param new_ids  Group levels of `newdata` in
#'   [dynamite::predict.dynamitefit()].
#' @param new_levels \[`character(1)`]\cr
#'   Defines if and how to sample the random intercepts for observations whose
#'   group level was not present in the original data. The options are
#'     * `"none"` (the default) which will signal an error if new levels
#'       are encountered.
#'     * `"bootstrap"` which will randomly draw from the posterior samples of
#'       the random intercepts across all original levels.
#'     * `"gaussian"` which will randomly draw from a gaussian
#'       distribution using the posterior samples of the random intercept
#'       standard deviation (and correlation matrix if applicable).
#'     * `"original"` which will randomly match each new level to one of
#'       the original levels. The posterior samples of the random intercept of
#'       the matched levels will then be used for the new levels.
#' @return An n_draws x n_groups x n_intercepts array of random intercepts.
#' @noRd
generate_random_intercept <- function(nu, sigma_nu, corr_matrix_nu, n_draws,
                                      n_id, orig_ids, new_ids, new_levels) {
  is_orig <- which(orig_ids %in% new_ids)
  is_new <- !new_ids %in% orig_ids
  out <- NULL
  if (identical(new_levels, "none")) {
    out <- nu[, is_orig, , drop = FALSE]
  } else {
    M <- nrow(sigma_nu)
    out <- array(0.0, c(n_draws, n_id, M))
    out[, which(!is_new), ] <- nu[, is_orig, , drop = FALSE]
    if (any(is_new)) {
      n_new <- sum(is_new)
      out[ , which(is_new), ] <- switch(new_levels,
        `bootstrap` = {
          idx <- sample.int(n_draws * n_id, n_draws * n_new, TRUE)
          array(matrix(nu, ncol = M)[idx, ], c(n_draws, n_new, M))
        },
        `gaussian` = {
          x <- array(0, c(n_new, M, n_draws))
          zeros <- rep(0, M)
          if (is.null(corr_matrix_nu)) {
            # easy to optimise...
            for(i in seq_len(n_draws)) {
              s <- diag(sigma_nu[, i]^2, M)
              x[, , i] <- MASS::mvrnorm(n_new, zeros, s)
            }
          } else {
            # Could also keep the Cholesky L from the sampling phase if this is
            # too slow, or switch algorithm. But probably no need as this is
            # only done once
            for(i in seq_len(n_draws)) {
              s <- diag(sigma_nu[, i])
              x[, , i] <- MASS::mvrnorm(
                n_new, zeros, s %*% corr_matrix_nu[, ,i] %*% s
              )
            }
          }
          aperm(x, c(3, 1, 2))
        },
        `original` = {
          match_ids <- sample(orig_ids, n_new, replace = TRUE)
          nu[, match_ids, , drop = FALSE]
        }
      )
    }
  }
  out
}

#' Prepare Environments to Evaluate Predictions or Fitted Values
#'
#' @inheritParams predict
#' @inheritParams clear_nonfixed
#' @param eval_type \[`character(1)`]\cr
#'   Either `"predict"` or `"fitted"`
#' @param predict_type \[`character()`]\cr
#'   Either `"response"`, `"mean"`, or `"link"`.
#' @noRd
prepare_eval_envs <- function(object, newdata, eval_type, predict_type,
                              resp_stoch, n_id, n_draws,
                              new_levels, group_var) {
  samples <- rstan::extract(object$stanfit)
  model_vars <- object$stan$model_vars
  newdata_names <- names(newdata)
  n_resp <- length(resp_stoch)
  eval_envs <- vector(mode = "list", length = n_resp)
  idx_draws <- seq_len(n_draws)
  nu_channels <- attr(object$dformulas$stoch, "random")$channels
  M <- length(nu_channels)
  if (!is.null(group_var) && M > 0) {
    orig_ids <- unique(object$data[[group_var]])
    new_ids <-  unique(newdata[[group_var]])
    n_all_draws <- ndraws(object)
    sigma_nus <- glue::glue("sigma_nu_{nu_channels}")
    sigma_nu <- t(do.call("cbind",
      samples[sigma_nus])[idx_draws, , drop = FALSE])
    nus <- glue::glue("nu_{nu_channels}")
    nu_samples <- array(unlist(samples[nus]),
      c(n_all_draws, n_id, M))[idx_draws, , , drop = FALSE]
    if (attr(object$dformulas$stoch, "random")$correlated) {
      corr_matrix_nu <- aperm(
        samples[["corr_matrix_nu"]][idx_draws, , , drop = FALSE])
    } else {
      corr_matrix_nu <- NULL
    }
    nu_samples <- generate_random_intercept(
      nu = nu_samples,
      sigma_nu = sigma_nu,
      corr_matrix_nu = corr_matrix_nu,
      n_draws = n_draws,
      n_id = n_id,
      orig_ids = orig_ids,
      new_ids = new_ids,
      new_levels = new_levels
    )
    dimnames(nu_samples)[[3]] <- nus
  }
  for (j in seq_len(n_resp)) {
    resp <- resp_stoch[j]
    resp_family <- object$dformulas$stoch[[j]]$family
    alpha <- paste0("alpha_", resp)
    beta <- paste0("beta_", resp)
    delta <- paste0("delta_", resp)
    phi <- paste0("phi_", resp)
    sigma <- paste0("sigma_", resp)
    e <- new.env()
    e$type <- predict_type
    e$n_id <- n_id
    e$n_draws <- n_draws
    e$k <- n_draws * n_id
    e$newdata <- newdata
    e$J_fixed <- model_vars[[j]]$J_fixed
    e$K_fixed <- model_vars[[j]]$K_fixed
    e$J_varying <- model_vars[[j]]$J_varying
    e$K_varying <- model_vars[[j]]$K_varying
    e$resp <- resp_stoch[j]
    e$phi <- samples[[phi]][idx_draws]
    e$sigma <- samples[[sigma]][idx_draws]
    if (resp %in% nu_channels) {
      e$nu <- matrix(nu_samples[, , paste0("nu_", resp)], ncol = n_id)
    }
    if (is_categorical(resp_family)) {
      resp_levels <- attr(object$stan$responses, "resp_class")[[resp]] |>
        attr("levels")
      e$resp_levels <- resp_levels
      e$S <- length(resp_levels)
      if (model_vars[[j]]$has_fixed_intercept) {
        e$alpha <- samples[[alpha]][idx_draws, , drop = FALSE]
      }
      if (model_vars[[j]]$has_varying_intercept) {
        e$alpha <- samples[[alpha]][idx_draws, , , drop = FALSE]
      }
      e$beta <- samples[[beta]][idx_draws, , , drop = FALSE]
      e$delta <- samples[[delta]][idx_draws, , , , drop = FALSE]
      e$xbeta <- matrix(0.0, e$k, e$S)
    } else {
      resp_levels <- NULL
      if (model_vars[[j]]$has_fixed_intercept) {
        e$alpha <- array(samples[[alpha]][idx_draws], c(n_draws, 1L))
      }
      if (model_vars[[j]]$has_varying_intercept) {
        e$alpha <- samples[[alpha]][idx_draws, , drop = FALSE]
      }
      e$beta <- samples[[beta]][idx_draws, , drop = FALSE]
      e$delta <- samples[[delta]][idx_draws, , , drop = FALSE]
      e$xbeta <- numeric(e$k)
    }
    e$call <- generate_sim_call(
      resp,
      resp_levels,
      resp_family,
      eval_type,
      model_vars[[j]]$has_fixed,
      model_vars[[j]]$has_varying,
      model_vars[[j]]$has_fixed_intercept,
      model_vars[[j]]$has_varying_intercept,
      model_vars[[j]]$has_random_intercept,
      model_vars[[j]]$has_offset
    )
    eval_envs[[j]] <- e
  }
  eval_envs
}

#' Generate a Quoted Expression to Evaluate Predictions or Fitted Values
#'
#' @param resp \[`character(1)`]\cr Name of the response.
#' @param resp_levels \[`character()`]\cr Levels of a categorical response.
#' @param family \[`dynamitefamily`]\cr Family of the response.
#' @param type  \[`character(1)`]\cr Either `"predict"` or `"fitted"`.
#' @param has_fixed \[logical(1)]\cr
#'   Does the channel have time-invariant predictors?
#' @param has_varying \[logical(1)]\cr
#'   Does the channel have time-varying predictors?
#' @param has_fixed_intercept \[logical(1)]\cr
#'   Does the channel have a time-invariant intercept?
#' @param has_varying_intercept \[logical(1)]\cr
#'   Does the channel have a time-varying intercept?
#' @param has_random_intercept \[logical(1)]\cr
#'   Does the channel have a random intercept?
#' @param has_offset \[logical(1)]\cr
#'   Does the channel have an offset?
#' @noRd
generate_sim_call <- function(resp, resp_levels, family, type,
                              has_fixed, has_varying,
                              has_fixed_intercept, has_varying_intercept,
                              has_random_intercept, has_offset) {
  if (is_categorical(family)) {
    glue::glue(
      "{{\n",
      paste0(
        "for (j in seq_len(n_id)) {{\n",
        "  idx_draw <- seq.int((j - 1L) * n_draws + 1L, j * n_draws)\n",
        "  for (s in seq_len(S - 1)) {{\n",
        "    xbeta[idx_draw, s + 1] <- ",
        "{ifelse_(!has_fixed_intercept && !has_varying_intercept, '0', '')}",
        "{ifelse_(has_fixed_intercept, 'alpha[, s]', '')}",
        "{ifelse_(has_varying_intercept, 'alpha[, a_time, s]', '')}",
        "{ifelse_(has_fixed, ",
        "' + .rowSums(x = model_matrix[idx_draw, J_fixed, drop = FALSE] ",
                        " * beta[, , s], ",
                        " m = n_draws, n = K_fixed)', '')}",
        "{ifelse_(has_varying, ",
        "' + .rowSums(x = model_matrix[idx_draw, J_varying, drop = FALSE] ",
        " * delta[, time, , s],",
        "m = n_draws, n = K_varying)', '')}",
        "}}\n",
        "}}\n"
      ),
      eval(str2lang(glue::glue("{type}_categorical"))),
      "}}"
    ) |> str2lang()
  } else {
    glue::glue(
      "{{\n",
      paste0(
        "for (j in seq_len(n_id)) {{\n",
        "  idx_draw <- seq.int((j - 1L) * n_draws + 1L, j * n_draws)\n",
        "  xbeta[idx_draw] <- ",
        "{ifelse_(!has_fixed_intercept && !has_varying_intercept, '0', '')}",
        "{ifelse_(has_fixed_intercept, 'alpha', '')}",
        "{ifelse_(has_varying_intercept, 'alpha[, a_time]', '')}",
        "{ifelse_(has_random_intercept, '+ nu[, j]', '')}",
        "{ifelse_(has_fixed, ",
        "' + .rowSums(x = model_matrix[idx_draw, J_fixed, drop = FALSE] * beta,
                        m = n_draws, n = K_fixed)', '')}",
        "{ifelse_(has_varying, ",
          "' + .rowSums(x = model_matrix[idx_draw, J_varying, drop = FALSE] ",
                        " * delta[, time, ],",
                        "m = n_draws, n = K_varying)', '')}",
        "}}\n"
      ),
      ifelse_(identical(type, "predict"),
        paste0(
          "if (type == 'link') {{",
          "  data.table::set(x = newdata, i = idx_data, j = '{resp}_link',
                             value = xbeta[idx_pred])",
          "}}"
        ),
        ""
      ),
      eval(str2lang(glue::glue("{type}_{family}"))),
      "}}"
    ) |> str2lang()
  }
}

# Fitted expressions ------------------------------------------------------

fitted_gaussian <- "
  data.table::set(x = newdata, i = idx, j = '{resp}_fitted', value = xbeta)
"

fitted_categorical <- "
    resp_cols <- c({paste0('\"', resp, '_fitted_', resp_levels, '\"',
                    collapse = ', ')})
    maxs <- apply(xbeta, 1, max)
    mval <- exp(xbeta - (maxs + log(rowSums(exp(xbeta - maxs)))))
    for (s in 1:S) {{
      data.table::set(
        x = newdata,
        i = idx,
        j = resp_cols[s],
        value = mval[, s]
      )
    }}
"

fitted_bernoulli <- "
  data.table::set(
    x = newdata,
    i = idx,
    j = '{resp}_fitted',
    value = plogis(xbeta)
  )
"

fitted_binomial <- "
  data.table::set(
    x = newdata,
    i = idx,
    j = '{resp}_fitted',
    value = plogis(xbeta)
  )
"

fitted_poisson <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  data.table::set(x = newdata, i = idx, j = '{resp}_fitted', value = exp_xbeta)
"

fitted_negbin <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  data.table::set(x = newdata, i = idx, j = '{resp}_fitted', value = exp_xbeta)
"

fitted_exponential <- "
  data.table::set(x = newdata, i = idx, j = '{resp}_fitted', value = exp(xbeta))
"

fitted_gamma <- "
  data.table::set(x = newdata, i = idx, j = '{resp}_fitted', value = exp(xbeta))
"

fitted_beta <- "
  data.table::set(
    x = newdata,
    i = idx,
    j = '{resp}_fitted',
    value = plogis(xbeta)
  )
"
# Predict expressions -----------------------------------------------------

predict_gaussian <- "
  if (type == 'mean') {{
    data.table::set(
      x = newdata,
      i = idx_data,
      j = '{resp}_mean',
      value = xbeta[idx_pred]
    )
  }}
  data.table::set(
    x = newdata,
    i = idx_data,
    j = '{resp}',
    value = rnorm(k, xbeta, sigma)[idx_pred]
  )
"

predict_categorical <- "
  if (type == 'link') {{
    resp_cols <- c({paste0('\"', resp, '_link_', resp_levels, '\"',
                           collapse = ', ')})
    for (s in 1:S) {{
      data.table::set(
        x = newdata,
        i = idx_data,
        j = resp_cols[s],
        value = xbeta[idx_pred, s]
      )
    }}
  }}
  if (type == 'mean') {{
    resp_cols <- c({paste0('\"', resp, '_mean_', resp_levels, '\"',
                           collapse = ', ')})
    maxs <- apply(xbeta, 1, max)
    mval <- exp(xbeta - (maxs + log(rowSums(exp(xbeta - maxs)))))
    for (s in 1:S) {{
      data.table::set(
        x = newdata,
        i = idx_data,
        j = resp_cols[s],
        value = mval[idx_pred, s]
      )
    }}
  }}
  data.table::set(
    x = newdata,
    i = idx_data,
    j = '{resp}',
    value = max.col(xbeta - log(-log(runif(S * k))))[idx_pred]
  )
"

predict_binomial <- "
  prob <- plogis(xbeta)
  if (type == 'mean') {{
    data.table::set(
      x = newdata,
      i = idx_data,
      j = '{resp}_mean',
      value = prob[idx_pred]
    )
  }}
  data.table::set(
    x = newdata,
    i = idx_data,
    j = '{resp}',
    value = rbinom(k, trials, prob)[idx_pred]
  )
"

predict_bernoulli <- "
  prob <- plogis(xbeta)
  if (type == 'mean') {{
    data.table::set(
      x = newdata,
      i = idx_data,
      j = '{resp}_mean',
      value = prob[idx_pred]
    )
  }}
  data.table::set(
    x = newdata,
    i = idx_data,
    j = '{resp}',
    value = rbinom(k, 1, prob)[idx_pred]
  )
"

predict_poisson <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  if (type == 'mean') {{
    data.table::set(
      x = newdata,
      i = idx_data,
      j = '{resp}_mean',
      value = exp_xbeta[idx_pred]
    )
  }}
  data.table::set(
    x = newdata,
    i = idx_data,
    j = '{resp}',
    value = rpois(k, exp_xbeta)[idx_pred]
  )
"

predict_negbin <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  if (type == 'mean') {{
    data.table::set(
      x = newdata,
      i = idx_data,
      j = '{resp}_mean',
      value = exp_xbeta[idx_pred]
    )
  }
  data.table::set(
    x = newdata,
    i = idx_data,
    j = '{resp}',
    value = rnbinom(k, size = phi, mu = exp_xbeta)[idx_pred]
  )
"

predict_exponential <- "
  if (type == 'mean') {{
    data.table::set(
      x = newdata,
      i = idx_data,
      j = '{resp}_mean',
      value = exp(xbeta[idx_pred])
    )
  }
  data.table::set(
    x = newdata,
    i = idx_data,
    j = '{resp}',
    value = rexp(k, rate = exp(-xbeta))[idx_pred]
  )
"

predict_gamma <- "
  if (type == 'mean') {{
    data.table::set(
      x = newdata,
      i = idx_data,
      j = '{resp}_mean',
      value = exp(xbeta[idx_pred])
    )
  }
  data.table::set(
    x = newdata,
    i = idx_data,
    j = '{resp}',
    value = rgamma(k, shape = phi, rate = phi * exp(-xbeta))[idx_pred]
  )
"

predict_beta <- "
  mu <- {'plogis(xbeta)'}
  if (type == 'mean') {{
    data.table::set(
      x = newdata,
      i = idx_data,
      j = '{resp}_mean',
      value = mu[idx_pred])
    )
  }
  data.table::set(
    x = newdata,
    i = idx_data,
    j = '{resp}',
    value = rbeta(k,  mu * phi, (1 - mu) * phi)[idx_pred]
  )
"

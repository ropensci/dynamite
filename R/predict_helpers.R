#' Check Validity of `n_draws` Argument for Prediction
#'
#' @param n_draws \[`integer(1)`]\cr
#'   Number of draws to use for `fitted` or `predict`.
#' @param full_draws \[`integer(1)`]\cr
#'   Number of draws available in the `dynamitefit` object.
#' @noRd
check_ndraws <- function(n_draws, full_draws) {
  n_draws <- ifelse_(is.null(n_draws), full_draws, n_draws)
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

#' Check Validity of `newdata` Argument for Prediction
#'
#' @inheritParams predict.dynamitefit
#' @noRd
check_newdata <- function(object, newdata) {
  if (is.null(newdata)) {
    data.table::copy(object$data)
  } else if (data.table::is.data.table(newdata) || is.data.frame(newdata)) {
    data.table::as.data.table(data.table::copy(newdata))
  } else if (!is.data.frame(newdata)) {
    stop_("Argument {.arg newdata} must be a {.cls data.frame} object.")
  }
}

#' Parse And Prepare `newdata` for Prediction
#'
#' @param dformulas \[`dynamiteformula`]\cr The model formulas.
#' @param newdata \[`data.frame`]\cr The data to be used for prediction.
#' @param data \[`data.frame`]\cr The original data used to fit the model.
#' @param type \[`character`]\cr Either `"response"`, `"mean"`, `"link"`.
#' @param eval_type \[`character(1)`]\cr Either `"predict"`, `"fitted"`, or
#'   `loglik.`
#' @param families_stoch \[`list`]\cr of `dynamitefamily` object.
#' @param resp_stoch \[`character`]\cr A vector of response variables
#' @param categories \[`list`]\cr Response variable categories
#'   for categorical families.
#' @param clear_names \[character()]\cr Vector of column names to remove
#'   from `newdata`.
#' @param new_levels \[`character(1)`]\cr Either `"none"`, `"bootstrap"`,
#'  `"gaussian"`, or `"original"`.
#' @param group_var \[`character(1)`,`NULL`]\cr Group variable name or `NULL`
#'   if there is only one group.
#' @param time_var \[`character(1)`]\cr Time index variable name.
#' @noRd
parse_newdata <- function(dformulas, newdata, data, type, eval_type,
                          resp_stoch, categories, clear_names, new_levels,
                          group_var, time_var) {
  if (!group_var %in% names(newdata)) {
    orig <- sort(data[[group_var]])[1L]
    data.table::set(x = newdata, j = group_var, value = orig)
  }
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
  missing_resp <- resp_stoch[!resp_stoch %in% names(newdata)]
  stopifnot_(
    identical(length(missing_resp), 0L),
    "Can't find response variable{?s} {.var {missing_resp}} in {.var newdata}."
  )
  # check and add missing factor levels
  factor_cols <- setdiff(
    names(which(vapply(data, is.factor, logical(1L)))),
    c(time_var, group_var)
  )
  cols <- intersect(names(newdata), factor_cols)
  for (i in cols) {
    l_orig <- levels(data[[i]])
    l_new <- levels(newdata[[i]])
    if (!setequal(l_orig, l_new)) {
      stopifnot_(
        all(l_new %in% l_orig),
        c(
          "{.cls factor} variable {.var {i}} in {.arg newdata} has new levels:",
          `x` = "Level{?s} {.val {setdiff(l_new, l_orig)}} {?is/are}
                 not present in the original data."
        )
      )
      newdata[[i]] <- factor(newdata[[i]], levels = l_orig)
    }
  }
  newdata <- fill_time_predict(
    newdata,
    group_var,
    time_var,
    time_scale = original_times[2L] - original_times[1L]
  )
  data.table::setDT(newdata, key = c(group_var, time_var))
  clear_names <- intersect(names(newdata), clear_names)
  if (length(clear_names) > 0L) {
    # TODO no need check length when data.table package is updated
    newdata[, (clear_names) := NULL]
  }
  drop_unused(dformulas$all, newdata, group_var, time_var)
  type <- ifelse_(eval_type %in% c("fitted", "loglik"), eval_type, type)
  if (identical(type, "loglik")) {
    cg <- attr(dformulas$stoch, "channel_groups")
    n_cg <- n_unique(cg)
    for (i in seq_len(n_cg)) {
      cg_idx <- which(cg == i)
      y <- ifelse_(
        length(cg_idx) > 1L,
        paste(c(resp_stoch[cg_idx], "loglik"), collapse = "_"),
        paste0(resp_stoch[cg_idx[1L]], "_loglik")
      )
      newdata[, (y) := NA_real_]
    }
  }
  for (i in seq_along(resp_stoch)) {
    y <- resp_stoch[i]
    if (type %in% c("mean", "link", "fitted")) {
      # create a separate column for each level of
      # a categorical response variables
      pred_col <- ifelse_(
        is_categorical(dformulas$stoch[[i]]$family),
        glue::glue("{y}_{type}_{categories[[y]]}"),
        glue::glue("{y}_{type}")
      )
      newdata[, (pred_col) := NA_real_]
    }
    data.table::set(
      x = newdata,
      j = glue::glue("{y}_store"),
      value = newdata[[y]]
    )
  }
  newdata
}

#' Parse `funs` Argument for Prediction
#'
#' @inheritParams predict.dynamitefit
#' @noRd
parse_funs <- function(object, type, funs, categories) {
  stopifnot_(
    !is.null(object$group_var),
    "Argument {.arg funs} requires data with multiple individuals."
  )
  stopifnot_(
    is.list(funs),
    "Argument {.arg funs} must be a {.cls list}."
  )
  funs_names <- names(funs)
  stopifnot_(
    !is.null(funs_names),
    "Argument {.arg funs} must be named."
  )
  valid_names <- get_responses(object$dformulas$all)
  stopifnot_(
    all(funs_names %in% get_responses(object$dformulas$all)),
    "The names of {.arg funs} must be response variables of the model."
  )
  out <- list()
  idx <- 1L
  suffix <- ifelse_(
    identical(type, "response"),
    "",
    paste0("_", type)
  )
  resp_stoch <- get_responses(object$dformulas$stoch)
  family_stoch <- get_families(object$dformulas$stoch)
  for (i in seq_along(funs)) {
    stopifnot_(
      is.list(funs[[i]]),
      "Each element of {.arg funs} must be a {.cls list}."
    )
    fun_names <- names(funs[[i]])
    stopifnot_(
      !is.null(fun_names),
      "Each element of {.arg funs} must be named."
    )
    for (j in seq_along(funs[[i]])) {
      stopifnot_(
        is.function(funs[[i]][[j]]),
        "Each element of {.arg funs} must contain only functions."
      )
      resp_idx <- funs_names[i] %in% resp_stoch
      target <- funs_names[i]
      name <- paste0(fun_names[j], "_", funs_names[i])
      if (length(resp_idx) > 0L) {
        category <- ifelse_(
          !identical(type, "response") &&
            is_categorical(family_stoch[[resp_idx]]),
          paste0("_", categories[[funs_names[i]]]),
          ""
        )
        target <- paste0(target, suffix, category)
        name <- paste0(name, category)
      }
      for (k in seq_along(name)) {
        out[[idx]] <- list(
          fun = funs[[i]][[j]],
          name = name[k],
          target = target[k]
        )
        idx <- idx + 1L
      }
    }
  }
  out
}

#' Adds NA Gaps to Fill Missing Time Points in a Data Frame For Predictions
#'
#' Note that if `impute` is `none` and model contains lagged predictors,
#' predictions will eventually fail.
#'
#' @inheritParams dynamite
#' @noRd
fill_time_predict <- function(data, group_var, time_var, time_scale) {
  # time_duplicated <- data[,
  #  any(duplicated(time_var)),
  #  by = group_var,
  #  env = list(time_var = time_var)
  # ]$V1
  time <- sort(unique(data[[time_var]]))
  time_ivals <- diff(time)
  time_scale <- min(diff(time))
  full_time <- seq(time[1L], time[length(time)], by = time_scale)
  data_groups <- as.integer(data[[group_var]])
  group <- unique(data_groups)
  n_group <- length(group)
  time_duplicated <- logical(n_group)
  time_groups <- list(
    has_missing = logical(n_group),
    has_gaps = logical(n_group)
  )
  time_missing <- logical(n_group)
  for (i in seq_len(n_group)) {
    idx_group <- which(data_groups == group[i])
    sub <- data[idx_group, ]
    time_duplicated[i] <- any(duplicated(sub[[time_var]]))
    time_groups$has_missing[i] <- !identical(sub[[time_var]], full_time)
    time_groups$has_gaps[i] <- length(idx_group) !=
      (diff(range(sub[[time_var]])) + 1L) * time_scale
  }
  d <- which(time_duplicated)
  stopifnot_(
    all(!time_duplicated),
    c(
      "Each time index must correspond to a single observation per group:",
      `x` = "{cli::qty(length(d))}Group{?s} {.var {d}} of {.var {group_var}}
             {cli::qty(length(d))}{?has/have} duplicate observations."
    )
  )
  if (length(time) > 1L) {
    original_order <- colnames(data)
    # time_groups <- data[,
    # {
    #  has_missing = !identical(time_var, full_time)
    #  has_gaps = .N != (diff(range(time_var)) + 1L) * time_scale
    #  list(has_missing, has_gaps)
    # },
    # by = group_var
    # env = list(
    #  time_var = time_var,
    #  group_var = group_var,
    #  time_scale = time_scale
    # )
    # ]
    if (any(time_groups$has_missing)) {
      if (any(time_groups$has_gaps)) {
        warning_(c(
          "Time index variable {.var {time_var}} of {.arg newdata} has gaps:",
          `i` = "Filling the {.arg newdata} to regular time points. This will
                 lead to propagation of NA values if the model contains
                 exogenous predictors and {.arg impute} is {.val none}."
        ))
      }
      full_data_template <- data.table::as.data.table(expand.grid(
        time = full_time,
        group = unique(data[[group_var]])
      ))
      names(full_data_template) <- c(time_var, group_var)
      data <- data.table::merge.data.table(
        full_data_template,
        data,
        by = c(time_var, group_var),
        all.x = TRUE
      )
    }
  }
  data
}

#' Impute Predictor Values in `newdata`
#'
#' @inheritParams predict
#' @param predictors \[`character()`]\cr A vector of predictor column names.
#' @param group_var \[`character(1)`]\cr Grouping variable name.
#' @noRd
impute_newdata <- function(newdata, impute, predictors, group_var) {
  switch(impute,
    `locf` = {
      newdata[,
        (predictors) := lapply(.SD, locf),
        .SDcols = predictors,
        by = group_var
        # by = group_var,
        # env = list(locf = locf)
      ]
    },
    `nocb` = {
      newdata[,
        (predictors) := lapply(.SD, nocb),
        .SDcols = predictors,
        by = group_var
        # by = group_var,
        # env = list(locf = locf)
      ]
    },
    `none` = {
      newdata
    }
  )
}

#' Assign NA Values to Time Indices Beyond Fixed Time Points
#'
#' @inheritParams parse_newdata
#' @param newdata_null \[logical(1)]\cr
#'   Was `newdata` `NULL` when `predict` was called?
#' @param resp_stoch \[character()]\cr Vector of response variable names.
#' @param fixed \[integer(1)]\cr The number of fixed time points.
#' @noRd
clear_nonfixed <- function(newdata, newdata_null, resp_stoch, eval_type,
                           group_var, time_var, clear_names,
                           fixed, global_fixed) {
  if (newdata_null && identical(eval_type, "predicted")) {
    if (global_fixed) {
      clear_idx <- newdata[,
        .I[seq.int(fixed + 1L, .N)],
        by = group_var
        # by = group_var,
        # env = list(fixed = fixed)
      ]$V1
    } else {
      clear_idx <- newdata[,
        .I[seq.int(fixed + which(apply(!is.na(.SD), 1L, any))[1L], .N)],
        .SDcols = resp_stoch,
        by = group_var
        # by = group_var,
        # env = list(fixed = fixed, any = any)
      ]$V1
    }
    # newdata[clear_idx, c(resp_stoch) := NA, env = list(clear_idx = clear_idx)]
    newdata[clear_idx, c(resp_stoch) := NA]
  }
}

#' Generate Random Effects for New Group Levels
#'
#' @inheritParams predict
#' @param nu Posterior draws of the random effects.
#' @param sigma_nu Posterior draws of the standard deviations of the random
#'   effects.
#' @param corr_matrix_nu Posterior draws of the correlation matrix of
#'   random effects.
#' @param n_group \[`integer(1)`]\cr Number of groups.
#' @param orig_ids \[`character()`]\cr Group levels of the original data.
#' @param new_ids \[`character()`]\cr Group levels of `newdata` in
#'   [dynamite::predict.dynamitefit()].
#' @param new_levels \[`character(1)`]\cr
#'   Defines if and how to sample the random effects for groups not present in
#'   the original data. The options are
#'     * `"none"` (the default) which will signal an error if new levels
#'       are encountered.
#'     * `"bootstrap"` which will randomly draw from the posterior samples of
#'       the random effects across all original levels.
#'     * `"gaussian"` which will randomly draw from a gaussian
#'       distribution using the posterior samples of the random effect
#'       standard deviations (and correlation matrix if applicable).
#'     * `"original"` which will randomly match each new level to one of
#'       the original levels. The posterior samples of the random effects of
#'       the matched levels will then be used for the new levels.
#' @return An n_draws x n_group x n_intercepts array of random intercepts.
#' @noRd
generate_random_effect <- function(nu, sigma_nu, corr_matrix_nu, n_draws,
                                   n_group, orig_ids, new_ids, new_levels) {
  is_orig <- which(orig_ids %in% new_ids)
  is_new <- !new_ids %in% orig_ids
  out <- NULL
  if (identical(new_levels, "none")) {
    out <- nu[, is_orig, , drop = FALSE]
  } else {
    M <- nrow(sigma_nu)
    out <- array(0.0, c(n_draws, n_group, M))
    out[, which(!is_new), ] <- nu[, is_orig, , drop = FALSE]
    if (any(is_new)) {
      n_new <- sum(is_new)
      out[, which(is_new), ] <- switch(new_levels,
        `bootstrap` = {
          idx <- sample.int(n_draws * length(orig_ids), n_draws * n_new, TRUE)
          array(matrix(nu, ncol = M)[idx, ], c(n_draws, n_new, M))
        },
        `gaussian` = {
          x <- array(0.0, c(n_new, M, n_draws))
          if (is.null(corr_matrix_nu)) {
            for (i in seq_len(n_draws)) {
              x[, , i] <- matrix(
                rnorm(n_new * M, mean = 0.0, sd = sigma_nu[, i]^2),
                nrow = n_new,
                ncol = M,
                byrow = TRUE
              )
            }
          } else {
            # Could also keep the Cholesky L from the sampling phase if this is
            # too slow, or switch algorithm. But probably no need as this is
            # only done once
            for (i in seq_len(n_draws)) {
              s <- diag(sigma_nu[, i])
              e <- eigen(s %*% corr_matrix_nu[, , i] %*% s)
              x[, , i] <- matrix(rnorm(n_new * M), nrow = n_new) %*%
                diag(sqrt(pmax(e$values, 0)), M) %*%
                t(e$vectors)
            }
          }
          aperm(x, c(3L, 1L, 2L))
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

#' Expand the Prediction Data Frame to Include All Covariates
#'
#' @param x Object returned by `predict.dynamitefit`.
#' @noRd
expand_predict_output <- function(simulated, observed, df) {
  add_cols <- setdiff(names(observed), names(simulated))
  observed <- observed[
    rep(seq_len(nrow(observed)), n_unique(simulated$.draw)),
  ]
  for (col in add_cols) {
    data.table::set(x = simulated, j = col, value = observed[[col]])
  }
  if (df) {
    data.table::setDF(simulated)
  }
  simulated
}

#' Prepare Environments to Evaluate Predictions or Fitted Values
#'
#' @inheritParams predict_dynamitefit
#' @inheritParams clear_nonfixed
#' @noRd
prepare_eval_envs <- function(object, simulated, observed,
                              type, eval_type, n_draws,
                              new_levels, group_var) {
  samples <- rstan::extract(object$stanfit)
  channel_vars <- object$stan$channel_vars
  channel_group_vars <- object$stan$channel_group_vars
  cg <- attr(object$dformulas$all, "channel_groups")
  n_cg <- n_unique(cg)
  eval_envs <- vector(mode = "list", length = n_cg)
  idx_draws <- seq_len(n_draws)
  rand <- which_random(object$dformulas$all)
  n_group <- n_unique(observed[[group_var]])
  k <- 0L # index of channel_vars
  l <- 0L # index of channel_group_vars
  orig_ids <- unique(object$data[[group_var]])
  new_ids <- unique(observed[[group_var]])
  extra_levels <- unique(new_ids[!new_ids %in% orig_ids])
  has_lfactor <- attr(object$dformulas$stoch, "lfactor")$P > 0
  stopifnot_(identical(length(extra_levels), 0L) || !has_lfactor,
    c(
      "Grouping variable {.var {group_var}} contains unknown levels:",
      `x` = "Level{?s} {.val {as.character(extra_levels)}}
             {?is/are} not present in the original data.",
      `i` = "Models with latent factors do not support new levels because of
             identifiability constraints."
    )
  )
  if (length(rand) > 0L) {
    n_all_draws <- ndraws(object)
    Ks <- object$stan$model_vars$Ks
    nu_channels <- names(Ks[Ks > 0L])
    sigma_nus <- glue::glue("sigma_nu_{nu_channels}")
    sigma_nu <- t(
      do.call("cbind", samples[sigma_nus])[idx_draws, , drop = FALSE]
    )
    nus <- glue::glue("nu_{nu_channels}")
    M <- nrow(sigma_nu)
    nu_samples <- array(
      unlist(samples[nus]),
      c(n_all_draws, length(orig_ids), M)
    )[idx_draws, , , drop = FALSE]
    corr_matrix_nu <- onlyif(
      attr(object$dformulas$stoch, "random_spec")$correlated,
      aperm(
        samples[["corr_matrix_nu"]][idx_draws, , , drop = FALSE]
      )
    )
    nu_samples <- generate_random_effect(
      nu = nu_samples,
      sigma_nu = sigma_nu,
      corr_matrix_nu = corr_matrix_nu,
      n_draws = n_draws,
      n_group = n_group,
      orig_ids = orig_ids,
      new_ids = new_ids,
      new_levels = new_levels
    )
    dimnames(nu_samples)[[3L]] <- make.unique(rep(nus, times = Ks[Ks > 0L]))
  }
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    family <- object$dformulas$all[[cg_idx[1L]]]$family
    if (is_deterministic(family)) {
      eval_envs[[i]] <- list()
      next
    }
    l <- l + 1L
    e <- new.env()
    e$family <- family
    e$out <- simulated
    e$type <- type
    e$n_group <- n_group
    e$n_draws <- n_draws
    e$k <- n_draws * n_group
    if (is_multivariate(family) || is_categorical(family)) {
      k <- k + length(cg_idx)
      resp <- get_responses(object$dformulas$all[cg_idx])
      resp_levels <- ifelse_(
        is_categorical(family),
        attr(
          attr(object$stan$responses, "resp_class")[[resp]],
          "levels"
        ),
        resp
      )
      prepare_eval_env_multivariate(
        e = e,
        resp = resp,
        resp_levels = resp_levels,
        cvars = channel_vars[cg_idx],
        samples = samples,
        nu_samples = nu_samples,
        has_random_effects = cg_idx %in% rand,
        idx = idx_draws,
        type = type,
        eval_type = eval_type
      )
    } else {
      j <- cg_idx[1L]
      k <- k + 1L
      resp <- object$dformulas$all[[j]]$response
      prepare_eval_env_univariate(
        e = e,
        resp = resp,
        cvars = channel_vars[[k]],
        samples = samples,
        nu_samples = nu_samples,
        has_random_effects = j %in% rand,
        idx = idx_draws,
        type = type,
        eval_type = eval_type
      )
    }
    eval_envs[[i]] <- e
  }
  eval_envs
}

#' Prepare a Evaluation Environment for a Univariate Channel
#'
#' @noRd
prepare_eval_env_univariate <- function(e, resp, cvars, samples, nu_samples,
                                        has_random_effects,
                                        idx, type, eval_type) {
  alpha <- paste0("alpha_", resp)
  beta <- paste0("beta_", resp)
  delta <- paste0("delta_", resp)
  phi <- paste0("phi_", resp)
  sigma <- paste0("sigma_", resp)
  lambda <- paste0("lambda_", resp)
  psi <- paste0("psi_", resp)
  e$J_fixed <- cvars$J_fixed
  e$K_fixed <- cvars$K_fixed
  e$J_varying <- cvars$J_varying
  e$K_varying <- cvars$K_varying
  e$J_random <- cvars$J_random
  e$K_random <- cvars$K_random
  e$has_random_intercept <- cvars$has_random_intercept
  e$has_lfactor <- cvars$has_lfactor
  e$resp <- resp
  e$phi <- onlyif(
    !is.null(samples[[phi]]),
    rep_len(c(samples[[phi]][idx]), length.out = e$k)
  )
  e$sigma <-  onlyif(
    !is.null(samples[[sigma]]),
    rep_len(c(samples[[sigma]][idx]), length.out = e$k)
  )
  if (has_random_effects) {
    nus <- make.unique(rep(paste0("nu_", resp), e$K_random))
    e$nu <- nu_samples[, , nus, drop = FALSE]
  }
  if (cvars$has_fixed_intercept) {
    e$alpha <- array(samples[[alpha]][idx], c(e$n_draws, 1L))
  }
  if (cvars$has_varying_intercept) {
    e$alpha <- samples[[alpha]][idx, , drop = FALSE]
  }
  if (cvars$has_lfactor) {
    e$lambda <- samples[[lambda]][idx, , drop = FALSE]
    e$psi <- samples[[psi]][idx, , drop = FALSE]
  }
  e$beta <- samples[[beta]][idx, , drop = FALSE]
  e$delta <- samples[[delta]][idx, , , drop = FALSE]
  e$xbeta <- numeric(e$k)
  e$call <- generate_sim_call_univariate(
    resp = resp,
    family = e$family,
    type = type,
    eval_type = eval_type,
    has_fixed = cvars$has_fixed,
    has_varying = cvars$has_varying,
    has_random = cvars$has_random,
    has_fixed_intercept = cvars$has_fixed_intercept,
    has_varying_intercept = cvars$has_varying_intercept,
    has_random_intercept = cvars$has_random_intercept,
    has_offset = cvars$has_offset,
    has_lfactor = cvars$has_lfactor
  )
}

#' Generate an Expression to Evaluate Predictions for a Univariate Channel
#'
#' @param resp \[`character(1)`]\cr Name of the response.
#' @param resp_levels \[`character()`]\cr Levels of a categorical response.
#' @param family \[`dynamitefamily`]\cr Family of the response.
#' @param eval_type  \[`character(1)`]\cr Either `"predict"`, `"fitted"` or
#'   `"loglik"`.
#' @param has_fixed \[logical(1)]\cr
#'   Does the channel have time-invariant predictors?
#' @param has_varying \[logical(1)]\cr
#'   Does the channel have time-varying predictors?
#' @param has_random \[logical(1)]\cr
#'   Does the channel have random effects?
#' @param has_fixed_intercept \[logical(1)]\cr
#'   Does the channel have a time-invariant intercept?
#' @param has_varying_intercept \[logical(1)]\cr
#'   Does the channel have a time-varying intercept?
#' @param has_random_intercept \[logical(1)]\cr
#'   Does the channel have a random intercept?
#' @param has_offset \[logical(1)]\cr
#'   Does the channel have an offset?
#' @param has_lfactor \[logical(1)]\cr
#'   Does the channel have a latent factor term?
#' @noRd
generate_sim_call_univariate <- function(resp, family, type, eval_type,
                                         has_fixed, has_varying, has_random,
                                         has_fixed_intercept,
                                         has_varying_intercept,
                                         has_random_intercept,
                                         has_offset, has_lfactor) {

  out <- paste0(
    "{\n",
    "idx_draw <- seq.int(1L, n_draws) - n_draws\n",
    "for (j in seq_len(n_group)) {\n",
    "  idx_draw <- idx_draw + n_draws\n",
    "  xbeta[idx_draw] <- ",
    ifelse_(!has_fixed_intercept && !has_varying_intercept, "0", ""),
    ifelse_(has_fixed_intercept, "alpha", ""),
    ifelse_(has_varying_intercept, "alpha[, a_time]", ""),
    ifelse_(has_random_intercept, " + nu[, j, 1]", ""),
    ifelse_(has_lfactor, " + lambda[, j] * psi[, time]", ""),
    ifelse_(
      has_fixed,
      " + .rowSums(
        x = model_matrix[idx_draw, J_fixed, drop = FALSE] * beta,
        m = n_draws,
        n = K_fixed
      )",
      ""
    ),
    ifelse_(
      has_varying,
      " + .rowSums(
        x = model_matrix[idx_draw, J_varying, drop = FALSE] *
              delta[, time, ],
        m = n_draws,
        n = K_varying
      )",
      ""
    ),
    ifelse_(
      has_random,
      " + .rowSums(
        x = model_matrix[idx_draw, J_random, drop = FALSE] *
              nu[, j, seq.int(1 + has_random_intercept, K_random)],
        m = n_draws,
        n = K_random - has_random_intercept
      )",
      ""
    ),
    "}\n",
    ifelse_(
      identical(type, "link") && identical(eval_type, "predicted"),
      glue::glue(
        "data.table::set(",
        "  x = out,",
        "  i = idx_data,",
        "  j = '{resp}_link',",
        "  value = xbeta[idx_out]",
        ")",
      ),
      ""
    ),
    "\n",
    glue::glue(predict_expr[[eval_type]][[family$name]]),
    "\n",
    ifelse_(
      identical(type, "mean") && identical(eval_type, "predicted"),
      glue::glue(predict_expr$mean[[family$name]]),
      ""
    ),
    "}"
  )
  str2lang(out)
}

#' Prepare a Evaluation Environment for a Multivariate Channel
#'
#' @noRd
prepare_eval_env_multivariate <- function(e, resp, resp_levels, cvars,
                                          samples, nu_samples,
                                          has_random_effects,
                                          idx, type, eval_type) {
  dims <- ifelse_(
    is_multivariate(e$family) && !is_multinomial(e$family),
    seq_len(length(resp)),
    seq.int(2L, length(resp_levels))
  )
  d <- ifelse_(
    is_multivariate(e$family) && !is_multinomial(e$family),
    length(dims),
    length(resp_levels)
  )
  e$dims <- dims
  e$d <- d
  e$resp <- resp
  e$L <- samples[[paste(c("L", resp), collapse = "_")]]
  e$link_cols <- ifelse_(
    is_multivariate(e$family) && !is_multinomial(e$family),
    paste0(resp, "_link"),
    paste0(resp, "_link_", resp_levels)
  )
  e$mean_cols <- paste0(resp, "_mean_", resp_levels)
  e$fitted_cols <- paste0(resp, "_fitted_", resp_levels)
  e$sigma <- matrix(0.0, e$n_draws, d)
  has_fixed <- logical(d)
  has_varying <- logical(d)
  has_random <- logical(d)
  has_fixed_intercept <- logical(d)
  has_varying_intercept <- logical(d)
  has_random_intercept <- logical(d)
  has_offset <- logical(d)
  has_lfactor <- logical(d)
  for (i in dims) {
    yi <- ifelse_(
      is_categorical(e$family),
      paste0(resp, "_", resp_levels[i]),
      resp[i]
    )
    j <- ifelse_(is_categorical(e$family) || is_multinomial(e$family), 1L, i)
    alpha <- paste0("alpha_", yi)
    beta <- paste0("beta_", yi)
    delta <- paste0("delta_", yi)
    phi <- paste0("phi_", yi)
    nu <- paste0("nu_", yi)
    sigma <- paste0("sigma_", yi)
    lambda <- paste0("lambda_", yi)
    psi <- paste0("psi_", yi)
    J_fixed <- paste0("J_fixed_", yi)
    K_fixed <- paste0("K_fixed_", yi)
    J_varying <- paste0("J_varying_", yi)
    K_varying <- paste0("K_varying_", yi)
    J_random <- paste0("J_random_", yi)
    K_random <- paste0("K_random_", yi)
    has_fixed[i] <- cvars[[j]]$has_fixed
    has_varying[i] <- cvars[[j]]$has_varying
    has_random[i] <- cvars[[j]]$has_random
    has_fixed_intercept[i] <- cvars[[j]]$has_fixed_intercept
    has_varying_intercept[i] <- cvars[[j]]$has_varying_intercept
    has_random_intercept[i] <- cvars[[j]]$has_random_intercept
    has_offset[i] <- cvars[[j]]$has_offset
    has_lfactor[i] <- cvars[[j]]$has_lfactor
    # Note: these will be the same for all yi for categorical and multinomial
    e[[J_fixed]] <- cvars[[j]]$J_fixed
    e[[K_fixed]] <- cvars[[j]]$K_fixed
    e[[J_varying]] <- cvars[[j]]$J_varying
    e[[K_varying]] <- cvars[[j]]$K_varying
    e[[J_random]] <- cvars[[j]]$J_random
    e[[K_random]] <- cvars[[j]]$K_random
    # No phi parameters yet
    # onlyif(!is.null(samples[[phi]]), e[[phi]] <- c(samples[[phi]][idx]))
    onlyif(!is.null(samples[[sigma]]), e$sigma[, i] <- c(samples[[sigma]][idx]))
    # index j here, not from cvars
    if (has_random_effects[j]) {
      nus <- make.unique(rep(paste0("nu_", yi), e[[K_random]]))
      e[[nu]] <- nu_samples[, , nus, drop = FALSE]
    }
    if (has_fixed_intercept[i]) {
      e[[alpha]] <- array(samples[[alpha]][idx], c(e$n_draws, 1L))
    }
    if (has_varying_intercept[i]) {
      e[[alpha]] <- samples[[alpha]][idx, , drop = FALSE]
    }
    if (has_lfactor[i]) {
      e[[lambda]] <- samples[[lambda]][idx, , drop = FALSE]
      e[[psi]] <- samples[[psi]][idx, , drop = FALSE]
    }
    e[[beta]] <- samples[[beta]][idx, , drop = FALSE]
    e[[delta]] <- samples[[delta]][idx, , , drop = FALSE]
    e$xbeta <- matrix(0.0, nrow = e$k, ncol = d)
  }
  e$call <- generate_sim_call_multivariate(
    d = d,
    dims = dims,
    resp = resp,
    resp_levels = resp_levels,
    family = e$family,
    type = type,
    eval_type = eval_type,
    has_fixed = has_fixed,
    has_varying = has_varying,
    has_random = has_random,
    has_fixed_intercept = has_fixed_intercept,
    has_varying_intercept = has_varying_intercept,
    has_random_intercept = has_random_intercept,
    has_offset = has_offset,
    has_lfactor
  )
}

#' Generate an Expression to Evaluate Predictions for a Multivariate Channel
#'
#' @noRd
generate_sim_call_multivariate <- function(d, dims, resp, resp_levels,
                                           family, type, eval_type,
                                           has_fixed, has_varying, has_random,
                                           has_fixed_intercept,
                                           has_varying_intercept,
                                           has_random_intercept,
                                           has_offset, has_lfactor) {
  init_text <- paste0(
    "idx_draw <- seq.int(1L, n_draws) - n_draws\n",
    "for (j in seq_len(n_group)) {\n",
    "  idx_draw <- idx_draw + n_draws\n"
  )
  xbeta_text <- character(d + 1)
  for (i in dims) {
    yi <- ifelse_(
      is_categorical(family),
      paste0(resp, "_", resp_levels[i]),
      resp[i]
    )
    xbeta_text[i] <- glue::glue(
      "\n\n  xbeta[idx_draw, {i}] <- ",
      ifelse_(!has_fixed_intercept[i] && !has_varying_intercept[i], "0", ""),
      ifelse_(has_fixed_intercept[i], "alpha_{yi}", ""),
      ifelse_(has_varying_intercept[i], "alpha_{yi}[, a_time]", ""),
      ifelse_(has_random_intercept[i], "+ nu_{yi}[, j, 1]", ""),
      ifelse_(has_lfactor[i], " + lambda_{yi}[, j] * psi_{yi}[, time]", ""),
      ifelse_(
        has_fixed[i],
        " + .rowSums(
          x = model_matrix[idx_draw, J_fixed_{yi}, drop = FALSE] * beta_{yi},
          m = n_draws,
          n = K_fixed_{yi}
        )",
        ""
      ),
      ifelse_(
        has_varying[i],
        " + .rowSums(
          x = model_matrix[idx_draw, J_varying_{yi}, drop = FALSE] *
                delta_{yi}[, time, ],
          m = n_draws,
          n = K_varying_{yi}
        )",
        ""
      ),
      ifelse_(
        has_random[i],
        " + .rowSums(
          x = model_matrix[idx_draw, J_random_{yi}, drop = FALSE] *
            nu_{yi}[, j, seq.int(1 + {has_random_intercept[i]}, K_random_{yi})],
          m = n_draws,
          n = K_random_{yi} - {has_random_intercept[i]}
        )",
        ""
      )
    )
  }
  xbeta_text[d + 1] <- "}"
  link_text <- character(0L)
  if (identical(type, "link") && identical(eval_type, "predicted")) {
    link_text <- glue::glue(
      "for (i in seq_len(d)) {{\n",
      "  data.table::set(\n",
      "    x = out,\n",
      "    i = idx_data,\n",
      "    j = link_cols[i],\n",
      "    value = xbeta[idx_out, i]\n",
      "  )\n",
      "}}"
    )
  }
  eval_type_text <- glue::glue(predict_expr[[eval_type]][[family$name]])
  type_text <- ifelse_(
    identical(type, "mean") && identical(eval_type, "predicted"),
    glue::glue(predict_expr$mean[[family$name]]),
    ""
  )
  out <- paste(
    c(
      "{",
      init_text,
      xbeta_text,
      link_text,
      eval_type_text,
      type_text,
      "}"
    ),
    collapse = "\n"
  )
  str2lang(out)
}

predict_expr <- list()

# Fitted expressions ------------------------------------------------------

predict_expr$fitted <- list()

predict_expr$fitted$gaussian <- "
  data.table::set(x = out, i = idx, j = '{resp}_fitted', value = xbeta)
"

predict_expr$fitted$mvgaussian <- "
  for (i in seq_len(d)) {{
    data.table::set(
      x = out,
      i = idx,
      j = paste0(resp[i], '_fitted'),
      value = xbeta[, i]
    )
  }}
"

predict_expr$fitted$categorical <- "
  mval <- exp(xbeta - log_sum_exp_rows(xbeta, k, d))
  for (s in 1:d) {{
    data.table::set(
      x = out,
      i = idx,
      j = fitted_cols[s],
      value = mval[, s]
    )
  }}
"

predict_expr$fitted$multinomial <- "
  mval <- exp(xbeta - log_sum_exp_rows(xbeta, k, d))
  for (s in 1:d) {{
    data.table::set(
      x = out,
      i = idx,
      j = fitted_cols[s],
      value = trials * mval[, s]
    )
  }}
"

predict_expr$fitted$bernoulli <- "
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_fitted',
    value = plogis(xbeta)
  )
"

predict_expr$fitted$binomial <- "
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_fitted',
    value = plogis(xbeta)
  )
"

predict_expr$fitted$poisson <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  data.table::set(x = out, i = idx, j = '{resp}_fitted', value = exp_xbeta)
"

predict_expr$fitted$negbin <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  data.table::set(x = out, i = idx, j = '{resp}_fitted', value = exp_xbeta)
"

predict_expr$fitted$exponential <- "
  data.table::set(x = out, i = idx, j = '{resp}_fitted', value = exp(xbeta))
"

predict_expr$fitted$gamma <- "
  data.table::set(x = out, i = idx, j = '{resp}_fitted', value = exp(xbeta))
"

predict_expr$fitted$beta <- "
  data.table::set(x = out, i = idx, j = '{resp}_fitted', value = plogis(xbeta))
"

predict_expr$fitted$student <- "
  data.table::set(x = out, i = idx, j = '{resp}_fitted', value = xbeta)
"

# Predicted expressions ---------------------------------------------------

predict_expr$predicted <- list()

predict_expr$predicted$gaussian <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = rnorm(k_obs, xbeta[idx_out], sigma[idx_out])
  )
"

predict_expr$predicted$mvgaussian <- "
  error <- matrix(0.0, k, d)
  idx_group <- seq.int(1, k, by = n_draws) - 1L
  u <- matrix(0.0, n_group, d)
  for (l in seq_len(n_draws)) {{
    idx_group <- idx_group + 1L
    u[] <- rnorm(n_group * d)
    error[idx_group, ] <- u %*% t(sigma[l, ] %*% L[l, , ])
  }}
  for (i in seq_len(d)) {{
    data.table::set(
      x = out,
      i = idx_data,
      j = resp[i],
      value = xbeta[idx_out, i] + error[idx_out, i]
    )
  }}
"

predict_expr$predicted$categorical <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = max.col(xbeta - log(-log(runif(k * d))))[idx_out]
  )
"

predict_expr$predicted$multinomial <- "
  pred <- matrix(0L, k, d)
  n <- max(trials)
  for (j in seq_len(n)) {{
    rows <- which(j <= trials)
    cols <- max.col(xbeta[rows, ] - log(-log(runif(d * length(rows)))))
    trial_idx <- cbind(rows, cols)
    pred[trial_idx] <- pred[trial_idx] + 1L
  }}
  for (s in seq_len(d)) {{
    data.table::set(
      x = out,
      i = idx_data,
      j = resp_levels[s],
      value = pred[idx_out, s]
    )
  }}
"

predict_expr$predicted$binomial <- "
  prob <- plogis(xbeta[idx_out])
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = rbinom(k_obs, trials[idx_out], prob)
  )
"

predict_expr$predicted$bernoulli <- "
  prob <- plogis(xbeta[idx_out])
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = rbinom(k_obs, 1, prob)
  )
"

predict_expr$predicted$poisson <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = rpois(k_obs, exp_xbeta[idx_out])
  )
"

predict_expr$predicted$negbin <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = rnbinom(k_obs, size = phi[idx_out], mu = exp_xbeta[idx_out])
  )
"

predict_expr$predicted$exponential <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = rexp(k_obs, rate = exp(-xbeta[idx_out]))
  )
"

predict_expr$predicted$gamma <- "
  phi_obs <- phi[idx_out]
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = rgamma(
      k_obs,
      shape = phi_obs,
      rate = phi_obs * exp(-xbeta[idx_out])
    )
  )
"

predict_expr$predicted$beta <- "
  mu <- plogis(xbeta[idx_out])
  phi_obs <- phi[idx_out]
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = rbeta(k_obs,  mu * phi_obs, (1 - mu) * phi_obs)
  )
"

predict_expr$predicted$student <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}',
    value = xbeta[idx_out] + sigma[idx_out] * rt(k_obs, phi[idx_out])
  )
"

# Mean expressions --------------------------------------------------------

predict_expr$mean <- list()

predict_expr$mean$gaussian <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}_mean',
    value = xbeta[idx_out]
  )
"

predict_expr$mean$mvgaussian <- "
  for (i in seq_len(d)) {{
    data.table::set(
      x = out,
      i = idx_data,
      j = paste0(resp[i], '_mean'),
      value = xbeta[idx_out, i]
    )
  }}
"

predict_expr$mean$categorical <- "
  mval <- exp(xbeta - log_sum_exp_rows(xbeta, k, d))
  for (s in 1:d) {{
    data.table::set(
      x = out,
      i = idx_data,
      j = mean_cols[s],
      value = mval[idx_out, s]
    )
  }}
"

predict_expr$mean$multinomial <- "
  mval <- exp(xbeta - log_sum_exp_rows(xbeta, k, d))
  for (s in 1:d) {{
    data.table::set(
      x = out,
      i = idx_data,
      j = mean_cols[s],
      value = trials[idx_out] * mval[idx_out, s]
    )
  }}
"

predict_expr$mean$binomial <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}_mean',
    value = prob[idx_out]
  )
"

predict_expr$mean$bernoulli <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}_mean',
    value = prob[idx_out]
  )
"

predict_expr$mean$poisson <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}_mean',
    value = exp_xbeta[idx_out]
  )
"

predict_expr$mean$negbin <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}_mean',
    value = exp_xbeta[idx_out]
  )
"

predict_expr$mean$exponential <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}_mean',
    value = exp(xbeta[idx_out])
  )
"

predict_expr$mean$gamma <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}_mean',
    value = exp(xbeta[idx_out])
  )
"

predict_expr$mean$beta <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}_mean',
    value = mu[idx_out]
  )
"

predict_expr$mean$student <- "
  data.table::set(
    x = out,
    i = idx_data,
    j = '{resp}_mean',
    value = xbeta[idx_out]
  )
"

# Log-likelihood expressions ----------------------------------------------

predict_expr$loglik <- list()

predict_expr$loglik$gaussian <- "
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = dnorm(y, xbeta, sigma, log = TRUE)
  )
"

predict_expr$loglik$mvgaussian <- "
  ll <- numeric(k)
  for (l in seq_len(n_draws)) {{
    idx_group <- seq.int(l, k, by = n_draws)
    sigma_chol <- diag(sigma[l, ]) %*% L[l, , ]
    log_det <- 2.0 * sum(log(diag(sigma_chol)))
    diffs <- t(y[idx_group, ] - xbeta[idx_group, ])
    z <- forwardsolve(sigma_chol, diffs)
    quad <- colSums(z^2)
    ll[idx_group] <- -0.5 * (d * log(2 * pi) + log_det + quad)
  }}
  for (i in seq_len(d)) {{
    data.table::set(
      x = out,
      i = idx,
      j = paste(c(resp, 'loglik'), collapse = '_'),
      value = ll
    )
  }}
"

predict_expr$loglik$categorical <- "
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = xbeta[cbind(seq_along(y), y)] - log_sum_exp_rows(xbeta, k, d)
  )
"

predict_expr$loglik$multinomial <- "
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = lgamma(trials + 1) +
      .rowSums(
        y * (xbeta - log_sum_exp_rows(xbeta, k, S)) - lgamma(y + 1), k, d
      )
  )
"

predict_expr$loglik$binomial <- "
  prob <- plogis(xbeta)
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = dbinom(y, trials, prob, log = TRUE)
  )
"

predict_expr$loglik$bernoulli <- "
  prob <- plogis(xbeta)
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = dbinom(y, 1, prob, log = TRUE)
  )
"

predict_expr$loglik$poisson <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = dpois(y, exp_xbeta, log = TRUE)
  )
"

predict_expr$loglik$negbin <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = dnbinom(y, size = phi, mu = exp_xbeta, log = TRUE)
  )
"

predict_expr$loglik$exponential <- "
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = dexp(y, rate = exp(-xbeta), log = TRUE)
  )
"

predict_expr$loglik$gamma <- "
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = dgamma(y, shape = phi, rate = phi * exp(-xbeta), log = TRUE)
  )
"

predict_expr$loglik$beta <- "
  mu <- plogis(xbeta)
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = dbeta(y,  mu * phi, (1 - mu) * phi, log = TRUE)
  )
"

predict_expr$loglik$student <- "
  data.table::set(
    x = out,
    i = idx,
    j = '{resp}_loglik',
    value = dt((y - xbeta)/sigma, phi, log = TRUE)
  )
"

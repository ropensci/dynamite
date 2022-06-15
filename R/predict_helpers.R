# TODO documentation
check_newdata <- function(newdata, data, type, families_stoch, resp_stoch,
                          categories, group_var, time_var) {
  if (!is.null(group_var)) {
    if (!(group_var %in% names(newdata))) {
      stop_("Can't find grouping variable {.var {group_var}} in {.var newdata}.")
    }
    group <- newdata[[group_var]]
    if (is.factor(group)) {
      # TODO is this necessary? only length of unique values matters
      group <- droplevels(group)
    }
    group <- unique(group)
    # TODO doesn't really matter at least at the moment
    if (!all(group %in% data[[group_var]])) {
      new_levels <- unique(group[!group %in% data[[group_var]]])
      stop_(c(
        "Grouping variable {.var {group_var}} contains unknown levels:",
        `x` = "Level{?s} {.val {as.character(new_levels)}}
               {?is/are} not present in the original data."
      ))
    }
  }
  if (!(time_var %in% names(newdata))) {
    stop_("Can't find time index variable {.var {time_var}} in {.var newdata}.")
  }
  time <- unique(newdata[[time_var]])
  if (!all(time %in% data[[time_var]])) {
    new_times <- unique(time[!time %in% data[[time_var]]])
    stop_(c(
      "Time index variable {.var {time_var}} contains unknown time points:",
      `x` = "Time point{?s} {.val {as.character(new_times)}}
             {?is/are} not present in the original data."
    ))
  }
  for (resp in resp_stoch) {
    if (is.null(newdata[[resp]])) {
      stop_("Can't find response variable {.var {resp}} in {.var newdata}.")
    }
  }
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
    if (identical(type, "response")) {
      newdata[, (glue::glue("{resp}_new")) := newdata[[resp]]]
    }
    if (!identical(type, "fitted")) {
      newdata[, (glue::glue("{resp}_store")) := newdata[[resp]]]
    }
  }
  data.table::setkeyv(newdata, c(group_var, time_var))
}

# TODO test imputation
impute_newdata <- function(newdata, impute, predictors, group_var) {
  if (identical(impute, "locf")) {
    newdata[, (predictors) := lapply(.SD, locf),
            .BY = group_var, .SDcols = predictors]
  }
}

clear_nonfixed <- function(newdata, newdata_null, resp_stoch,
                           group_var, clear_names, fixed, n_id, n_time) {
  # TODO maybe not needed
  # non_na <- newdata |>
  #   dplyr::group_by(.data[[group_var]]) |>
  #   dplyr::summarise(
  #     obs = stats::complete.cases(
  #       dplyr::across(
  #         dplyr::all_of(resp_stoch)
  #       )
  #     ), .groups = "keep")
  # fixed_obs <- non_na |>
  #   dplyr::summarise(
  #     first_obs = which(.data$obs)[1],
  #     horizon = all(.data$obs[.data$first_obs:(.data$first_obs + n_fixed - 1)])
  #   )
  # lacking_obs <- is.na(fixed_obs$horizon) | (fixed_obs$horizon < n_fixed)
  # if (any(lacking_obs)) {
  #   groups_lacking <- unique(fixed_obs[lacking_obs, group_var])
  #   stop_("Insufficient non-NA observations in groups: ", cs(groups_lacking))
  # }
  if (newdata_null) {
    # predict_idx <- unlist(lapply(seq_len(n_id), function(i) {
    #   (fixed_obs$first_obs[i] + n_fixed):n_time + (i - 1) * n_time
    # }))
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

prepare_eval_envs <- function(object, newdata, eval_type, predict_type,
                              resp_stoch, n_id, n_draws) {
  samples <- rstan::extract(object$stanfit)
  model_vars <- object$stan$model_vars
  specials <- evaluate_specials(object$dformulas$stoch, newdata)
  newdata_names <- names(newdata)
  n_resp <- length(resp_stoch)
  eval_envs <- vector(mode = "list", length = n_resp)
  for (j in seq_len(n_resp)) {
    resp <- resp_stoch[j]
    resp_family <- object$dformulas$stoch[[j]]$family
    alpha <- paste0("alpha_", resp)
    beta <- paste0("beta_", resp)
    delta <- paste0("delta_", resp)
    phi <- paste0("phi_", resp)
    sigma <- paste0("sigma_", resp)
    nu <- paste0("nu_", resp)
    e <- new.env()
    e$type <- predict_type
    e$n_draws <- n_draws
    e$k <- n_draws * n_id
    e$newdata <- newdata
    e$J_fixed <- model_vars[[j]]$J_fixed
    e$K_fixed <- model_vars[[j]]$K_fixed
    e$J_varying <- model_vars[[j]]$J_varying
    e$K_varying <- model_vars[[j]]$K_varying
    e$resp <- resp_stoch[j]
    e$phi <- samples[[phi]]
    e$sigma <- samples[[sigma]]
    e$offset <- specials[[j]]$offset
    e$trials <- specials[[j]]$trials
    if (is_categorical(resp_family)) {
      resp_levels <- attr(object$stan$responses, "resp_class")[[resp]] |>
        attr("levels")
      e$alpha <- samples[[alpha]]
      e$beta <- samples[[beta]]
      e$delta <- samples[[delta]]
      e$nu <- samples[[nu]]
      e$resp_levels <- resp_levels
      e$S <- length(resp_levels)
    } else {
      resp_levels <- NULL
      e$alpha <- samples[[alpha]]
      e$beta <- samples[[beta]]
      e$delta <- samples[[delta]]
      if (model_vars[[j]]$has_random_intercept) {
        e$nu <- c(t(samples[[nu]][1L:n_draws, ]))
      }
    }
    e$call <- generate_sim_call(resp, resp_levels, resp_family, eval_type,
                                model_vars[[j]]$has_fixed,
                                model_vars[[j]]$has_varying,
                                model_vars[[j]]$has_fixed_intercept,
                                model_vars[[j]]$has_varying_intercept,
                                model_vars[[j]]$has_random_intercept,
                                model_vars[[j]]$has_offset)
    eval_envs[[j]] <- e
  }
  eval_envs
}

generate_sim_call <- function(resp, resp_levels, family, type,
                              has_fixed, has_varying,
                              has_fixed_intercept, has_varying_intercept,
                              has_random_intercept, has_offset) {
  if (is_categorical(family)) {
    glue::glue(
      "{{\n",
      paste0(
        "  xbeta <- matrix(0.0, k, S)\n",
        "  for (j in seq_len(n_draws)) {{\n",
        "{ifelse_(has_fixed,
                  'X_f <- model_matrix[j, J_fixed, drop = FALSE]', '')}\n",
        "{ifelse_(has_varying,
                  'X_v <- model_matrix[j, J_varying, drop = FALSE]', '')}\n",
        "    idx_draw <- seq.int((j - 1L) * n_id + 1L, j * n_id)\n",
        "    for (s in seq_len(S - 1)) {{\n",
        "      xbeta[idx_draw, s + 1] <- ",
        "{ifelse_(!has_fixed_intercept && !has_varying_intercept, '0', '')}",
        "{ifelse_(has_fixed_intercept, 'alpha[j, s]', '')}",
        "{ifelse_(has_varying_intercept, 'alpha[j, a_time, s]', '')}",
        "{ifelse_(has_fixed, ' + X_f %*% beta[j, , s]', '')}",
        "{ifelse_(has_varying, ' + X_v %*% delta[j, time, , s]', '')}",
        "    }}",
        "  }}\n"
      ),
      eval(str2lang(glue::glue("{type}_categorical"))),
      "}}"
    ) |> str2lang()
  } else {
    glue::glue(
      "{{\n",
      paste0(
        "  xbeta <- numeric(k)\n",
        "  for (j in seq_len(n_draws)) {{\n",
        "{ifelse_(has_fixed,
                  'X_f <- model_matrix[j, J_fixed, drop = FALSE]', '')}\n",
        "{ifelse_(has_varying,
                  'X_v <- model_matrix[j, J_varying, drop = FALSE]', '')}\n",
        "    idx_draw <- seq.int((j - 1L) * n_id + 1L, j * n_id)\n",
        "    xbeta[idx_draw] <- ",
        "{ifelse_(!has_fixed_intercept && !has_varying_intercept, '0', '')}",
        "{ifelse_(has_fixed_intercept, 'alpha[j,]', '')}",
        "{ifelse_(has_varying_intercept, 'alpha[j, a_time]', '')}",
        "{ifelse_(has_fixed, ' + X_f %*% beta[j, ]', '')}",
        "{ifelse_(has_varying, ' + X_v %*% delta[j, time, ]', '')}",
        "  }}\n"
      ),
      ifelse_(identical(type, "predict"),
        paste0(
          "  if (type == 'link') {{",
          "    data.table::set(x = newdata, i = idx_pred, j = '{resp}_link',
                               value = xbeta)",
          "  }}"
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
  data.table::set(x = newdata, i = idx, j = '{resp}_fitted',
                  value = plogis(xbeta))
"

fitted_binomial <- "
  data.table::set(x = newdata, i = idx, j = '{resp}_fitted',
                  value = plogis(xbeta))
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
# Predict expressions -----------------------------------------------------

predict_gaussian <- "
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_mean',
                    value = xbeta)
  }}
  data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                  value = rnorm(k, xbeta, sigma))
"

predict_categorical <- "
    if (type == 'link') {{
      resp_cols <- c({paste0('\"', resp, '_link_', resp_levels, '\"',
                             collapse = ', ')})
      for (s in 1:S) {{
        data.table::set(
          x = newdata,
          i = idx_pred,
          j = resp_cols[s],
          value = xbeta[, s]
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
          i = idx_pred,
          j = resp_cols[s],
          value = mval[, s]
        )
      }}
    }}
    data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                    value = max.col(xbeta - log(-log(runif(S * k)))))
"

predict_binomial <- "
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_mean',
                    value = plogis(xbeta))
  }}
  data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                  value = rbinom(k, trials, plogis(xbeta)))
"

predict_bernoulli <- "
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_mean',
                    value = plogis(xbeta))
  }}
  data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                  value = rbinom(k, 1, plogis(xbeta)))
"

predict_poisson <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_mean',
                    value = exp_xbeta)
  }}
  data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                  value = rpois(k, exp_xbeta))
"

predict_negbin <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_mean',
                    value = exp_xbeta)
  }
  data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                  value = rnbinom(k, size = phi, mu = exp_xbeta))
"

predict_exponential <- "
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_mean',
                    value = exp(xbeta))
  }
  data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                  value = rexp(k, rate = exp(-xbeta)))
"

predict_gamma <- "
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_mean',
                    value = exp(xbeta))
  }
  data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                  value = rgamma(k, shape = phi, rate = phi * exp(-xbeta)))
"

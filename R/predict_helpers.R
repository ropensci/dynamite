# TODO documentation
check_newdata <- function(newdata, data, group_var, time_var) {
  if (is.null(newdata)) {
    newdata <- data
  } else {
    if (!(group_var %in% names(newdata))) {
      stop_("Grouping variable '", group_var, "' not found in 'newdata'")
    }
    group <- newdata[[group_var]]
    if (is.factor(group)) {
      # TODO is this necessary? only length of unique values matters
      group <- droplevels(group)
    }
    group <- unique(group)
    # TODO doesn't really matter at least at the moment
    if (!all(group %in% data[[group_var]])) {
      stop_("Grouping variable '", group_var, "' ",
            "contains new levels not found in the original data")
    }
    if (!(time_var %in% names(newdata))) {
      stop_("Time index variable '", time_var, "' not found in 'newdata'")
    }
    time <- unique(newdata[[time_var]])
    if (!all(time %in% data[[time_var]])) {
      stop_("Timing variable '", time_var, "' ",
            "contains time points not found in the original data")
    }
  }
  newdata
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

softmax <- function(x) {
  exp(x - log_sum_exp(x))
}

# Fitted expressions ------------------------------------------------------

fitted_gaussian <- quote({
  c(model_matrix %*% matrix(samples, nrow = ncol(model_matrix), byrow = TRUE))
})

fitted_categorical <- quote({
  n_draws <- nrow(samples)
  n_id <- nrow(model_matrix)
  sim <- matrix(NA, n_draws * n_id, dim(samples)[4] + 1)
  for (k in seq_len(n_draws)) {
    idx_k <- ((k - 1) * n_id + 1):(k * n_id)
    xbeta <- cbind(model_matrix %*% samples[k, , , ], 0)
    maxs <- apply(xbeta, 1, max)
    sim[idx_k, ] <- exp(xbeta - (maxs + log(rowSums(exp(xbeta - maxs)))))
  }
})

fitted_bernoulli <- quote({
  plogis(c(model_matrix %*% t(samples)))
})

fitted_binomial <- quote({
  fitted_bernoulli(model_matrix, samples)
})

fitted_poisson <- quote({
  exp(c(model_matrix %*% t(samples)))
})

fitted_negbin <- quote({
  exp(c(model_matrix %*% t(samples)))
})

# Predict expressions -----------------------------------------------------

generate_predict_call <- function(resp, family, has_fixed, has_varying,
                                  has_fixed_intercept, has_varying_intercept,
                                  has_random_intercept, has_offset) {
  if (is_categorical(family)) {
    glue::glue(predict_categorical) |> str2lang()
  } else {
    glue::glue(
      "{{\n",
      paste0(
        "{ifelse_(!has_fixed_intercept && !has_varying_intercept, '  0', '')}",
        "{ifelse_(has_fixed_intercept, '  xbeta <- alpha', '')}",
        "{ifelse_(has_varying_intercept, '  xbeta <- alpha[, a_time]', '')}",
        "{ifelse_(has_random_intercept, ' + nu', '')}",
        "{ifelse_(has_fixed,",
          "' + .rowSums(x = model_matrix[, J_fixed, drop = FALSE] * beta, m = k, n = J_fixed)', '')}",
          "{ifelse_(has_varying,",
        "' + .rowSums(x = model_matrix[, J_varying, drop = FALSE] * delta[, time, ], m = k, n = J_varying)', '')}\n"
      ),
      paste0(
      "  if (type == 'link') {{",
      "    data.table::set(x = newdata, i = idx_pred, j = '{resp}_store', value = xbeta)",
      "  }}"
      ),
      eval(str2lang(paste0("predict_", family))),
      "}}"
    ) |> str2lang()
  }
}

predict_gaussian <- "
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_store', value = xbeta)
  }}
  data.table::set(x = newdata, i = idx_pred, j = '{resp}', value = rnorm(k, xbeta, sigma))
"

predict_categorical <- "
  xbeta <- alpha[idx_par, a_time, , drop = FALSE] {ifelse_(has_fixed || has_varying, '+', '')}
    cbind(
      {ifelse_(has_fixed,
         rowSums(model_matrix[, J_fixed, drop = FALSE] * beta)', '')}
      {ifelse_(has_fixed && has_varying), '+', ''}
      {ifelse_(has_varying,
         rowSums(model_matrix[, J_varying, drop = FALSE] * delta[, time, , , drop = FALSE])', '')}
      {ifelse_(has_fixed || has_varying, ', 0', '0')}
    )
  maxs <- apply(xbeta, 1, max)
  if (type == 'link') {{
    data.table::set(
      x = newdata,
      i = idx_pred,
      j = '{resp}_{resp_levels}',
      value = xbeta)
  }}
  if (type == 'mean') {{
    data.table::set(
      x = newdata,
      i = idx_pred,
      j = '{resp}_{resp_levels}',
      value = exp(xbeta - (maxs + log(rowSums(exp(xbeta - maxs)))))
    )
  }}
  data.table::set(x = newdata, i = idx_pred, j = '{resp}', value = max.col(xbeta - log(-log(runif(S * k)))))
"

predict_binomial <- "
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_store', value = plogis(xbeta))
  }}
  data.table::set(x = newdata, i = idx_pred, j = '{resp}', value = rbinom(k, trials, plogis(xbeta)))
"

predict_bernoulli <- "
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_store', value = plogis(xbeta))
  }}
  data.table::set(x = newdata, i = idx_pred, j = '{resp}', value = rbinom(k, 1, plogis(xbeta)))
"

predict_poisson <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_store',
                    value = exp_xbeta)
  }}
  data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                  value = rpois(k, exp_xbeta))
"

predict_negbin <- "
  exp_xbeta <- {ifelse_(has_offset, 'exp(xbeta + offset)', 'exp(xbeta)')}
  if (type == 'mean') {{
    data.table::set(x = newdata, i = idx_pred, j = '{resp}_store', value = exp_xbeta)
  }
  data.table::set(x = newdata, i = idx_pred, j = '{resp}',
                  value = rnbinom(n_id, size = phi, mu = exp_xbeta))
"

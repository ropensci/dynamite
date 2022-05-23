# TODO documentation
check_newdata <- function(newdata, data, group_var, time_var) {
  if (is.null(newdata)) {
    newdata <- data |>
      dplyr::arrange(
        dplyr::across(
          dplyr::all_of(c(group_var, time_var))
        )
      )
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
    newdata <- newdata |>
      dplyr::arrange(
        dplyr::across(
          dplyr::all_of(c(group_var, time_var))
        )
      )
    # TODO just use the original time points starting from start_time
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

fitted_gaussian <- function(model_matrix, samples) {
  c(model_matrix %*% matrix(samples, nrow = ncol(model_matrix), byrow = TRUE))
}

fitted_categorical <- function(model_matrix, samples) {
  n_draws <- nrow(samples)
  n_id <- nrow(model_matrix)
  sim <- matrix(NA, n_draws * n_id, dim(samples)[4] + 1)
  for (k in seq_len(n_draws)) {
    idx_k <- ((k - 1) * n_id + 1):(k * n_id)
    xbeta <- cbind(model_matrix %*% samples[k, , , ], 0)
    maxs <- apply(xbeta, 1, max)
    sim[idx_k, ] <- exp(xbeta - (maxs + log(rowSums(exp(xbeta - maxs)))))
  }
  sim
}

fitted_bernoulli <- function(model_matrix, samples) {
  plogis(c(model_matrix %*% t(samples)))
}

fitted_binomial <- function(model_matrix, samples) {
  fitted_bernoulli(model_matrix, samples)
}

fitted_poisson <- function(model_matrix, samples) {
  exp(c(model_matrix %*% t(samples)))
}

fitted_negbin <- function(model_matrix, samples) {
  exp(c(model_matrix %*% t(samples)))
}

predict_gaussian <- function(model_matrix, samples, specials,
                             resp, time, type, n_draws, J_fixed, J_varying) {
  beta <- samples[[paste0("beta_", resp)]]
  delta <- samples[[paste0("delta_", resp)]]
  n_id <- nrow(model_matrix) / n_draws
  sim <- sim_r <- numeric(nrow(model_matrix))
  for (k in seq_len(n_draws)) {
    idx_k <- ((k - 1) * n_id + 1):(k * n_id)
    xbeta <- model_matrix[idx_k, J_fixed, drop = FALSE] %*% beta[k, ] +
      model_matrix[idx_k, J_varying, drop = FALSE] %*% delta[k, time, ]
    sim[idx_k] <- xbeta
    sigma <- samples[[paste0("sigma_", resp)]][k]
    sim_r[idx_k] <- rnorm(n_id, xbeta, sigma)
  }
  list(response = sim_r, mean_or_link = sim)
}

predict_categorical <- function(model_matrix, samples, specials,
                                resp, time, type, n_draws) {
  beta <- samples[[paste0("beta_", resp)]]
  n_id <- nrow(model_matrix) / n_draws
  S <- dim(samples[[paste0("beta_", resp)]])[4] + 1
  sim <- matrix(NA, nrow(model_matrix), S)
  sim_r <- numeric(nrow(model_matrix))
  for (k in seq_len(n_draws)) {
    idx_k <- ((k - 1) * n_id + 1):(k * n_id)
    xbeta <- cbind(model_matrix[idx_k, , drop = FALSE] %*% beta[k, time, , ], 0)
    maxs <- apply(xbeta, 1, max)
    if (type == "link") {
      sim[idx_k, ] <- xbeta
    }
    if (type == "mean") {
      sim[idx_k, ] <- exp(xbeta - (maxs + log(rowSums(exp(xbeta - maxs)))))
    }
    sim_r[idx_k] <- max.col(xbeta - log(-log(runif(S * n_id))))
  }
  list(response = sim_r, mean_or_link = sim)
}

predict_bernoulli <- function(model_matrix, samples, specials,
                              resp, time, type, n_draws) {
  predict_binomial(model_matrix, samples, specials, resp, time, type, n_draws)
}

predict_binomial <- function(model_matrix, samples, specials,
                             resp, time, type, n_draws) {
  beta <- samples[[paste0("beta_", resp)]]
  n_id <- nrow(model_matrix) / n_draws
  sim <- sim_r <- numeric(nrow(model_matrix))
  trials <- specials$trials
  if (is.null(trials)) {
    trials <- rep(1, n_id)
  }
  for (k in seq_len(n_draws)) {
    idx_k <- ((k - 1) * n_id + 1):(k * n_id)
    xbeta <- model_matrix[idx_k, , drop = FALSE] %*% beta[k, time, ]
    if (type == "link") {
      sim[idx_k, ] <- xbeta
    }
    if (type == "mean") {
      sim[idx_k] <- plogis(xbeta)
    }
    sim_r[idx_k] <- rbinom(n_id, trials, plogis(xbeta))
  }
  list(response = sim_r, mean_or_link = sim)
}

predict_poisson <- function(model_matrix, samples, specials,
                            resp, time, type, n_draws) {
  beta <- samples[[paste0("beta_", resp)]]
  n_id <- nrow(model_matrix) / n_draws
  sim <- sim_r <- numeric(nrow(model_matrix))
  offset <- specials$offset
  has_offset <- is.null(offset)
  for (k in seq_len(n_draws)) {
    idx_k <- ((k - 1) * n_id + 1):(k * n_id)
    xbeta <- model_matrix[idx_k, , drop = FALSE] %*% beta[k, time, ]
    exp_xbeta <- if (has_offset) exp(xbeta + offset) else exp(xbeta)
    if (type == "link") {
      sim[idx_k, ] <- xbeta
    }
    if (type == "mean") {
      sim[idx_k] <- exp_xbeta
    }
    sim_r[idx_k] <- rpois(n_id, exp_xbeta)
  }
  list(response = sim_r, mean_or_link = sim)
}

predict_negbin <- function(model_matrix, samples, specials,
                           resp, time, type, n_draws) {
  beta <- samples[[paste0("beta_", resp)]]
  phi <- samples[[paste0("phi_", resp)]]
  n_id <- nrow(model_matrix) / n_draws
  sim <- sim_r <- numeric(nrow(model_matrix))
  for (k in seq_len(n_draws)) {
    idx_k <- ((k - 1) * n_id + 1):(k * n_id)
    xbeta <- model_matrix[idx_k, , drop = FALSE] %*% beta[k, time, ]
    exp_xbeta <- exp(xbeta)
    if (type == "link") {
      sim[idx_k, ] <- xbeta
    }
    if (type == "mean") {
      sim[idx_k] <- exp_xbeta
    }

    sim_r[idx_k] <- rnbinom(n_id, size = phi[k], mu = exp_xbeta)
  }
  list(response = sim_r, mean_or_link = sim)
}

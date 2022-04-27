log_sum_exp <- function (x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
}

softmax <- function (x) {
    exp(x - log_sum_exp(x))
}

fitted_gaussian <- function(model_matrix, samples) {
    c(model_matrix %*% matrix(samples, nrow = ncol(model_matrix), byrow = TRUE))
}


fitted_categorical <- function(model_matrix, samples) {
    n_draws <- nrow(samples)
    n_id <- nrow(model_matrix)
    sim <- matrix(NA, n_draws * n_id, dim(samples)[4] + 1)
    for(k in seq_len(n_draws)) {
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

predict_gaussian <- function(model_matrix, samples, specials, resp, time, type, n_draws) {
    beta <- samples[[paste0("beta_", resp)]]
    n_id <- nrow(model_matrix) / n_draws
    sim <- sim_r <- numeric(nrow(model_matrix))
    for(k in seq_len(n_draws)) {
        idx_k <- ((k - 1) * n_id + 1):(k * n_id)
        xbeta <- model_matrix[idx_k, , drop = FALSE] %*% beta[k, time, ]
        sim[idx_k] <- xbeta
        sigma <- samples[[paste0("sigma_", resp)]][k]
        sim_r[idx_k] <- rnorm(n_id, xbeta, sigma)
    }
    list(response = sim_r, mean_or_link = sim)
}

predict_categorical <- function(model_matrix, samples, specials, resp, time, type, n_draws) {
    beta <- samples[[paste0("beta_", resp)]]
    n_id <- nrow(model_matrix) / n_draws
    S <- dim(samples[[paste0("beta_", resp)]])[4] + 1
    sim <- matrix(NA, nrow(model_matrix), S)
    sim_r <- numeric(nrow(model_matrix))
    for(k in seq_len(n_draws)) {
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

predict_bernoulli <- function(model_matrix, samples, specials, resp, time, type, n_draws) {
    predict_binomial(model_matrix, samples, specials, resp, time, type, n_draws)
}

predict_binomial <- function(model_matrix, samples, specials, resp, time, type, n_draws) {
    beta <- samples[[paste0("beta_", resp)]]
    n_id <- nrow(model_matrix) / n_draws
    sim <- sim_r <- numeric(nrow(model_matrix))
    trials <- specials$trials
    if (is.null(trials)) {
        trials <- rep(1, n_id)
    }
    for(k in seq_len(n_draws)) {
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

predict_poisson <- function(model_matrix, samples, specials, resp, time, type, n_draws) {
    beta <- samples[[paste0("beta_", resp)]]
    n_id <- nrow(model_matrix) / n_draws
    sim <- sim_r <- numeric(nrow(model_matrix))
    offset <- specials$offset
    has_offset <- is.null(offset)
    for(k in seq_len(n_draws)) {
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

predict_negbin <- function(model_matrix, samples, specials, resp, time, type, n_draws) {
    beta <- samples[[paste0("beta_", resp)]]
    phi <- samples[[paste0("phi_", resp)]]
    n_id <- nrow(model_matrix) / n_draws
    sim <- sim_r <- numeric(nrow(model_matrix))
    for(k in seq_len(n_draws)) {
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

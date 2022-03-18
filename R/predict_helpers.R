log_sum_exp <- function (x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
}

softmax <- function (x) {
    exp(x - log_sum_exp(x))
}

fitted_gaussian <- function(model_matrix, samples) {
    c(model_matrix %*% t(samples))
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
predict_gaussian <- function(model_matrix, samples, resp, time, type) {
    beta <- samples[[paste0("beta_", resp)]]
    n_draws <- nrow(beta)
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
predict_categorical <- function(model_matrix, samples, resp, time, type) {
    beta <- samples[[paste0("beta_", resp)]]
    n_draws <- nrow(beta)
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

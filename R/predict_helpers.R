log_sum_exp <- function (x) {
    m <- max(x)
    m + log(sum(exp(x - m)))
}

softmax <- function (x) {
    exp(x - log_sum_exp(x))
}

# categorical_rng <- function(y, x, beta, J) {
#     sample(y, 1, prob = softmax(x[J] %*% beta))
# }
# gaussian_rng <- function(x, beta, sigma, J) {
#     rnorm(1, x[J] %*% beta, sigma)
# }
# create_predict_functions <- function(x) {
#     resp_all <- get_resp(x)
#     rng <- vector("list", length(resp_all))
#     names(rng) <- resp_all
#     families <- lapply(x, function(x) x$family$name)
#     for (i in seq_along(resp_all)) {
#         rng[[i]] <- get(paste0(families[i], "_rng"))
#     }
#     rng
# }


# categorical_mean <- function(y, x, beta, J) {
#     softmax(x[J] %*% beta)
# }
# gaussian_mean <- function(x, beta, J) {
#     x[J] %*% beta
# }
# create_mean_functions <- function(x) {
#     resp_all <- get_resp(x)
#     rng <- vector("list", length(resp_all))
#     names(rng) <- resp_all
#     families <- lapply(x, function(x) x$family$name)
#     for (i in seq_along(resp_all)) {
#         rng[[i]] <- get(paste0(families[i], "_mean"))
#     }
#     rng
# }


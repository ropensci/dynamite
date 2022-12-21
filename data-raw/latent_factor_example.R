# code to create `latent_factor_example` object

set.seed(123)
N <- 40L
T_ <- 20L
D <- 10
B <- t(splines::bs(1:T_, df = D, intercept = TRUE))
a <- cumsum(rnorm(D))
psi <- numeric(T_)
lambda_i <- rnorm(N, 1, 0.2)
alpha_i <- rnorm(N, 0, 0.5)
for(t in 1:T_) {
  psi[t] <- B[, t] %*% a
}
y <- matrix(0, N, T_)
for(t in 1:T_) {
  y[, t] <- rnorm(N, alpha_i + lambda_i * psi[t], 0.2)
}
latent_factor_example <- data.frame(
  y = c(y),
  id = seq_len(N),
  time = rep(seq_len(T_), each = N)
)
usethis::use_data(latent_factor_example, overwrite = TRUE)

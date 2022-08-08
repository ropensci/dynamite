# code to create `gaussian_example` object

set.seed(123)
N <- 50L
T_ <- 30L
K_fixed <- 1L
K_varying <- 3L
K <- 4L
tau <- c(0.2, 0.4, 0.1)
sigma <- 0.2
beta <- 2.0
Bs <- t(splines::bs(seq.int(2L, T_), df = 20L, degree = 3L, intercept = TRUE))
D <- nrow(Bs)
sigma_nu <- 0.1
nu <- rnorm(N, 0.0, sigma_nu)
a <- array(0.0, c(K_varying, D))
# Splines start from t = 2, first time point is fixed
delta <- array(NA, c(T_, K_varying))

for (k in seq_len(K_varying)) {
  a[k, ] <- cumsum(rnorm(D, 0, tau[k]))
  for (t in seq.int(2L, T_)) {
    delta[t, k] <- a[k, ] %*% Bs[, t - 1]
  }
}
x <- matrix(rnorm(T_ * N), N, T_)
z <- matrix(rbinom(T_ * N, 1.0, 0.7), N, T_)
y <- matrix(NA, N, T_)
y[, 1L] <- rnorm(N)
for (t in seq.int(2L, T_)) {
  m <- nu + beta * z[, t] + delta[t, 1L] +
    delta[t, 2L] * x[, t] + delta[t, 3L] * y[, t - 1L]
  y[, t] <- rnorm(N, m, sigma)
}

gaussian_example <- data.frame(
  y = c(y), x = c(x), z = c(z), id = seq_len(N),
  time = rep(seq_len(T_), each = N))

usethis::use_data(gaussian_example, overwrite = TRUE)

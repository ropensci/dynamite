## code to create `gaussian_example` dataset

set.seed(123)
N <- 50
T <- 30
K_fixed <- 1
K_varying <- 3
K <- 4
tau <- c(0.2, 0.4, 0.1)
sigma <- 0.2
beta <- 2
Bs <- t(splines::bs(2:T, df = 20, degree = 3, intercept = TRUE))
D <- nrow(Bs)
sigma_nu <- 0.1
nu <- rnorm(N, 0, sigma_nu)

a <- array(0, c(K_varying, D))
# Splines start from t=2, first time point is fixed
delta <- array(NA, c(T, K_varying))

for(k in 1:K_varying) {
  a[k, ] <- cumsum(rnorm(D, 0, tau[k]))
  for(t in 2:T){
    delta[t, k] = a[k, ] %*% Bs[, t - 1]
  }
}
x <- matrix(rnorm(T * N), N, T)
z <- matrix(rbinom(T * N, 1, 0.7), N, T)
y <- matrix(NA, N, T)
y[, 1] <- rnorm(N)
for(t in 2:T) {
  m <- nu + beta * z[, t] + delta[t, 1] + delta[t, 2] * x[, t] + delta[t, 3] * y[, t - 1]
  y[, t] <- rnorm(N, m, sigma)
}

gaussian_example <- data.frame(
  y = c(y), x = c(x), z = c(z), id = 1:N,
  time = rep(1:T, each = N))
usethis::use_data(gaussian_example, overwrite = TRUE)

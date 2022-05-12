## code to create `gaussian_example` dataset

set.seed(123)
N <- 100
T <- 30
K_fixed <- 1
K_varying <- 3
K <- 4
tau <- c(0.2, 0.4, 0.1)
sigma <- 0.2
beta <- 2
Bs <- t(splines::bs(1:T, df = 20, degree = 3, intercept = TRUE))
D <- nrow(Bs)

a <- array(0, c(K_varying, D))
delta <- array(0, c(T, K_varying))

for(k in 1:K_varying) {
  a[k, ] <- cumsum(c(k/4, rnorm(D - 1, 0, tau[k])))
  for(t in 1:T){
    delta[t, k] = a[k, ] %*% Bs[, t]
  }
}

x <- matrix(rnorm(T * N), N, T)
z <- matrix(rbinom(T * N, 1, 0.7), N, T)
y <- matrix(NA, N, T)
y[, 1] <- beta * z[, 1] + delta[1, 1] + delta[1, 2] * x[, 1] + rnorm(N)
for(t in 2:T) {
  m <- beta * z[, t] + delta[t, 1] + delta[t, 2] * x[, t] + delta[t, 3] * y[, t - 1]
  y[, t] <- rnorm(N, m, sigma)
}

gaussian_example <- data.frame(y = c(y), x = c(x), z = c(z), id = 1:N,
  time = rep(1:T, each = N))
usethis::use_data(gaussian_example, overwrite = TRUE)

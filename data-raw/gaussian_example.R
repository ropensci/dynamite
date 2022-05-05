## code to prepare `gaussian_example` dataset

set.seed(123)
N <- 100
T <- 30
K <- 3
tau <- c(0.2, 0.4, 0.1)
sigma <- 0.2

X <- array(0, c(T, N, K))
X[, , 1] <- 1
X[, , 2] <- rnorm(T * N)
Bs <- t(splines::bs(1:T, df = 20, degree = 3, intercept = TRUE))
D <- nrow(Bs)

a <- array(0, c(K, D))

beta <- array(0, c(T, K))
a <- array(0, c(K, D))

for(k in 1:K) {
  a[k, ] <- cumsum(c(k/4, rnorm(D - 1, 0, tau[k])))
  for(t in 1:T){
    beta[t, k] = a[k, ] %*% Bs[, t]
  }
}

y <- matrix(NA, N, T)
y[,1] <- X[t, , 1:2] %*% beta[t, 1:2] + rnorm(N)
for(t in 2:T) {
  X[t, ,3] <- y[, t-1]
  y[, t] <- rnorm(N, X[t, , ] %*% beta[t, ], sigma)
}


gaussian_example <- data.frame(y = c(y), x = c(t(X[,,2])), id = 1:N,
  time = rep(1:T, each = N))
usethis::use_data(gaussian_example, overwrite = TRUE)

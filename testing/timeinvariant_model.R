logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

library(rstan)

T <- 60
N <- 40
yg <- array(0, c(1, T, N))
yc <- array(0, c(2, T, N))
S <- c(2, 3)

alpha_1 <- c(0.2, 0)
alpha_2 <- c(0.4, -0.2, 0)
alpha_3 <- 0.1
K <- 4
Kj <- c(4, 4, 4)
x <- array(0, c(T, N, K))
ca <- list(y1 = contr.sum, y2 = contr.sum)
beta_1 <- rbind(matrix(c(1, -0.5, 0.2, 0.5), 1, 4), 0)
beta_2 <- rbind(matrix(c(
  0.4, 0.2, -0.3, -0.5,
  0.7, 0.1, -0.5, -0.2), 2, 4), 0)
beta_3 <- matrix(c(0.2, 0.4, -0.2, 0.3), 1, 4)
for(i in 1:N) {
  yc[1, 1, i] <- sample(1:2, size = 1, prob = c(0.6, 0.4))
  yc[2, 1, i] <- sample(1:3, size = 1, prob = c(0.4, 0.4, 0.2))
  yg[1, 1, i] <- rnorm(1)
}
for(t in 2:T) {
  y1 <- factor(yc[1, t-1,], levels = 1:2)
  y2 <- factor(yc[2, t-1,], levels = 1:3)
  y3 <- yg[1, t-1, ]
  x[t, , ] <- model.matrix(~ y1 + y2 + y3, contrasts.arg = ca)[, -1] # intercept handled separately
  for(i in 1:N) {
    yc[1, t, i] <- sample(1:2, size = 1, prob = softmax(alpha_1 + beta_1 %*% x[t, i,]))
    yc[2, t, i] <- sample(1:3, size = 1, prob = softmax(alpha_2 + beta_2 %*% x[t, i,]))
    yg[1, t, i] <- rnorm(1, alpha_3 + beta_3 %*% x[t, i,], 0.2)
  }
}
Xidx <- matrix(1:4,4,3)
d <- list(N = N, T = T, Pc = 2, S = S, yc = yc, K = K, x = x, Kj = Kj, Xidx = Xidx,
  Pg = 1, yg = yg)

model <- stan_model("inst/stan/time_constant_model.stan")

fit <- sampling(model, data = d, chains = 1, refresh=10)

pi_logit <- apply(extract(fit, "pi")[[1]], 2:3, mean)
softmax(pi_logit[1, 1:2])
softmax(pi_logit[2, ])

print(fit, "alpha_c")
alpha_1
alpha_2

print(fit, "beta_c")
beta_1

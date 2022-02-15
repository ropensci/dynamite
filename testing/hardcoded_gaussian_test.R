library(rstan)
rstan_options(javascript = FALSE)

logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
}

softmax <- function (x) {
    exp(x - logsumexp(x))
}


set.seed(1)
# number of time points
T <- 100
# number of individuals
N <- 100
K <- 1

internal_knots <- 80#min(5, T - 3)
knots <- seq(1, T, length.out = internal_knots + 2)[-c(1, internal_knots + 2)]
Bsplines <- t(splines::bs(1:T, knots = knots, degree = 3, intercept = TRUE))
D <- nrow(Bsplines)
matplot(t(Bsplines),type="l")

a_raw <- array(0, c(K, D))
a_raw[] <- rnorm(length(a_raw))
tau <- 0.2
sigma <- 0.1
beta <- array(0, c(T, K))
a <- array(0, c(K, D))

for(k in 1:K) {
    a[k, 1] <- a_raw[k, 1]
    for (i in 2:D) {
        a[k, i] = a[k, i-1] + tau[k] * a_raw[k, i]
    }
    for(t in 1:T){
        beta[t, k] = a[k, ] %*% Bsplines[, t]
    }
}

ts.plot(beta)

y <- matrix(NA, T, N)
X <- array(0, c(T, N, K))
X[] <- rnorm(length(X))
for(t in 1:T) {
    for(i in 1:N){
        y[t, i] <- rnorm(1, X[t, i, ] %*% t(beta[t, ]), sigma)
    }
}

ts.plot(y[,1:5])
####
model_centered <- stan_model("hardcoded_gaussian_case_centered.stan")
model_noncentered <- stan_model("hardcoded_gaussian_case_noncentered.stan")

d <- list(N = N, T = T, y = y, X = X,  K = K, B = Bsplines, D = D)

fit_c <- sampling(model_centered, data = d,
    chains = 4, cores = 4,
    refresh = 10, iter= 2000,
    seed = 1)
fit_nc <- sampling(model_noncentered, data = d,
    chains = 4, cores = 4,
    refresh = 10, iter= 2000,
    seed = 1)

print(fit_c, pars = c("sigma", "tau"))
print(fit_nc, pars = c("sigma", "tau"))

beta_estimates <- extract(fit_nc, "beta", permuted = TRUE)[[1]]

s <- 2
k <- 1
ts.plot(cbind(
    colMeans(beta_estimates[, , s, k]),
    apply(beta_estimates[, , s, k], 2, quantile, 0.025),
    apply(beta_estimates[, , s, k], 2, quantile, 0.975)), col = 1)
lines(beta[, s, k], col = 2)

# all in once:
ts.plot(cbind(c(aperm(a, 3:1)), get_posterior_mean(fit_nc, "a")[,1]), col = 1:2)

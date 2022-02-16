library(rstan)
rstan_options(javascript = FALSE)

logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
}

softmax <- function (x) {
    exp(x - logsumexp(x))
}


create_data <- function(T, N, S, K, n, tau = runif(K)) {

    X <- array(0, c(T, N, K))
    X[] <- rnorm(length(X), 1)

    knots <- seq(1, T, length.out = n + 2)[-c(1, n + 2)]
    Bsplines <- t(splines::bs(1:T, knots = knots, degree = 3, intercept = TRUE))
    D <- nrow(Bsplines)
    #matplot(t(Bsplines),type="l")

    a_raw <- array(0, c(S-1, K, D))
    a_raw[] <- rnorm(length(a_raw))
    beta <- array(0, c(T, S, K))
    a <- array(0, c(S-1, K, D))
    for(s in 1:(S-1)) {
        for(k in 1:K) {
            a[s, k, 1] <- a_raw[s, k, 1]
            for (i in 2:D) {
                a[s, k, i] = a[s, k, i-1] + a_raw[s, k, i] * tau[k]
            }
            for(t in 1:T){
                beta[t, s, k] = a[s, k, ] %*% Bsplines[, t]
            }
        }
    }
    y <- matrix(NA, T, N)
    for(t in 1:T) {
        for(i in 1:N){
            y[t, i] <- sample(1:S, size = 1, prob = softmax(X[t, i, ] %*% t(beta[t, , ])))
        }
    }

    TraMineR::seqdplot(TraMineR::seqdef(t(y)))
    list(N = N, T = T, y = y, X = X,  K = K, S = S, B = Bsplines, D = D, knots = knots,
        tau = tau, beta = beta, a = a)
}

model_centered <- stan_model("hardcoded_categorical_case_centered.stan")
model_noncentered <- stan_model("hardcoded_categorical_case_noncentered.stan")


set.seed(1)
d <- create_data(T = 100, N = 100, S = 3, K = 1, n = 10, tau = 0.5)

fit_c <- sampling(model_centered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

fit_nc <- sampling(model_noncentered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

# centered 10x faster, noncentered has bit better ESS
print(fit_c, pars = "tau")
print(fit_nc, pars = "tau")
ess_c <- summary(fit_c, pars = c("tau", "beta", "a"))[[1]][, "n_eff"]
ess_nc <- summary(fit_nc, pars = c("tau", "beta", "a"))[[1]][, "n_eff"]
summary(ess_c / ess_nc)

beta_estimates <- extract(fit_c, "beta", permuted = TRUE)[[1]]

ts.plot(cbind(
    colMeans(beta_estimates[, , 1, 1]),
    colMeans(beta_estimates[, , 2, 1]),
    d$beta[, 1, 1],
    d$beta[, 2, 1]), col = 1:2,lty=c(1,1,2,2))

ts.plot(cbind(c(aperm(d$a, 3:1)), get_posterior_mean(fit_c, "a")[,1]), col = 1:2)


# smaller tau
set.seed(1)
d <- create_data(T = 100, N = 100, S = 3, K = 1, n = 10, tau = 0.1)

fit_c <- sampling(model_centered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

fit_nc <- sampling(model_noncentered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

# centered 4x faster, noncentered better ESS?
print(fit_c, pars = "tau")
print(fit_nc, pars = "tau")
ess_c <- summary(fit_c, pars = c("tau", "beta"))[[1]][, "n_eff"]
ess_nc <- summary(fit_nc, pars = c("tau", "beta"))[[1]][, "n_eff"]
summary(ess_c / ess_nc)


# less time points
set.seed(1)
d <- create_data(T = 20, N = 100, S = 3, K = 1, n = 10, tau = 0.1)

fit_c <- sampling(model_centered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

fit_nc <- sampling(model_noncentered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

print(fit_c, pars = "tau")
print(fit_nc, pars = "tau")
ess_c <- summary(fit_c, pars = c("tau", "beta"))[[1]][, "n_eff"]
ess_nc <- summary(fit_nc, pars = c("tau", "beta"))[[1]][, "n_eff"]
summary(ess_c / ess_nc)


# less time points and individuals, centered version performs poorly
set.seed(1)
d <- create_data(T = 20, N = 10, S = 3, K = 1, n = 10, tau = 0.1)

fit_c <- sampling(model_centered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

fit_nc <- sampling(model_noncentered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

print(fit_c, pars = "tau")
print(fit_nc, pars = "tau")
ess_c <- summary(fit_c, pars = c("tau", "beta"))[[1]][, "n_eff"]
ess_nc <- summary(fit_nc, pars = c("tau", "beta"))[[1]][, "n_eff"]
summary(ess_c / ess_nc)



# larger data, centered took 385 seconds, noncentered took 6583 seconds!
set.seed(1)
d <- create_data(T = 100, N = 1000, S = 3, K = 1, n = 20, tau = 0.1)

fit_c <- sampling(model_centered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

fit_nc <- sampling(model_noncentered, data = d,
    chains = 1, cores = 1,
    refresh = 1000, iter= 2000,
    seed = 1)

print(fit_c, pars = "tau")
print(fit_nc, pars = "tau")
ess_c <- summary(fit_c, pars = c("tau", "beta"))[[1]][, "n_eff"]
ess_nc <- summary(fit_nc, pars = c("tau", "beta"))[[1]][, "n_eff"]
summary(ess_c / ess_nc)




## Without intercept use sum-to-zero-constraint


create_data0 <- function(T, N, S, K, n, tau = runif(K)) {

    X <- array(0, c(T, N, K))
    X[] <- rnorm(length(X))

    knots <- seq(1, T, length.out = n + 2)[-c(1, n + 2)]
    B <- splines::bs(1:T, knots = knots, degree = 3, intercept = FALSE)
    # sum-to-zero constraint
    Z <- qr.Q(qr(colSums(B)), complete = TRUE)[, 2:ncol(B)]
    Bsplines <- t(B %*% Z)
    D <- nrow(Bsplines)

    #matplot(t(Bsplines),type="l")

    a_raw <- array(0, c(S-1, K, D))
    a_raw[] <- rnorm(length(a_raw))
    beta <- array(0, c(T, S, K))
    beta0 <- cbind(matrix(rnorm(K * (S - 1)), S - 1, K), 0) # intercept for beta
    a <- array(0, c(S-1, K, D))
    for(s in 1:(S-1)) {
        for(k in 1:K) {
            a[s, k, 1] <- a_raw[s, k, 1]
            for (i in 2:D) {
                a[s, k, i] = a[s, k, i-1] + a_raw[s, k, i] * tau[k]
            }
            for(t in 1:T){
                beta[t, s, k] = beta0[s, k] + a[s, k, ] %*% Bsplines[, t]
            }
        }
    }
    y <- matrix(NA, T, N)
    for(t in 1:T) {
        for(i in 1:N){
            y[t, i] <- sample(1:S, size = 1, prob = softmax(X[t, i, ] %*% t(beta[t, , ])))
        }
    }

    TraMineR::seqdplot(TraMineR::seqdef(t(y)))
    list(N = N, T = T, y = y, X = X,  K = K, S = S, B = Bsplines, D = D, knots = knots,
        tau = tau, beta0 = beta0, beta = beta, a = a)
}

set.seed(1)
d <- create_data0(T = 50, N = 200, S = 3, K = 1, n = 10, tau = 0.1)


model_beta0 <- stan_model("hardcoded_categorical_case_centered_intercept.stan")

fit_beta0 <- sampling(model_beta0, data = d,
    chains = 1, cores = 1,
    refresh = 100, iter= 2000,
    seed = 1)

print(fit_beta0, "tau")
d$tau

print(fit_beta0, "beta0")
d$beta0

b <- extract(fit_beta0, "beta")[[1]]
mean(b[,,1,1])
mean(b[,,2,1])


## generate data as before but fit sum-to-zero version
set.seed(234)
d <- create_data(T = 100, N = 200, S = 3, K = 1, n = 10, tau = 0.5)

fit_c <- sampling(model_centered, data = d,
    chains = 1, cores = 1,
    refresh = 500, iter= 2000,
    seed = 1)

d2 <- d
# sum-to-zero constraint
Z <- qr.Q(qr(colSums(d$B)), complete = TRUE)[, 2:ncol(d$B)]
d2$B <- t(d$B %*% Z)
d2$D <- nrow(d2$B)
fit_beta0 <- sampling(model_beta0, data = d2,
    chains = 1, cores = 1,
    refresh = 500, iter= 2000,
    seed = 1)

print(fit_c,"tau")
print(fit_beta0,"tau")

b <- apply(extract(fit_c, "beta")[[1]], 2:3, mean)
b0 <- apply(extract(fit_beta0, "beta")[[1]], 2:3, mean)

ts.plot(cbind(b[,1:2], b0[,1:2], d$beta[, 1:2, 1]), col = rep(1:3, each = 2), lty = 1:2)

print(fit_beta0,"beta0")
colMeans(d$beta)
colMeans(b0)
colMeans(b)


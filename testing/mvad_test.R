library(TraMineR)
library(dplyr)
library(tidyr)
data(mvad)
n_ids <- 250
d <- pivot_longer(mvad[1:n_ids,], 15:86, "time")
d$time <- rep(1:72, length = nrow(d))
fit <- btvcm:::btvcmfit(
    obs(value ~ 1, family = categorical()) + lags() +
        splines(df = 5, intercept = TRUE, shrinkage = TRUE),
    d, id, chains = 1, refresh = 1)

print(fit$stanfit, "tau_1")
print(fit$stanfit, "alpha_1")
b <- apply(rstan::extract(fit$stanfit, "beta_1")[[1]], 2:4, mean)
ts.plot(b[, , 1])
ts.plot(b[, , 2])
ts.plot(b[, , 3])
ts.plot(b[, , 4])

d2 <- d %>% select(id, time, value)
# add fake series with random transitions
d2 <- rbind(d2, data.frame(id=713, time = 1:72, value = sample(unique(d$value), size = 72, replace = TRUE)))
fit2 <- btvcm:::btvcmfit(
    obs(value ~ 1, family = categorical()) + lags() +
        splines(intercept = TRUE, shrinkage = FALSE),
    d2, id, chains = 2, refresh = 100, iter = 4000)

print(fit2$stanfit, "tau_1")
print(fit2$stanfit, "alpha_1")
b <- apply(rstan::extract(fit3$stanfit, "beta_1")[[1]], 2:4, mean)
ts.plot(b[, , 1])
ts.plot(b[, , 2])
ts.plot(b[, , 3])
ts.plot(b[, , 4])


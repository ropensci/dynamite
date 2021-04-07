library(rstan)
rstan_options(javascript = FALSE)
# Simulate data
library(seqHMM)
set.seed(1)
init <- c(1, 0, 0)
A <- matrix(c(0.92, 0.08, 0,
  0, 0.98, 0.02,
  0, 0, 1), 3, 3, TRUE)

B <- matrix(c(0.6, 0.2, 0.1, 0.1, 0,
  0.2, 0.3, 0.3, 0.1, 0.1,
  0.1, 0.1, 0.1, 0.4, 0.3), 3, 5, TRUE)
sim <- simulate_hmm(100, init, A, B, sequence_length = 50)

TraMineR::seqIplot(sim$states)
TraMineR::seqdplot(sim$observation)

y <- t(data.matrix(sim$observations)) # convert to integers 1,2,3,4,5
data <- list(N = ncol(y), T = nrow(y), S = max(y), y = y, K = 3)

model_hmm <- stan_model("inst/stan/hmm.stan")
fit_hmm <- sampling(model_hmm, data = data,
  chains = 4, cores = 4, refresh = 100, iter=2000,
 seed = 1)
print(fit_hmm, c("A","B"))

model_hmm2 <- stan_model("inst/stan/hmm_common_states.stan")
fit_hmm2 <- sampling(model_hmm2, data = data,
  chains = 4, cores = 4, refresh = 100, iter=2000,
  seed = 1)
print(fit_hmm2, c("A","B"))

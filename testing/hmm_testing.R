library(rstan)
rstan_options(javascript = FALSE)
# Simulate data
library(seqHMM)
set.seed(1)
T <- 20
N <- 100
S <- 5
K <- 3
init <- c(1, 0, 0)
A <- matrix(c(0.92, 0.08, 0,
  0, 0.98, 0.02,
  0, 0, 1), K, K, TRUE)

B <- matrix(c(0.6, 0.2, 0.1, 0.1, 0,
  0.2, 0.3, 0.3, 0.1, 0.1,
  0.1, 0.1, 0.1, 0.4, 0.3), K, S, TRUE)
sim <- simulate_hmm(N, init, A, B, sequence_length = T)

TraMineR::seqIplot(sim$states)
TraMineR::seqdplot(sim$observation)

y <- t(data.matrix(sim$observations)) # convert to integers 1,2,3,4,5
data <- list(N = ncol(y), T = nrow(y), S = max(y), y = y, K = K, dirichlet_alpha = 1)

model_hmm <- stan_model("inst/stan/hmm.stan")
fit_hmm <- sampling(model_hmm, data = data,
  chains = 4, cores = 4, refresh = 100, iter=2000,
 seed = 1, pars = c("A", "B"))
print(fit_hmm, c("A","B"))

d <- data
d$first_B <- d$last_B <- rep(1, 5)
model_hmm2 <- stan_model("inst/stan/hmm_common_states.stan")
fit_hmm2 <- sampling(model_hmm2, data = data,
  chains = 4, cores = 4, refresh = 100, iter=2000,
  seed = 1, pars = c("A", "B"))
print(fit_hmm2, c("A","B"))



## One latent state process for all individuals
states <- as.numeric(sim$states[4,])
states <- rep(1:K, times = c(10,6,4))
y <- matrix(0, T, N)
for(t in 1:T) {
  y[t, ] <- sample(1:S, N, TRUE, B[states[t],])
}

TraMineR::seqdplot(seqdef(t(y)))

# for some reason now we get multimodality issues whereas with wrong model above
# everything worked ok...
# because we always start with state one, use first time points as prior for B[1]
# similarly we likely end up with last state in the end so prior for B[K]
# This helps with multimodality
first_B <- 5*(1 + table(factor(y[1:2, ], levels = 1:5)))
last_B <- 5*(1 + table(factor(y[(T-1):T, ], levels = 1:5)))
# some scaling depending on S, T, K or N maybe needed to keep prior tight but not too tight...

d <- list(N = ncol(y), T = nrow(y), S = max(y), y = y, K = 3, dirichlet_alpha = 1,
  first_B = first_B, last_B = last_B)
fit <- sampling(model_hmm2, data = d,
  chains = 4, cores = 4, refresh = 100, iter=2000,
  seed = 1)
print(fit, c("A","B"))
pairs(fit, pars = c("lp__", "A"))
pairs(fit, pars = c("lp__", "B[1,1]", "B[1,2]","B[2,1]", "B[2,2]"))

# Multimodality in B?? (was before prior for B[1] and B[K])
lapply(1:4, function(x) matrix(get_posterior_mean(fit, "B")[,x],3,5,TRUE))
get_posterior_mean(fit, "A") # maybe here as well, modes are just close to each other


## Using more general code
# number of symbols
S <- max(y)
# covariates, dummy coding so S - 1 covariates
K <- S - 1
x <- array(0, c(T, N, K))
for(t in 2:T) { # first time point fixed, so we actually start from T = 2 (remove t=1 later)
  x[t, , ] <- model.matrix( ~ factor(y[t - 1, ], levels = 1:S))[, -1] # intercept handled separately
}

d <- list(N = N, T = T - 1, Pc = 1, S = array(max(y), 1), yc = array(y[-1, ], c(1, T - 1, N)),
  K = K, x = x[-1, ,],
  Kj = array(0, 1), Xidx = matrix(0, K, 1),
  Kt = array(K, 1), Xidx_t = matrix(1:K, K, 1),
  Pg = 0, yg = array(0, c(0, T - 1, N)), n_states = 3)

model <- stan_model("inst/stan/hmm_states_model.stan")

fit <- sampling(model, data = d,
  chains = 1, cores = 1, refresh = 1, iter=2000,
  seed = 1)


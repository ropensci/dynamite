## code to create `multichannel_example` object

library(dynamite)
set.seed(1)
n_id <- 50L
n_time <- 20L

# first time point
d <- data.frame(
  g = rnorm(n_id),
  p = rpois(n_id, 4.0),
  b = rbinom(n_id, 1.0, 0.4),
  time = 1L,
  id = seq_len(n_id)
)

# rest of the time points
d <- dplyr::right_join(
  d,
  data.frame(
    time = rep(seq_len(n_time), each = n_id),
    id = seq_len(n_id)
  )
)

# formula for dynamite
f <- obs(g ~ lag(g) + lag(logp), family = "gaussian") +
  obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
  obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
  aux(numeric(logp) ~ log(p + 1))

# true values used for generating the data
alpha_g <- 0
beta_g <- c(0.8, 0.2)
a_g <- alpha_g + beta_g[1] * mean(d$g[d$time==1]) +
  beta_g[2] * mean(log(1 + d$p[d$time==1]))
sigma_g <- 1

alpha_p <- 0.1
beta_p <- c(-0.3, 0.7, 1)
a_p <- alpha_p + beta_p[1] * mean(d$g[d$time == 1]) +
  beta_p[2] * mean(log(1 + d$p[d$time == 1])) +
  beta_p[3] * mean(d$b[d$time == 1])

a_b <- 0.2
beta_b <- c(0.5, 0.2, 0.4, -0.6, 0.2)

# Run the model with the fixed parameters to obtain proper dynamitefit object
init <- dplyr::lst(a_g, beta_g, sigma_g, beta_p, a_p, beta_b, a_b)

fit <- dynamite(
  dformula = f,
  data = d,
  group = "id",
  time = "time",
  verbose = FALSE,
  chains = 1,
  iter = 1,
  algorithm = "Fixed_param",
  init = list(init)
)

# simulate data for time > 1
multichannel_example <- predict(fit, n_draws = 1) |>
  dplyr::mutate(g = g_new, p = p_new, b = b_new) |>
  dplyr::select(id, time, g, p, b)

usethis::use_data(multichannel_example, overwrite = TRUE)

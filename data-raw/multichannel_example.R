## code to create `multichannel_example` dataset
## code to create `multichannel_example` dataset

library(dynamite)
set.seed(1)
n_id <- 500
n_time <- 30
D <- 10
# first time point
d <- data.frame(g = rnorm(n_id),
                p = rpois(n_id, 5),
                b = rbinom(n_id, 1, 0.4),
                time = 1, id = 1:n_id)
# rest of the time points
d <- dplyr::right_join(d,
                       data.frame(
                         time = rep(1:n_time, each = n_id),
                         id = 1:n_id))
# known covariate
d$x <- rnorm(nrow(d))

f <-
  obs(g ~ -1 + lag(g) + lag(logp) + varying(~x), family = gaussian()) +
  obs(p ~ x + lag(g) + lag(logp) + lag(b), family = poisson()) +
  obs(b ~ lag(b) * lag(logp) + lag(b) * x + lag(b) * lag(g), family = bernoulli()) +
  aux(numeric(logp) ~ log(p + 1)) + splines(df = D)

# true values used for generating the data
tau_alpha_g <- 0.5
beta_g <- c(0.5, 0.5)
tau_g <- 0.5
sigma_g <- 0.5
a_g <- 0.1
omega_raw_alpha_g <- cumsum(rnorm(D - 1, 0, tau_alpha_g))
omega_g <- cumsum(c(-1, rnorm(D - 1, 0, tau_g)))

a_p <- -0.5
beta_p <- c(1, 0.2, 0.5, 0.1)

a_b <- 0.2
beta_b <- c(0.1, 0.2, 2, 0.4, -0.6, 0.2, -0.2)

# Run the model with the fixed parameters to obtain proper dynamitefit object
init <- dplyr::lst(beta_g, tau_g, omega_g, a_g, omega_raw_alpha_g,
                   tau_alpha_g, sigma_g, beta_p, a_p, beta_b, a_b)

fit <- dynamite(f, d, "id", "time",
                chains = 1, iter = 1,
                algorithm = "Fixed_param", init = list(init))

multichannel_example <- predict(fit, n_draws = 1) |>
  dplyr::mutate(g = g_new, p = p_new, b = b_new) |>
  dplyr::select(id, time, x, g, p, b)

usethis::use_data(multichannel_example, overwrite = TRUE)




## code to create `multichannel_example` dataset

set.seed(1)
n_id <- 50
n_time <- 20
d <- data.frame(g = rnorm(n_id),
                p = rpois(n_id, 10),
                b = rbinom(n_id, size = 1, prob = 0.3),
                time = 1, id = 1:n_id)
d <- dplyr::right_join(d,
                       data.frame(
                         time = rep(1:n_time, each = n_id),
                         id = 1:n_id))
d$x <- runif(nrow(d))


f <-
  obs(g ~ -1 + lag(g) + lag(logp) + varying(~x), family = gaussian()) +
  obs(p ~ x + lag(g) + lag(logp) + lag(b), family = poisson()) +
  obs(b ~ lag(logp) + lag(b) + varying(~ -1 + x + lag(g)), family = bernoulli()) +
  aux(logp ~ log(p + 1) + past(0)) +
  splines(df = 5)

get_priors(f, d, "id", "time")
init <- list(a_g = 2, tau_alpha_g = 0.5, beta_g = c(0.8, 0.2),
             delta_g = 2, tau_g = 0.2, sigma_g = 1,
             a_p = 2, beta_p = c(0.4, 0.2, 0.3, 1),
             a_b = 0, beta_b = c(0.1, 0.8), delta_b = c(0.1, -0.1),
             tau_b = c(0.5, 0.2))

fit <- dynamite(f,
                d, "id", "time",
                chains = 1, algorithm = "Fixed_param", init = list(init))

head(fit$data)

fulldata <- predict(fit, n_draws = 1)

library(ggplot2)
ggplot(fulldata, aes(time, g_new, group = id)) +
  geom_line()
ggplot(fulldata, aes(time, p_new, group = id)) +
  geom_line()
ggplot(fulldata, aes(time, b_new, group = id)) +
  geom_jitter(height=0.05,width=0.1, alpha=0.5)

## code to create `multichannel_example` dataset
library(dynamite)
set.seed(123)
n_id <- 100
n_time <- 40
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
  obs(b ~ lag(logp) + lag(b) + varying(~ -1 + x + lag(g)), family = bernoulli()) +
  aux(numeric(logp) ~ log(p + 1)) + splines(df = D)

# true values used for generating the data
beta_g <- c(0.7, 0.2)
a_g <- 0.1
tau_alpha_g <- 1
omega_raw_alpha_g <- cumsum(rnorm(D - 1, 0, tau_alpha_g))
ts.plot(splines::bs(1:n_time, df = D - 1) %*% omega_raw_alpha_g)
tau_g <- 0.2
omega_g <- cumsum(c(-1, rnorm(D - 1, 0, tau_g)))
sigma_g <- 0.1
beta_p <- c(0.4, 0.2, 0.3, 1)
a_p <- -0.1
a_b <- -0.4
beta_b <- c(0.1, 0.8)
tau_b <- c(0.5, 0.2)
omega_b <- rbind(cumsum(c(0.2, rnorm(D - 1, 0, tau_b[1]))),
                 cumsum(c(-0.1, rnorm(D - 1, 0, tau_b[2]))))

# Run the model with the fixed parameters to obtain proper dynamitefit object
init <- dplyr::lst(beta_g, tau_g, omega_g, a_g, omega_raw_alpha_g,
                   tau_alpha_g, sigma_g, beta_p, a_p, beta_b, omega_b,
                   tau_b, a_b)

fit <- dynamite(f, d, "id", "time",
                chains = 1, iter = 1,
                algorithm = "Fixed_param", init = list(init))

#
fulldata <- predict(fit, n_draws = 1)

library(ggplot2)
ggplot(fulldata, aes(time, g_new, group = id)) +
  geom_line()
ggplot(fulldata, aes(time, p_new, group = id)) +
  geom_line()
ggplot(fulldata, aes(time, b_new, group = id)) +
  geom_jitter(height = 0.05,width = 0.05, alpha = 0.5)


fulldata <- fulldata |>
  dplyr::mutate(g = g_new, p = p_new, b = b_new) |>
  dplyr::select(id, time, x, g, p, b)

fit2 <- dynamite(f,  fulldata, "id", "time",
                 chains = 1, cores = 1, init = 0, refresh = 10)

fit
fit2
plot_deltas(fit2, scales = "free")

plot_betas(fit2)

pred<-predict(fit2,type="mean")

pred |> dplyr::filter(id == 1 & time > 1) |>
  dplyr::group_by(time) |>
  dplyr::summarise(mean = mean(g_mean),
            ymin = quantile(g_mean, 0.05),
            ymax = quantile(g_mean, 0.95),
            obs = g[1]) |>
  ggplot(aes(time, mean)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax)) +
  geom_line() +
  geom_point(aes(y = obs), colour = "tomato")

pred |> dplyr::filter(id == 1 & time > 1) |>
  dplyr::group_by(time) |>
  dplyr::summarise(mean = mean(p_mean),
                   ymin = quantile(p_mean, 0.05),
                   ymax = quantile(p_mean, 0.95),
                   obs = p[1]) |>
  ggplot(aes(time, mean)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax)) +
  geom_line() +
  geom_point(aes(y = obs), colour = "tomato")

pred |> dplyr::filter(id == 1 & time > 1) |>
  dplyr::group_by(time) |>
  dplyr::summarise(mean = mean(b_mean),
                   ymin = quantile(b_mean, 0.05),
                   ymax = quantile(b_mean, 0.95),
                   obs = b[1]) |>
  ggplot(aes(time, mean)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax)) +
  geom_line() +
  geom_point(aes(y = obs), colour = "tomato")


pred<-predict(fit2)

pred |> dplyr::filter(id == 1 & time > 1) |>
  dplyr::group_by(time) |>
  dplyr::summarise(mean = mean(g_new),
                   ymin = quantile(g_new, 0.05),
                   ymax = quantile(g_new, 0.95),
                   obs = g[1]) |>
  ggplot(aes(time, mean)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax)) +
  geom_line() +
  geom_point(aes(y = obs), colour = "tomato")

pred |> dplyr::filter(id == 1 & time > 1) |>
  dplyr::group_by(time) |>
  dplyr::summarise(mean = mean(p_new),
                   ymin = quantile(p_new, 0.05),
                   ymax = quantile(p_new, 0.95),
                   obs = p[1]) |>
  ggplot(aes(time, mean)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax)) +
  geom_line() +
  geom_point(aes(y = obs), colour = "tomato")

pred |> dplyr::filter(id == 1 & time > 1) |>
  dplyr::group_by(time) |>
  dplyr::summarise(mean = as.integer(mean(b_new) > 0.5),
                   ymin = as.integer(quantile(b_new, 0.05) > 0.5),
                   ymax = as.integer(quantile(b_new, 0.95) > 0.5),
                   obs = b[1]) |>
  ggplot(aes(time, mean)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax),alpha=0.1) +
  geom_line() +
  geom_point(aes(y = obs), colour = "tomato")

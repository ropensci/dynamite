set.seed(1)
n_id <- 50
n_time <- 20
nu_x <- rnorm(n_id)
d <- data.frame(y = sample(factor(c("a", "b", "c", "d")), size = n_id, replace = TRUE),
                x = rnorm(n_id, nu_x),
                time = 1, id = 1:n_id)
d <- dplyr::right_join(d,
                       data.frame(
                         time = rep(1:n_time, each = n_id),
                         id = 1:n_id))

f <-
  obs(x ~ lag(x) + varying(~ -1 + lag(y)), family = gaussian(), random_intercept = TRUE) +
  obs(y ~ lag(x) + lag(y), family = categorical()) +
  splines(df = 10)

get_priors(f, d, "id", "time")
init <- list(sigma_nu_x = 1, nu_x = nu_x, beta_x = 0.8, tau_x = c(0.1,0.5,0.2),
             a_x = -1, sigma_x = 0.5,
             beta_y = matrix(rnorm(12), 3, 4), a_y = c(0.1, -0.5, 0.2))

fit <- dynamite(f,
                d, "id", "time",
                chains = 1, algorithm = "Fixed_param", init = list(init))

head(fit$data)

head(predict(fit))

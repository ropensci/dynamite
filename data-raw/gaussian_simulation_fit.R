# Code to create `gaussian_simulation_fit` object used in the
# `dynamite_simulation` vignette
library(dynamite)

set.seed(1)
n_id <- 100L
n_time <- 20L
d <- data.frame(
  y = rnorm(n_id, 1, 0.5),
  time = 1,
  id = seq_len(n_id)
)

d <- dplyr::right_join(
  d,
  data.frame(
    time = rep(seq_len(n_time), each = n_id),
    id = seq_len(n_id)
  ),
  by = c("time", "id")
)
d$x <- rnorm(nrow(d))

f <- obs(y ~ 1 + varying(~ -1 + x + lag(y)), family = "gaussian") +
  splines(df = 10)

init <- list(
  omega_y = rbind(
    c(0.0,  0.2, -0.8, 1.0, 0.5, 0.0, 0.1, -0.2, -0.5,  0.1),
    c(0.3, -0.4,  0.7, 0.5, 0.0, 0.0, 0.0,  0.6,  0.9, -0.8)
  ),
  tau_y = c(1.0, 0.75),
  a_y = -1,
  sigma_y = 0.5
)

gaussian_simulation_fit <- dynamite(
  dformula = f,
  data = d,
  time = "time",
  group = "id",
  chains = 1,
  iter = 1,
  algorithm = "Fixed_param",
  init = list(init),
)

usethis::use_data(
  gaussian_simulation_fit,
  overwrite = TRUE,
  compress = "xz"
)

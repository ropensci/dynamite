#' @srrstats {RE5.0, BS7.3} Scaling is approximately linear in number of time
#'  points, and less than linear in terms of groups. However, this is the best
#'  case scenario, as the performance of the NUTS algorithm depends on the
#'  posterior geometry. Also, the implementation of different probability
#'  distributions in Stan vary in terms of performance, e.g., gamma,
#'  exponential and categorical distributions can be slower than other
#'  supported distributions because of the lack of glm style parameterization.
#'  Finally, the splines and other more complex terms can likely affect the
#'  scaling in various ways.
#'

test_that("scaling for gaussian model is linear in number of time points", {

  skip_if_not(run_extended_tests)
  set.seed(1)
  N <- 10L
  T_ <- 1000L
  x <- matrix(rnorm(N * T_), N, T_)
  y <- matrix(NA, N, T_)
  y <- 2 + 2 * x + rnorm(N * T_)

  d <- data.frame(
    y = c(y), x = c(x), id = seq_len(N),
    time = rep(seq_len(T_), each = N))

  code <- get_code(obs(y ~ x, family = "gaussian"),
    data = d,
    group = "id",
    time = "time")
  model <- rstan::stan_model(model_code = code)

  n <- seq(100, 1000, length = 5)
  times <- numeric(length(n))
  for (i in seq_along(n)) {
    data <- get_data(obs(y ~ x, family = "gaussian"),
      data = d |> dplyr::filter(time <= n[i]),
      group = "id",
      time = "time")
    fit <- rstan::sampling(model, data = data,
       init = 0, refresh = 0, chains = 1, iter = 5000, seed = 1)
    times[i] <- sum(rstan::get_elapsed_time(fit))
  }
  expect_equal(min(times / n), max(times / n), tolerance = 0.1)
})

test_that("scaling for gaussian model is less linear in number of groups", {

  skip_if_not(run_extended_tests)
  set.seed(1)
  N <- 5000L
  T_ <- 50L
  x <- matrix(rnorm(N * T_), N, T_)
  y <- matrix(NA, N, T_)
  y <- 2 + 2 * x + rnorm(N * T_)

  d <- data.frame(
    y = c(y), x = c(x), id = seq_len(N),
    time = rep(seq_len(T_), each = N))

  code <- get_code(obs(y ~ x, family = "gaussian"),
    data = d,
    group = "id",
    time = "time")
  model <- rstan::stan_model(model_code = code)

  n <- c(1000, 5000)
  times <- numeric(length(n))
  for (i in seq_along(n)) {
    data <- get_data(obs(y ~ x, family = "gaussian"),
      data = d |> dplyr::filter(id <= n[i]),
      group = "id",
      time = "time")
    fit <- rstan::sampling(model, data = data,
      init = 0, refresh = 0, chains = 1, iter = 5000, seed = 1)
    times[i] <- sum(rstan::get_elapsed_time(fit))
  }
  expect_lt(times[2] / n[2], times[1] / n[1])
})


test_that("scaling for gamma model is linear in number of time points", {

  skip_if_not(run_extended_tests)
  set.seed(1)
  N <- 10L
  T_ <- 1000L
  x <- matrix(rnorm(N * T_), N, T_)
  y <- matrix(NA, N, T_)
  y <- rgamma(N * T_, 5, exp(2 + 2))

  d <- data.frame(
    y = c(y), x = c(x), id = seq_len(N),
    time = rep(seq_len(T_), each = N))

  code <- get_code(obs(y ~ x, family = "gamma"),
    data = d,
    group = "id",
    time = "time")
  model <- rstan::stan_model(model_code = code)

  n <- seq(100, 1000, length = 5)
  times <- numeric(length(n))
  for (i in seq_along(n)) {
    data <- get_data(obs(y ~ x, family = "gamma"),
      data = d |> dplyr::filter(time <= n[i]),
      group = "id",
      time = "time")
    fit <- rstan::sampling(model, data = data,
      init = 0, refresh = 0, chains = 1, iter = 5000, seed = 1)
    times[i] <- sum(rstan::get_elapsed_time(fit))
  }
  expect_equal(min(times / n), max(times / n), tolerance = 0.1)
})

test_that("scaling for gamma model is linear in number of groups", {

  skip_if_not(run_extended_tests)
  set.seed(1)
  N <- 5000L
  T_ <- 50L
  x <- matrix(rnorm(N * T_), N, T_)
  y <- matrix(NA, N, T_)
  y <- rgamma(N * T_, 5, exp(2 + 2))

  d <- data.frame(
    y = c(y), x = c(x), id = seq_len(N),
    time = rep(seq_len(T_), each = N))

  code <- get_code(obs(y ~ x, family = "gamma"),
    data = d,
    group = "id",
    time = "time")
  model <- rstan::stan_model(model_code = code)

  n <- seq(100, 1000, length = 5)
  times <- numeric(length(n))
  for (i in seq_along(n)) {
    data <- get_data(obs(y ~ x, family = "gamma"),
      data = d |> dplyr::filter(id <= n[i]),
      group = "id",
      time = "time")
    fit <- rstan::sampling(model, data = data,
      init = 0, refresh = 0, chains = 1, iter = 5000, seed = 1)
    times[i] <- sum(rstan::get_elapsed_time(fit))
  }
  expect_equal(min(times / n), max(times / n), tolerance = 0.1)
})

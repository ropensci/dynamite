run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

data.table::setDTthreads(1) # For CRAN

test_that("multivariate gaussian fit and predict work", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  N <- 100
  T_ <- 50
  S <- crossprod(matrix(rnorm(4), 2, 2))
  L <- t(chol(S))
  y1 <- matrix(0, N, T_)
  y2 <- matrix(0, N, T_)
  x <- matrix(0, N, T_)
  for (t in 2:T_) {
    for (i in 1:N){
      mu <- c(0.7 * y1[i, t-1], 0.4 * y2[i, t-1] - 0.2 * y1[i, t-1])
      y <- mu + L %*% rnorm(2)
      y1[i, t] <- y[1L]
      y2[i, t] <- y[2L]
      x[i, t] <- rnorm(1, c(0.5 * y1[i, t-1]), 0.5)
    }
  }
  d <- data.frame(
    y1 = c(y1),
    y2 = c(y2),
    x = c(x),
    t = rep(1:T_, each = N),
    id = 1:N
  )
  f <- obs(c(y1, y2) ~ -1 + lag(y1) | -1 + lag(y1) + lag(y2), "mvgaussian") +
    obs(x ~ -1 + lag(y1), "gaussian")
  fit <- dynamite(
    dformula = f,
    data = d,
    time = "t",
    group = "id",
    chains = 1,
    iter = 2000,
    refresh = 0
  )
  expect_error(sumr <- summary(fit, type = "corr"), NA)
  expect_equal(sumr$mean, cov2cor(S)[2, 1], tolerance = 0.1)
  expect_error(sumr <- summary(fit, type = "sigma"), NA)
  expect_equal(sumr$mean, c(0.5, sqrt(diag(S))), tolerance = 0.1)
  expect_error(sumr <- summary(fit, type = "beta"), NA)
  expect_equal(sumr$mean, c(0.5, 0.7, -0.2, 0.4), tolerance = 0.1)
  expect_error(predict(fit, n_draws = 5), NA)
})

test_that("multinomial fit and predict work", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  n_id <- 100L
  n_time <- 20L

  d <- data.frame(
    y1 = sample(10, size = n_id, replace = TRUE),
    y2 = sample(12, size = n_id, replace = TRUE),
    y3 = sample(15, size = n_id, replace = TRUE),
    time = 1,
    id = seq_len(n_id)
  )
  d$n <- d$y1 + d$y2 + d$y3

  d <- dplyr::right_join(
    d,
    data.frame(
      time = rep(seq_len(n_time), each = n_id),
      id = seq_len(n_id)
    )
  )

  d$z <- rnorm(nrow(d))
  d$n[is.na(d$n)] <- sample(37, size = sum(is.na(d$n)), replace = TRUE)

  f <- obs(
    c(y1, y2, y3) ~ z + lag(y1) + lag(y2) + lag(y3) + trials(n),
    family = "multinomial"
  )

  init <- list(
    beta_y2 = c(1.2, 0.8, 0.2, 0.1),
    beta_y3 = c(1, 0.5, 0.6, 0.3),
    a_y2 = -0.1,
    a_y3 = 0.2
  )

  expect_error(
    fit <- dynamite(
      dformula = f,
      data = d,
      time = "time",
      group = "id",
      chains = 1,
      iter = 1,
      algorithm = "Fixed_param",
      init = list(init),
    ),
    NA
  )
  expect_error(
    pred <- predict(fit, type = "response"),
    NA
  )
  expect_identical(
    pred$y1_new + pred$y2_new + pred$y3_new,
    pred$n
  )
})


test_that("cumulative fit works", {
  skip_if_not(run_extended_tests)

  set.seed(0)

  n <- 100
  t <- 30
  x <- matrix(0, n, t)
  y <- matrix(0, n, t)
  p <- matrix(0, n, 4)
  alpha <- c(-1, 0, 1)

  for (i in seq_len(t)) {
    x[, i] <- rnorm(n)
    eta <- 0.6 * x[, i]
    p[, 1] <- 1 - plogis(eta - alpha[1])
    p[, 2] <- plogis(eta - alpha[1]) - plogis(eta - alpha[2])
    p[, 3] <- plogis(eta - alpha[2]) - plogis(eta - alpha[3])
    p[, 4] <- plogis(eta - alpha[3])
    y[, i] <- apply(p, 1, sample, x = 1:4, size = 1, replace = FALSE)
  }

  d <- data.frame(
    y = factor(c(y), levels = 1:4), x = c(x),
    time = rep(seq_len(t), each = n),
    id = rep(seq_len(n), t)
  )

  expect_error(
    fit_logit <- dynamite(
      dformula =
        obs(y ~ x, family = "cumulative", link = "logit"),
      data = d,
      time = "time",
      group = "id"
    ),
    NA
  )
})

test_that("non-glm categorical fit works", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac") # Seems to segfault on MacOS
  expect_error(
    mockthat::with_mock(
      stan_supports_categorical_logit_glm = function(...) FALSE,
      dynamite(
        dformula = obs(x ~ z + lag(x) + lag(y), family = "categorical") +
          obs(y ~ z + lag(x) + lag(y), family = "categorical"),
        data = categorical_example,
        time = "time",
        group = "id",
        iter = 2000,
        chains = 1,
        refresh = 0,
        verbose = FALSE,
        threads_per_chain = 2
      )
    ),
    NA
  )
})

test_that("get_parameter_dims works for dynamiteformula", {
  skip_if_not(run_extended_tests)

  expect_error(
    dims <- get_parameter_dims(
      x = obs(x ~ z + lag(x) + lag(y), family = "categorical") +
        obs(y ~ z + lag(x) + lag(y), family = "categorical"),
      data = categorical_example,
      time = "time",
      group = "id"
    ),
    NA
  )
  expect_identical(
    dims,
    get_parameter_dims(categorical_example_fit)
  )
})

test_that("predict with random variable trials works", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  N <- 100
  T_ <- 50
  y <- matrix(0, N, T_)
  n <- matrix(0, N, T_)
  for(t in 2:T_) {
    for(i in 1:N) {
      n[i, t] <- rpois(1, 5)
      lp <- -3.0 + 1.25 * y[i, t-1]
      y[i, t] <- rbinom(1, 1 + n[i, t], plogis(lp))
    }
  }
  d <- data.frame(
    n = c(n),
    y = c(y),
    t = rep(1:T_, each = N),
    id = 1:N
  )
  f <- obs(n ~ 1, family = "poisson") +
    aux(integer(nn) ~ n + 1) +
    obs(y ~ lag(y) + trials(nn), family = "binomial")
  fit <- dynamite(
    dformula = f,
    data = d,
    time = "t",
    group = "id",
    chains = 1,
    iter = 2000,
    refresh = 0
  )
  expect_error(sumr <- summary(fit), NA)
  expect_equal(sumr$mean, c(log(5), -3.0, 1.25), tolerance = 0.1)
  expect_error(predict(fit, n_draws = 5), NA)
})

test_that("shrinkage for splines is functional", {
  skip("Shrinkage feature removed at least for now.")

  set.seed(1)
  expect_error(
    fit <- dynamite(
      dformula =
        obs(
          y ~ -1 + z + varying(~ x + lag(y)) + random(~1),
          family = "gaussian"
        ) +
        random_spec() +
        splines(df = 20, shrinkage = TRUE),
      data = gaussian_example,
      time = "time",
      group = "id",
      iter = 2000,
      warmup = 1000,
      chains = 1,
      refresh = 0
    ),
    NA
  )
  expect_error(
    summary(fit, types = "xi"),
    NA
  )
})

test_that("update without recompile works", {
  skip_if_not(run_extended_tests)

  set.seed(0)
  gaussian_fit <- dynamite(
    dformula =
      obs(
        y ~ -1 + z + varying(~ x + lag(y)) + random(~1),
        family = "gaussian"
      ) +
      random_spec() +
      splines(df = 20),
    data = gaussian_example,
    time = "time",
    group = "id",
    iter = 2000,
    warmup = 1000,
    thin = 1,
    chains = 2,
    cores = 2,
    refresh = 0,
    save_warmup = FALSE,
    pars = c(
      "omega_alpha_1_y", "omega_raw_alpha_y", "nu_raw", "nu", "L",
      "sigma_nu", "a_y"
    ),
    include = FALSE
  )
  expect_error(
    fit <- update(
      gaussian_fit,
      data = gaussian_example,
      warmup = 1000,
      iter = 2000,
      thin = 1
    ),
    NA
  )
  # Internal update_ function
  expect_error(
    lfo(gaussian_fit, L = 20, chains = 4, verbose_stan = FALSE),
    NA
  )
})

test_that("custom stan model works", {
  skip_if_not(run_extended_tests)

  # The same as gaussian_example_fit, but we need to refit because
  # the results may not be reproducible across different platforms
  set.seed(1)
  initial_fit <- dynamite(
    dformula = obs(y ~ -1 + z + varying(~ x + lag(y)) +
      random(~1), family = "gaussian") +
      random_spec() +
      splines(df = 20),
    data = gaussian_example,
    time = "time",
    group = "id",
    iter = 2000,
    warmup = 1000,
    thin = 10,
    chains = 2,
    cores = 2,
    refresh = 0,
    save_warmup = FALSE,
    pars = c(
      "omega_alpha_1_y", "omega_raw_alpha_y", "nu_raw", "nu", "L",
      "sigma_nu", "a_y"
    ),
    include = FALSE
  )
  code <- get_code(initial_fit)
  set.seed(1)
  expect_error(
    custom_fit <- dynamite(
      dformula = obs(y ~ -1 + z + varying(~ x + lag(y)) +
        random(~1), family = "gaussian") +
        random_spec() +
        splines(df = 20),
      data = gaussian_example,
      time = "time",
      group = "id",
      iter = 2000,
      warmup = 1000,
      thin = 10,
      chains = 2,
      cores = 2,
      refresh = 0,
      save_warmup = FALSE,
      pars = c(
        "omega_alpha_1_y", "omega_raw_alpha_y", "nu_raw", "nu", "L",
        "sigma_nu", "a_y"
      ),
      include = FALSE,
      custom_stan_model = code
    ),
    NA
  )
  expect_equal(
    rstan::extract(initial_fit$stanfit, permuted = FALSE),
    rstan::extract(custom_fit$stanfit, permuted = FALSE)
  )
})

test_that("dynamice works", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  n <- 50
  p <- 1000
  y <- replicate(p, stats::arima.sim(list(ar = 0.7), n, sd = 0.1))
  d <- data.frame(y = c(y), time = 1:n, id = rep(seq_len(p), each = n))
  dmiss <- d
  dmiss$y[sample(seq_len(nrow(d)), size = 0.2 * nrow(d))] <- NA

  # Long format imputation
  expect_error(
    fit_long <- dynamice(
      obs(y ~ lag(y), "gaussian"),
      time = "time",
      group = "id",
      data = dmiss,
      chains = 1,
      refresh = 0,
      backend = "rstan",
      impute_format = "long",
      keep_imputed = FALSE,
      mice_args = list(m = 5, print = FALSE)
    ),
    NA
  )

  # Wide format imputation
  expect_error(
    fit_wide <- dynamice(
      obs(y ~ lag(y), "gaussian"),
      time = "time",
      group = "id", data = dmiss,
      chains = 1,
      refresh = 0,
      backend = "rstan",
      impute_format = "wide",
      keep_imputed = FALSE,
      mice_args = list(m = 5, print = FALSE)
    ),
    NA
  )
})

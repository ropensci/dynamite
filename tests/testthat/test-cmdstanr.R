#' @srrstats {G5.10} Extended tests can be switched on via setting the
#'   environment variable DYNAMITE_EXTENDED_TESTS to "true".

run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

data.table::setDTthreads(1) # For CRAN

test_that("stanc_options argument works", {
  skip_if_not(run_extended_tests)
  set.seed(1)

  fit <- dynamite(
    dformula = obs(y ~ -1 + varying(~x), family = "gaussian") +
      lags(type = "varying") +
      splines(df = 20),
    gaussian_example,
    "time",
    "id",
    #parallel_chains = 2,
    chains = 1,
    refresh = 0,
    backend = "cmdstanr",
    stanc_options = list("O0"),
    show_messages = FALSE,
    init = 0
  )
  expect_equal(
    summary(fit, parameter = "alpha_y")$mean[2],
    1.5,
    tolerance = 0.1,
    ignore_attr = TRUE
  )
})

test_that("cmdstanr backend works for categorical model", {
  skip_if_not(run_extended_tests)
  set.seed(1)

  fit_dynamite <- update(
    categorical_example_fit,
    stanc_options = list("O0"),
    backend = "cmdstanr",
    show_messages = FALSE
  )
  expect_equal(
    coef(fit_dynamite)$mean[1L], -0.5,
    tolerance = 0.1,
    ignore_attr = TRUE
  )
  expect_error(get_code(fit_dynamite), NA)
})

test_that("LOO and LFO works for AR(1) model estimated with cmdstanr", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  fit <- dynamite(
    dformula = obs(LakeHuron ~ 1, "gaussian") + lags(),
    data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
    time = "time",
    group = "id",
    chains = 1,
    iter_sampling = 1000,
    iter_warmup = 1000,
    refresh = 0,
    backend = "cmdstanr",
    stanc_options = list("O0"),
    show_messages = FALSE,
    init = 0
  )
  l <- loo(fit)
  expect_equal(
    l$estimates,
    structure(
      c(
        -107.877842970846, 2.86041434691809, 215.755685941693,
        7.36848739076899, 0.561813071004331, 14.736974781538
      ),
      dim = 3:2,
      dimnames = list(c("elpd_loo", "p_loo", "looic"), c("Estimate", "SE"))
    ),
    tolerance = 1
  )
  expect_error(plot(l), NA)

  l <- lfo(fit, L = 20)
  expect_equal(l$ELPD, -91.5427292742266, tolerance = 1)
  expect_equal(l$ELPD_SE, 7.57777033647258, tolerance = 1)
  expect_error(plot(l), NA)
})

test_that("within-chain parallelization with cmdstanr works", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac") # Seems to segfault on MacOS
  set.seed(1)
  f <- obs(g ~ lag(g) + lag(logp), family = "gaussian") +
    obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
    obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
    aux(numeric(logp) ~ log(p + 1) | init(0))
  fit_dynamite <- dynamite(
    dformula = f,
    data = multichannel_example,
    time = "time",
    group = "id",
    backend = "cmdstanr",
    show_messages = FALSE,
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1
  )
  expect_equal(
    coef(fit_dynamite)$mean[1L],
    0.003,
    tolerance = 0.1,
    ignore_attr = TRUE
  )
  expect_error(get_code(fit_dynamite), NA)
})

test_that("multivariate gaussian with threading produces a valid model", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac")
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
      mu <- c(0.7 * y1[i, t - 1], 0.4 * y2[i, t - 1] - 0.2 * y1[i, t - 1])
      y <- mu + L %*% rnorm(2)
      y1[i, t] <- y[1L]
      y2[i, t] <- y[2L]
      x[i, t] <- rnorm(1, c(0.5 * y1[i, t - 1]), 0.5)
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
  code <- get_code(
    x = f,
    data = d,
    time = "t",
    group = "id",
    backend = "cmdstanr",
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())
})

test_that("threading produces valid model code for other distributions", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac")

  set.seed(1)
  n_id <- 50L
  n_time <- 20L
  n <- n_id * n_time
  d <- data.frame(
    g = rgamma(n, 3, 0.5),
    p = rpois(n, 5.0),
    b = rbinom(n, 10, 0.4),
    s = rt(n, 15),
    bb = rbeta(n, 3, 6),
    e = rexp(n, 2),
    time = rep(seq_len(n_time), each = n_id),
    id = rep(seq_len(n_id)),
    trials = rep(10, n)
  )
  f <- obs(g ~ lag(g), family = "gamma") +
    obs(p ~ lag(g) + lag(b), family = "negbin") +
    obs(b ~ lag(b) + lag(b) * lag(g) + trials(trials), family = "binomial") +
    obs(s ~ lag(g), family = "student") +
    obs(bb ~ lag(b), family = "beta") +
    obs(e ~ lag(s), family = "exponential")
  code <- get_code(
    x = f,
    data = d,
    time = "time",
    group = "id",
    backend = "cmdstanr",
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())
})

test_that("threading produces a valid model for cumulative", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac")

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

  code <- get_code(
    x = obs(y ~ x, family = "cumulative", link = "logit"),
    data = d,
    time = "time",
    group = "id",
    backend = "cmdstanr",
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())

  code <- get_code(
    x = obs(y ~ -1 + x + varying(~1), family = "cumulative", link = "logit") +
      splines(df = 10),
    data = d,
    time = "time",
    group = "id",
    backend = "cmdstanr",
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())

  # no predictors
  code <- get_code(
    x = obs(y ~ -1 + varying(~1), family = "cumulative", link = "logit") +
      splines(df = 10),
    data = d,
    time = "time",
    group = "id",
    backend = "cmdstanr",
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())

  code <- get_code(
    x = obs(y ~ 1, family = "cumulative", link = "logit") +
      splines(df = 10),
    data = d,
    time = "time",
    group = "id",
    backend = "cmdstanr",
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())
})

test_that("threading produces a valid model for multinomial", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac")

  set.seed(1)
  n_id <- 100L
  n_time <- 20L
  d <- data.frame(
    y1 = sample(10, size = n_id * n_time, replace = TRUE),
    y2 = sample(15, size = n_id * n_time, replace = TRUE),
    y3 = sample(20, size = n_id * n_time, replace = TRUE),
    z = rnorm(n_id * n_time),
    time = seq_len(n_time),
    id = rep(seq_len(n_id), each = n_time)
  )
  d$n <- d$y1 + d$y2 + d$y3
  f <- obs(
    c(y1, y2, y3) ~ z + lag(y1) + lag(y2) + lag(y3) + trials(n),
    family = "multinomial"
  )
  code <- get_code(
    x = f,
    data = d,
    time = "time",
    group = "id",
    backend = "cmdstanr",
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())
})

test_that("syntax is correct for various models", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac")

  # categorical_example_fit
  code <- get_code(
    x = obs(x ~ z + lag(x) + lag(y), family = "categorical") +
      obs(y ~ z + lag(x) + lag(y), family = "categorical"),
    data = categorical_example,
    time = "time",
    group = "id",
    backend = "cmdstanr"
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())

  # gaussian_example_fit
  code <- get_code(
    x = obs(
      y ~ -1 + z + varying(~ x + lag(y)) + random(~1), family = "gaussian"
    ) +
      random_spec() +
      splines(df = 20),
    data = gaussian_example,
    time = "time",
    group = "id",
    backend = "cmdstanr"
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())

  # categorical_example_fit
  code <- get_code(
    x = obs(g ~ lag(g) + lag(logp), family = "gaussian") +
      obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
      obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
      aux(numeric(logp) ~ log(p + 1)| init(0)),
    data = multichannel_example,
    time = "time",
    group = "id",
    backend = "cmdstanr"
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())

  # ordered probit model
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

  code <- get_code(
    x = obs(y ~ x, family = "cumulative", link = "logit"),
    data = d,
    time = "time",
    group = "id",
    backend = "cmdstanr"
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(code)
  model <- with(e, {
    cmdstanr::cmdstan_model(file, compile = FALSE)
  })
  expect_true(model$check_syntax())
})

test_that("latent factor syntax is correct", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac")

  set.seed(123)
  N <- 40L
  T_ <- 20L
  D <- 10
  B <- t(splines::bs(1:T_, df = D, intercept = TRUE))
  z <- rnorm(D)
  a <- cumsum(z) + rnorm(D)
  b <- cumsum(z) + rnorm(D)
  psi <- matrix(0, 2, T_)
  lambda_yi <- rnorm(N, 1, 0.2)
  lambda_xi <- rnorm(N, 1, 0.2)
  for (t in 1:T_) {
    psi[1, t] <- B[, t] %*% a
    psi[2, t] <- B[, t] %*% b
  }
  y <- matrix(0, N, T_)
  x <- matrix(0, N, T_)
  for (t in 1:T_) {
    y[, t] <- rnorm(N, lambda_yi * psi[1, t], 0.2)
    x[, t] <- rnorm(N, lambda_xi * psi[2, t], 0.2)
  }
  latent_factor_example <- data.frame(
    y = c(y),
    x = c(y),
    id = seq_len(N),
    time = rep(seq_len(T_), each = N)
  )
  lfactor_opts <- expand.grid(
    nonzero_lambda = c(FALSE, TRUE),
    noncentered_psi = c(FALSE, TRUE),
    correlated = c(FALSE, TRUE)
  )
  for (i in seq_len(nrow(lfactor_opts))) {
    code <- get_code(
      x = obs(y ~ 1, family = "gaussian") +
        obs(x ~ 1, family = "gaussian") +
        lfactor(
          nonzero_lambda = lfactor_opts$nonzero_lambda[i],
          noncentered_psi = lfactor_opts$noncentered_psi[i],
          correlated = lfactor_opts$correlated[i]
        ) +
        splines(df = 10),
      data = latent_factor_example,
      group = "id",
      time = "time"
    )
    e <- new.env()
    e$file <- cmdstanr::write_stan_file(code)
    model <- with(e, {
      cmdstanr::cmdstan_model(file, compile = FALSE)
    })
    expect_true(model$check_syntax())
  }
})

test_that("dynamice with cmdstanr backend works", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac")

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
      dformula = obs(y ~ lag(y), "gaussian"),
      time = "time",
      group = "id",
      data = dmiss,
      chains = 1,
      refresh = 0,
      backend = "cmdstanr",
      impute_format = "long",
      keep_imputed = FALSE,
      mice_args = list(m = 3, print = FALSE)
    ),
    NA
  )
})

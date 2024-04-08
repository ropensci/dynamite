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
    parallel_chains = 2,
    chains = 2,
    refresh = 0,
    backend = "cmdstanr",
    stanc_options = list("O0"),
    show_messages = FALSE,
    init = 0
  )
  expect_equal(summary(fit, parameter = "alpha_y")$mean[2], 1.5,
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
  fit <- dynamite(obs(LakeHuron ~ 1, "gaussian") + lags(),
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
  expect_equal(l$estimates,
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
    aux(numeric(logp) ~ log(p + 1))
  fit_dynamite <- dynamite(
    f,
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
  out <- dynamite(
    dformula = f,
    data = d,
    time = "t",
    group = "id",
    backend = "cmdstanr",
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1,
    debug = list(no_compile = TRUE, model_code = TRUE)
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(out$model_code)
  model <- with(e, {cmdstanr::cmdstan_model(file, compile = FALSE)})
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
    time = rep(seq_len(n_time), each = n_id),
    id = rep(seq_len(n_id)),
    trials = rep(10, n)
  )
  f <- obs(g ~ lag(g), family = "gamma") +
    obs(p ~ lag(g) + lag(b), family = "negbin") +
    obs(b ~ lag(b) + lag(b) * lag(g) + trials(trials), family = "binomial") +
    obs(s ~ lag(g), family = "student") +
    obs(bb ~ lag(b), family = "beta")
  out <- dynamite(
    dformula = f,
    data = d,
    time = "time",
    group = "id",
    backend = "cmdstanr",
    threads_per_chain = 2,
    grainsize = 10,
    chains = 2,
    parallel_chains = 1,
    debug = list(no_compile = TRUE, model_code = TRUE)
  )
  e <- new.env()
  e$file <- cmdstanr::write_stan_file(out$model_code)
  model <- with(e, {cmdstanr::cmdstan_model(file, compile = FALSE)})
  expect_true(model$check_syntax())
})

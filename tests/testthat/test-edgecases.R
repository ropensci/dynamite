#' @srrstats {G5.8,G5.8a, G5.8b, G5.8c, G5.8d} Edge conditions are tested.

# Model edgecases ---------------------------------------------------------

set.seed(0)
timepoints <- 10
individuals <- 5
total_obs <- timepoints * individuals

test_data <- data.frame(
  time = 1:timepoints,
  group = gl(individuals, timepoints),
  offset = sample(50:100, size = total_obs, replace = TRUE),
  trials = sample(50:100, size = total_obs, replace = TRUE)
) |>
  dplyr::mutate(
    y1 = as.factor(sample(5, size = total_obs, replace = TRUE)),
    y2 = rnorm(n = total_obs, mean = 1, sd = 2),
    y3 = rbinom(n = total_obs, size = trials, prob = 0.75),
    y4 = rbinom(n = total_obs, size = 1, prob = 0.66),
    y5 = rnbinom(n = total_obs, size = 100, prob = 0.33),
    y6 = rpois(n = total_obs, lambda = log(offset) + 1),
    y7 = rexp(n = total_obs, rate = 0.1),
    y8 = rgamma(n = total_obs, shape = 2, rate = 2 * exp(-5)),
    y9 = rbeta(n = total_obs, 6, 4),
    x1 = sample(letters[1:4], size = total_obs, replace = TRUE),
    x2 = rnorm(total_obs),
    x3 = as.factor(sample(4, size = total_obs, replace = TRUE))
  )

debug <- list(no_compile = TRUE)

test_that("single channel models are valid", {
  expect_error(
    obs_categorical <- obs(y1 ~ x1, family = "categorical"), NA
  )
  expect_error(
    obs_gaussian <- obs(y2 ~ x2, family = "gaussian"), NA
  )
  expect_error(
    obs_binomial <- obs(y3 ~ x3 + trials(trials), family = "binomial"), NA
  )
  expect_error(
    obs_bernoulli <- obs(y4 ~ x1, family = "bernoulli"), NA
  )
  expect_error(
    obs_negbinom <- obs(y5 ~ x2, family = "negbin"), NA
  )
  expect_error(
    obs_poisson <- obs(y6 ~ x3 + offset(offset), family = "poisson"), NA
  )
  expect_error(
    obs_exponential <- obs(y7 ~ x3, family = "exponential"), NA
  )
  expect_error(
    obs_gamma <- obs(y8 ~ x1, family = "gamma"), NA
  )
  expect_error(
    obs_beta <- obs(y9 ~ x3, family = "beta"), NA
  )
  expect_error(
    dynamite(obs_categorical, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_gaussian, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_binomial, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_bernoulli, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_negbinom, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_poisson, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_exponential, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_gamma, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_beta, test_data, "group", "time", debug = debug), NA
  )
})

test_that("multichannel models are valid", {
  expect_error(
    obs_all <- obs(y1 ~ x1, family = "categorical") +
      obs(y2 ~ x2, family = "gaussian") +
      obs(y3 ~ x3 + trials(trials), family = "binomial") +
      obs(y4 ~ x1, family = "bernoulli") +
      obs(y5 ~ x2, family = "negbin") +
      obs(y6 ~ x3 + offset(offset), family = "poisson") +
      obs(y7 ~ x3, family = "exponential") +
      obs(y8 ~ x1, family = "gamma") +
      obs(y9 ~ x1, family = "beta"),
    NA
  )
  expect_error(
    dynamite(obs_all, test_data, "group", "time", debug = debug), NA
  )
})

test_that("binomial with trials only works", {
  expect_error(
    obs_binomial <- obs(y3 ~ trials(trials), family = "binomial"), NA
  )
  expect_error(
    dynamite(obs_binomial, test_data, "group", "time", debug = debug), NA
  )
})

test_that("intercepts are handled correctly", {
  expect_error(
    obs_all_alpha <- obs(y1 ~ -1 + x1, family = "categorical") +
      obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ -1 + x3 + varying(~x1) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      splines(df = 5),
    NA
  )
  expect_error(
    dynamite(obs_all_alpha, test_data, "group", "time", debug = debug), NA
  )
})

test_that("random intercepts are handled correctly", {
  expect_error(
    obs_all_alpha <- obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ -1 + x3 + varying(~x1) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      splines(df = 5) + random(c("y2", "y3")),
    NA
  )
  expect_equal(
    c("sigma_nu_y2", "sigma_nu_y3", "L"),
    get_priors(obs_all_alpha, test_data, "group", "time")$parameter[c(1, 6, 24)]
  )

  expect_error(
    obs_all_alpha <- obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ -1 + x3 + varying(~x1) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      splines(df = 5) + random(correlated = FALSE),
    NA
  )
  expect_error(
    obs_all_alpha <- obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ -1 + x3 + varying(~x1) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      splines(df = 5) + random(correlated = FALSE, noncentered = FALSE),
    NA
  )
  expect_error(
    obs_all_alpha <- obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ -1 + x3 + varying(~x1) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      splines(df = 5) + random(correlated = TRUE, noncentered = FALSE),
    NA
  )
  expect_error(
    obs_all_alpha <- obs(y2 ~ x2, family = "gaussian") +
      random(noncentered = FALSE),
    NA
  )
  expect_false(
    "L" %in% get_priors(obs_all_alpha, test_data, "group", "time")$parameter
  )
})

test_that("shrinkage is handled correctly", {
  expect_error(
    obs_all_alpha <- obs(y1 ~ -1 + varying(~x1), family = "categorical") +
      obs(x3 ~ varying(~ -1 + x1), family = "categorical") +
      obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ lag(x3) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      obs(y9 ~ -1 + x1 + varying(~x2), family = "beta") +
      splines(df = 5, shrinkage = TRUE),
    NA
  )

  expect_equal(
    unname(unlist(lapply(
      dynamite(obs_all_alpha, test_data,
        "group", "time",
        debug = debug
      )$stan$model_vars, "[[", "shrinkage"
    ))),
    rep(TRUE, 6)
  )
  expect_equal(get_priors(
    obs_all_alpha, test_data,
    "group", "time"
  )$parameter[58], "xi")
})

test_that("noncentered splines are handled correctly", {
  expect_error(
    obs_all_alpha <- obs(y1 ~ -1 + varying(~x1), family = "categorical") +
      obs(x3 ~ varying(~ -1 + x1), family = "categorical") +
      obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ lag(x3) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      obs(y9 ~ -1 + x1 + varying(~x2), family = "beta") +
      splines(df = 5, noncentered = rep(TRUE, 6)),
    NA
  )

  expect_equal(
    unname(unlist(lapply(
      dynamite(obs_all_alpha, test_data,
        "group", "time",
        debug = debug
      )$stan$model_vars, "[[", "noncentered"
    ))),
    rep(TRUE, 6)
  )
})

test_that("lower bounds for tau are handled correctly", {
  expect_error(
    obs_all_alpha <- obs(y1 ~ -+x1, family = "categorical") +
      obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ -1 + x3 + varying(~x1) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      obs(y9 ~ -1 + x1 + varying(~x2), family = "beta") +
      splines(df = 5, lb_tau = c(0, 2, 1, 0.4, 0.01)),
    NA
  )

  expect_equal(
    lapply(
      dynamite(obs_all_alpha, test_data,
        "group", "time",
        debug = debug
      )$stan$model_vars, "[[", "lb"
    ),
    list(y1 = 0, y2 = 2, y3 = 1, y4 = 0.4, y9 = 0.01)
  )
})

test_that("manual fixed() terms work", {
  expect_error(
    obs_fixed <- obs(y1 ~ fixed(~ x1 + lag(y2, 1)), family = "categorical"),
    NA
  )
})

test_that("latent factors are handled correctly", {
  expect_error(
    obs_all_lfactor <- obs(y2 ~ -1 + x2, family = "gaussian") +
      obs(y3 ~ -1 + x3 + varying(~x1) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      splines(df = 5) + lfactor(c("y2", "y3"),
        nonzero_lambda = c(TRUE, FALSE, FALSE)),
    NA
  )
  expect_equal(
    c("sigma_lambda_y2", "tau_psi_y2", "psi_y2",
      "sigma_lambda_y3", "psi_y3", "L_lf"),
    get_priors(obs_all_lfactor, test_data, "group", "time")$parameter[
      c(1:3, 6:7, 25)]
  )
})
# Lag edgecases -----------------------------------------------------------

test_that("lags are parsed", {
  expect_error(
    obs_a <- obs(y1 ~ x1 + lag(y2, 1), family = "categorical") +
      obs(y2 ~ x2 + lag(y1, 1), family = "gaussian"),
    NA
  )
  expect_error(
    obs_b <- obs(y1 ~ -1 + x1 + varying(~ lag(y2, 1)),
      family = "categorical"
    ) +
      obs(y2 ~ -1 + x2 + varying(~ lag(y1, 1)), family = "gaussian") +
      splines(),
    NA
  )
  expect_error(
    obs_c <- obs(y1 ~ x1, family = "categorical") +
      obs(y2 ~ x2, family = "gaussian") +
      lags(k = 1, type = "fixed"),
    NA
  )
  expect_error(
    obs_d <- obs(y1 ~ x1, family = "categorical") +
      obs(y2 ~ x2, family = "gaussian") +
      lags(k = 1, type = "varying") +
      splines(),
    NA
  )
  expect_error(
    dynamite(obs_a, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_b, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_c, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_d, test_data, "group", "time", debug = debug), NA
  )
})

test_that("lag without explicit shift equal shift of one", {
  expect_identical(extract_lags("lag(y)")$k, 1L)
})

test_that("lags() and lag() give equal results", {
  f1 <- expect_error(
    obs(y1 ~ x1, family = "categorical") +
      obs(y2 ~ x2, family = "gaussian") +
      lags(k = 1, type = "fixed"),
    NA
  )
  f2 <- expect_error(
    obs(y1 ~ x1 + lag(y1) + lag(y2), family = "categorical") +
      obs(y2 ~ x2 + lag(y1) + lag(y2), family = "gaussian"),
    NA
  )
  expect_identical(
    get_priors(f1, test_data, "group", "time"),
    get_priors(f2, test_data, "group", "time")
  )
})

test_that("higher order lags() and lag() give equal results", {
  f1 <- obs(y1 ~ x1, family = "categorical") +
    obs(y2 ~ x2, family = "gaussian") +
    lags(k = 1:2, type = "fixed")
  f2 <- obs(y1 ~ x1 + lag(y1, 1) + lag(y2, 1) +
    lag(y1, 2) + lag(y2, 2), family = "categorical") +
    obs(y2 ~ x2 + lag(y1, 1) + lag(y2, 1) +
      lag(y1, 2) + lag(y2, 2), family = "gaussian")
  expect_identical(
    get_priors(f1, test_data, "group", "time"),
    get_priors(f2, test_data, "group", "time")
  )
})

test_that("disjoint and intersecting lags() ang lag() give equal results", {
  f1 <- obs(y2 ~ x2, family = "gaussian") +
    obs(y7 ~ x2 + lag(y2), family = "exponential") +
    lags(k = 1L)
  f2 <- obs(y2 ~ x2, family = "gaussian") +
    obs(y7 ~ x2, family = "exponential") +
    lags(k = 1L)
  expect_identical(
    get_priors(f1, test_data, "group", "time"),
    get_priors(f2, test_data, "group", "time")
  )
})

test_that("lags() is ignored if all lags already exist", {
  f1 <- obs(y2 ~ x2 + lag(y2) + lag(y7), family = "gaussian") +
    obs(y7 ~ x2 + lag(y2) + lag(y7), family = "exponential") +
    lags(k = 1L)
  f2 <- obs(y2 ~ x2 + lag(y2) + lag(y7), family = "gaussian") +
    obs(y7 ~ x2 + lag(y2) + lag(y7), family = "exponential")
  expect_identical(
    get_priors(f1, test_data, "group", "time"),
    get_priors(f2, test_data, "group", "time")
  )
})

# Data edgecases ----------------------------------------------------------

test_that("data expansion to full time scale works", {
  set.seed(0)
  # groups
  mis_rows <- sample(2L:(nrow(test_data) - 1L), 10)
  test_data_mis <- test_data[-mis_rows, ]
  fit <- dynamite(
    obs(
      y2 ~ x1 + x2 + x3 + y1 + y3 + y4 + y5 + y6 + y7 + y8 + y9,
      family = "gaussian"
    ),
    data = test_data_mis,
    group = "group",
    time = "time",
    verbose = FALSE,
    debug = debug
  )
  expected_data <- test_data
  expected_data[mis_rows, seq(3, ncol(test_data))] <- NA
  expected_data$x1 <- factor(expected_data$x1)
  expected_data <- droplevels(expected_data)
  expected_data$trials <- NULL
  expected_data$offset <- NULL
  data.table::setDT(expected_data, key = c("group", "time"))
  expect_equal(fit$data, expected_data, ignore_attr = TRUE)
  # no groups
  test_data_single <- test_data[test_data[["group"]] == 1L, ]
  mis_rows_single <- c(3, 5, 6, 9)
  test_data_single_mis <- test_data_single[-mis_rows_single, ]
  fit_single <- dynamite(
    obs(
      y2 ~ x1 + x2 + x3 + y1 + y3 + y4 + y5 + y6 + y7 + y8 + y9,
      family = "gaussian"
    ),
    data = test_data_single_mis,
    time = "time",
    verbose = FALSE,
    debug = debug
  )
  expected_data_single <- test_data_single
  expected_data_single[mis_rows_single, seq(3, ncol(test_data_single))] <- NA
  expected_data_single$x1 <- factor(expected_data_single$x1)
  expected_data_single <- droplevels(expected_data_single)
  expected_data_single$trials <- NULL
  expected_data_single$offset <- NULL
  data.table::setDT(expected_data_single, key = c("time"))
  expect_equal(fit_single$data, expected_data_single, ignore_attr = TRUE)
})

test_that("no groups data preparation works", {
  test_data_single <- test_data |>
    dplyr::filter(.data$group == 1) |>
    dplyr::select(!"group")
  obs_all <- obs(y1 ~ x2 + lag(y1), family = "categorical") +
    obs(y2 ~ x2, family = "gaussian") +
    obs(y3 ~ x3 + trials(trials), family = "binomial") +
    obs(y4 ~ x1, family = "bernoulli") +
    obs(y5 ~ x2, family = "negbin") +
    obs(y6 ~ x3 + offset(offset), family = "poisson") +
    obs(y7 ~ x3, family = "exponential") +
    obs(y8 ~ x1, family = "gamma")
  expect_error(
    dynamite(obs_all, test_data_single, time = "time", debug = debug),
    NA
  )
})

# Deterministic edgecases -------------------------------------------------

test_that("deterministic channels are parsed", {
  expect_error(
    obs_det <- obs(y5 ~ x1 + lag(d, 1) + lag(y5, 1) + lag(x1, 1),
      family = "negbin"
    ) +
      aux(numeric(d) ~ lag(d, 1) + lag(f, 2) + x2 | init(0)) +
      aux(numeric(f) ~ lag(y5, 1) + x2 * 3 + 1 | init(c(0, 1))),
    NA
  )
  expect_error(
    dynamite(obs_det, test_data, "group", "time", debug = debug), NA
  )
})

test_that("deterministic simultaneity is supported", {
  expect_error(
    obs(y5 ~ x1 + lag(d, 1) + lag(y5, 1) + lag(x1, 1), family = "negbin") +
      aux(numeric(d) ~ y5 + 3),
    NA
  )
})

test_that("deterministic types are supported", {
  expect_error(
    aux(factor(a) ~ factor(c(1, 2, 3), levels = c(1, 2, 3))) +
      aux(numeric(b) ~ log(1.0)) +
      aux(integer(c) ~ 1L) +
      aux(logical(d) ~ TRUE),
    NA
  )
})

test_that("deterministic lags with zero observed lags is evaluated", {
  obs_zerolag <-
    obs(y2 ~ x1, family = "gaussian") +
    aux(numeric(d) ~ abs(y2) + lag(d) | init(0.5))
  expect_error(
    dynamite(obs_zerolag, test_data, "group", "time", debug = debug), NA
  )
})

test_that("past definition computed from data is supported", {
  expect_error(
    obs_past <- obs(y7 ~ lag(d) + lag(y7, 1), family = "exponential") +
      aux(numeric(d) ~ lag(d, 1) + lag(y3, 1) | past(log(abs(x2)))),
    NA
  )
  expect_error(
    fit <- dynamite(obs_past, test_data, "group", "time", debug = debug), NA
  )
})

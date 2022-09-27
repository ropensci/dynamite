#' @srrstats {G5.10} Extended tests can be switched on via setting the
#'   environment variable DYNAMITE_EXTENDED_TESTS to "true".
#' @srrstats {G5.5, G5.6b} Seeds are used appropriately in the tests.
#' @srrstats {G5.4, G5.4a, G5.4b, G5.4c, G5.6, G5.6a, BS7.0, BS7.1, BS7.2}
#'   Simple linear regression and GLM models are tested so that they match with
#'   lm and glm function outputs (within a tolerance due to MCMC, use of
#'   default priors, and discrepancy between ML estimate vs posterior mean).
#'   Further recovery and correctness tests are also implemented.
#' @srrstats {G5.7} Tested that the parameters of the true data generating
#'   process are recovered when increasing the data size.
run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

test_that("parameters for the linear regression are recovered as with lm", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  y <- 2 - 1 * x + rnorm(n, sd = 0.1)
  d <- data.frame(time = 1:n, y = y, x = x)

  fit_lm <- lm(y ~ x, data = d)
  priors <- get_priors(obs(y ~ x, family = "gaussian"),
    data = d, time = "time"
  )
  priors$prior <- c("normal(0, 5)", "normal(0, 1)", "exponential(1)")
  fit_dynamite <- dynamite(obs(y ~ x, family = "gaussian"),
    data = d, time = "time", priors = priors, chains = 1, refresh = 0
  )
  expect_equal(coef(fit_dynamite)$mean, coef(fit_lm),
    tolerance = 0.01,
    ignore_attr = TRUE
  )
})

test_that("parameters for the poisson glm are recovered as with glm", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  y <- rpois(n, exp(2 - 1 * x))
  d <- data.frame(time = 1:n, y = y, x = x)

  fit_glm <- glm(y ~ x, data = d, family = poisson)
  fit_dynamite <- dynamite(obs(y ~ x, family = "poisson"),
    data = d, time = "time", chains = 1, refresh = 0
  )
  expect_equal(coef(fit_dynamite)$mean, coef(fit_glm),
    tolerance = 0.01,
    ignore_attr = TRUE
  )
})

test_that("parameters for the binomial glm are recovered as with glm", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  n <- 100
  u <- sample(1:10, n, TRUE)
  x <- rnorm(n)
  y <- rbinom(n, u, plogis(1 - x))
  d <- data.frame(time = 1:n, y = y, x = x, u = u)

  fit_glm <- glm(cbind(y, u - y) ~ x, data = d, family = binomial)
  fit_dynamite <- dynamite(obs(y ~ x + trials(u), family = "binomial"),
    data = d, time = "time", chains = 1, refresh = 0
  )
  expect_equal(coef(fit_dynamite)$mean, coef(fit_glm),
    tolerance = 0.01,
    ignore_attr = TRUE
  )
})

test_that("parameters for the gamma glm are recovered as with glm", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  y <- rgamma(n, 2, 2 / exp(1 - 2 * x))
  d <- data.frame(time = 1:n, y = y, x = x)

  fit_glm <- glm(y ~ x, data = d, family = Gamma(link = "log"))
  fit_dynamite <- dynamite(obs(y ~ x, family = "gamma"),
    data = d, time = "time", chains = 1, refresh = 0
  )
  expect_equal(coef(fit_dynamite)$mean[1:2], coef(fit_glm),
    tolerance = 0.01,
    ignore_attr = TRUE
  )
})

test_that("parameters for an AR(1) model are recovered as with arima", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  fit <- dynamite(obs(LakeHuron ~ 1, "gaussian") + lags(),
    data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
    "id", "time", chains = 1, refresh = 0
  )
  fit_arima <- arima(LakeHuron, c(1, 0, 0))
  expect_equal(coef(fit)$mean[2], coef(fit_arima)[1],
    tolerance = 0.01,
    ignore_attr = TRUE
  )
  expect_equal(
    coef(fit)$mean[1],
    coef(fit_arima)[2] * (1 - coef(fit_arima)[1]),
    tolerance = 1, ignore_attr = TRUE
  )
})

test_that("parameters of a time-varying gaussian model are recovered", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  create_data <- function(N = 10L, T_ = 100L, D = 50L) {
    K_fixed <- 1L
    K_varying <- 2L
    tau <- c(0.2, 0.4)
    sigma <- 0.1
    beta <- 2.0
    Bs <-
      t(splines::bs(seq.int(1L, T_), df = D, degree = 3L, intercept = TRUE))
    D <- nrow(Bs)
    a <- array(0.0, c(K_varying, D))
    delta <- array(NA, c(T_, K_varying))
    for (k in seq_len(K_varying)) {
      a[k, ] <- cumsum(rnorm(D, 0, tau[k]))
      for (t in seq.int(1L, T_)) {
        delta[t, k] <- a[k, ] %*% Bs[, t]
      }
    }
    x <- matrix(rnorm(T_ * N), N, T_)
    z <- matrix(rbinom(T_ * N, 1.0, 0.7), N, T_)
    y <- matrix(NA, N, T_)
    y[, 1L] <- rnorm(N)
    for (t in seq.int(1L, T_)) {
      m <- beta * z[, t] + delta[t, 1L] +
        delta[t, 2L] * x[, t]
      y[, t] <- rnorm(N, m, sigma)
    }

    list(
      data = data.frame(
        y = c(y), x = c(x), z = c(z), id = seq_len(N),
        time = rep(seq_len(T_), each = N)
      ),
      true_values = c(delta = delta, tau = tau, beta = beta, sigma = sigma)
    )
  }
  d <- create_data()
  dformula <- obs(y ~ -1 + z + varying(~x), family = "gaussian") +
    splines(df = 50)
  # compile model only once
  code <- get_code(dformula,
    data = d$data,
    group = "id",
    time = "time"
  )
  model <- rstan::stan_model(model_code = code)

  # simulate multiple datasets
  n <- 10
  diffs <- matrix(NA, length(d$true_values), n)
  pars <- c("alpha_y", "delta_y", "tau_alpha_y", "tau_y", "beta_y", "sigma_y")
  for (i in seq_len(n)) {
    data <- get_data(dformula, group = "id", time = "time", data = d$data)
    diffs[, i] <- rstan::get_posterior_mean(
      rstan::sampling(model,
        data = data,
        refresh = 0, chains = 1, iter = 2000,
        pars = pars
      ),
      pars = pars
    ) - d$true_values
    d <- create_data()
  }
  # small MSE
  expect_lt(mean(diffs^2), 0.005)

  # test with a single large dataset
  d <- create_data(T_ = 500, N = 500, D = 100)
  data <- get_data(obs(y ~ -1 + z + varying(~x), family = "gaussian") +
    splines(df = 100), group = "id", time = "time", data = d$data)
  fit_long <- rstan::sampling(model,
    data = data,
    refresh = 0, chains = 1, iter = 2000,
    pars = pars
  )
  estimates <- c(rstan::get_posterior_mean(fit_long, pars = pars))
  expect_equal(c(estimates), d$true_values,
    ignore_attr = TRUE, tolerance = 0.1
  )
})


test_that("prior parameters are recovered with zero observations", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  d <- data.frame(y = rep(NA, 10), x = rnorm(10), id = 1, time = 1:10)
  p <- get_priors(obs(y ~ x, "gaussian"), d, "id", "time")
  p$prior[] <- c("normal(2, 0.1)", "normal(5, 0.5)", "exponential(10)")
  fit_prior <- dynamite(obs(y ~ x, "gaussian"),
    data = d,
    group = "id",
    time = "time",
    priors = p,
    iter = 55000,
    warmup = 5000,
    chains = 1,
    cores = 1,
    refresh = 0,
    save_warmup = FALSE
  )
  sumr <- summary(fit_prior) |>
    dplyr::select(parameter, mean, sd, q5, q95) |>
    as.data.frame()

  sigma_y <- sumr |>
    dplyr::filter(parameter == "alpha_y") |>
    dplyr::select(mean, sd, q5, q95)

  m <- 2 - d$x[1] * 5
  s <- sqrt(0.1^2 + d$x[1]^2 * 0.5^2)
  expect_equal(
    unlist(sumr[1, 2:5]),
    c(m, s, qnorm(c(0.05, 0.95), m, s)),
    tolerance = 0.1, ignore_attr = TRUE
  )
  expect_equal(
    unlist(sumr[2, 2:5]),
    c(5, 0.5, qnorm(c(0.05, 0.95), 5, 0.5)),
    tolerance = 0.1, ignore_attr = TRUE
  )
  expect_equal(
    unlist(sumr[3, 2:5]),
    c(0.1, 0.1, qexp(c(0.05, 0.95), 10)),
    tolerance = 0.1, ignore_attr = TRUE
  )
})

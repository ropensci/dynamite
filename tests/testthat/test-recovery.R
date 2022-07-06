#' @srrstats {G5.10} Extended tests can be switched on via setting the
#'   environment variable DYNAMITE_EXTENDED_TESTS to 1.
#' @srrstats {G5.6, G5.6a} Simple linear regression and GLM models are tested
#'   so that they match with lm and glm function (within a tolerance due to
#'   MCMC, use of default priors, and discrepancy between mode vs posterior
#'   mean).
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected
#' to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstats {G5.6b} *Parameter recovery tests should be run with multiple
#' random seeds when either data simulation or the algorithm contains a random component.
#' (When long-running, such tests may be part of an extended, rather than regular, test suite; see G4.10-4.12, below).*
#' @srrstats {G5.4} **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets
#' (potentially through comparisons using binding frameworks such as RStata)*
#' @srrstats {G5.4a} *For new methods, it can be difficult to separate out correctness of the method
#' from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented
#' against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
#' @srrstats {G5.4b} *For new implementations of existing methods, correctness tests should include tests against previous implementations.
#' Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, o
#' r use stored outputs from those where that is not possible.*
#' @srrstats {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original
#' implementations is not available*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.7} **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*


test_that("parameters for the linear regression are recovered as with lm", {
  # Need to skip on CRAN as the compilation can take too much time
  skip_on_cran()
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  y <- 2 - 1 * x + rnorm(n, sd = 0.1)
  d <- data.frame(time = 1:n, y = y, x = x)

  fit_lm <- lm(y ~ x, data = d)
  priors <- get_priors(obs(y ~ x, family = "gaussian"),
    data = d, time = "time")
  priors$prior <- c("normal(0, 5)", "normal(0, 1)", "exponential(1)")
  fit_dynamite <- dynamite(obs(y ~ x, family = "gaussian"),
    data = d, time = "time", priors = priors, chains = 1, refresh = 0)
  expect_equal(coef(fit_lm), coef(fit_dynamite)$mean, tolerance = 0.01,
    ignore_attr = TRUE)
})

test_that("parameters for the poisson glm are recovered as with glm", {
  # Need to skip on CRAN as the compilation can take too much time
  skip_on_cran()
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  y <- rpois(n, exp(2 - 1 * x))
  d <- data.frame(time = 1:n, y = y, x = x)

  fit_glm <- glm(y ~ x, data = d, family = poisson)
  fit_dynamite <- dynamite(obs(y ~ x, family = "poisson"),
    data = d, time = "time", chains = 1, refresh = 0)
  expect_equal(coef(fit_glm), coef(fit_dynamite)$mean, tolerance = 0.01,
    ignore_attr = TRUE)
})

test_that("parameters for the binomial glm are recovered as with glm", {
  # Need to skip on CRAN as the compilation can take too much time
  skip_on_cran()
  set.seed(1)
  n <- 100
  u <- sample(1:10, n, TRUE)
  x <- rnorm(n)
  y <- rbinom(n, u, plogis(1 - x))
  d <- data.frame(time = 1:n, y = y, x = x, u = u)

  fit_glm <- glm(cbind(y, u - y) ~ x, data = d, family = binomial)
  fit_dynamite <- dynamite(obs(y ~ x + trials(u), family = "binomial"),
    data = d, time = "time", chains = 1, refresh = 0)
  expect_equal(coef(fit_glm), coef(fit_dynamite)$mean, tolerance = 0.01,
    ignore_attr = TRUE)
})

test_that("parameters for the gamma glm are recovered as with glm", {
  # Need to skip on CRAN as the compilation can take too much time
  skip_on_cran()
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  y <- rgamma(n, 2, 2/exp(1 - 2 * x))
  d <- data.frame(time = 1:n, y = y, x = x)

  fit_glm <- glm(y ~ x, data = d, family = Gamma(link = "log"))
  fit_dynamite <- dynamite(obs(y ~ x, family = "gamma"),
    data = d, time = "time", chains = 1, refresh = 0)
  expect_equal(coef(fit_glm), coef(fit_dynamite)$mean[1:2], tolerance = 0.01,
    ignore_attr = TRUE)
})


run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "1")

test_that("parameters of a time-varying gaussian model are recovered", {

  skip_if_not(run_extended_tests)

  set.seed(1)
  N <- 100L
  T_ <- 500L
  K_fixed <- 1L
  K_varying <- 2L
  tau <- c(0.2, 0.4)
  sigma <- 0.1
  beta <- 2.0
  Bs <-
    t(splines::bs(seq.int(2L, T_), df = 100L, degree = 3L, intercept = TRUE))
  D <- nrow(Bs)
  sigma_nu <- 0.1
  nu <- rnorm(N, 0.0, sigma_nu)
  a <- array(0.0, c(K_varying, D))
  # Splines start from t = 2, first time point is fixed
  delta <- array(NA, c(T_, K_varying))

  for (k in seq_len(K_varying)) {
    a[k, ] <- cumsum(rnorm(D, 0, tau[k]))
    for (t in seq.int(2L, T_)) {
      delta[t, k] <- a[k, ] %*% Bs[, t - 1]
    }
  }
  x <- matrix(rnorm(T_ * N), N, T_)
  z <- matrix(rbinom(T_ * N, 1.0, 0.7), N, T_)
  y <- matrix(NA, N, T_)
  y[, 1L] <- rnorm(N)
  for (t in seq.int(2L, T_)) {
    m <- nu + beta * z[, t] + delta[t, 1L] +
      delta[t, 2L] * x[, t]
    y[, t] <- rnorm(N, m, sigma)
  }

  d <- data.frame(
    y = c(y), x = c(x), z = c(z), id = seq_len(N),
    time = rep(seq_len(T_), each = N))

  fit <- dynamite(
    dformula =
      obs(
        y ~ -1 + z + varying(~ x),
        family = "gaussian",
        random_intercept = TRUE
      ) +
      splines(df = 100),
    data = d,
    group = "id",
    time = "time",
    iter = 6000,
    warmup = 1000,
    chains = 2,
    cores = 2,
    refresh = 0,
    save_warmup = FALSE
  )
})

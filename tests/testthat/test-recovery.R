#' @srrstats {G5.10} Extended tests can be switched on via setting the
#'   environment variable DYNAMITE_EXTENDED_TESTS to 1.
#' @srrstats {G5.5, G5.6b} Seeds are used appropriately in the tests.
#' @srrstats {G5.6, G5.6a, BS7.0, BS7.1, BS7.2} Simple linear regression and
#'   GLM models are tested so that they match with lm and glm function
#'   (within a tolerance due to MCMC, use of default priors, and discrepancy
#'   between mode vs posterior mean). Furher recovery tests are also
#'   implemented.
#'
#' @srrstats {G5.4} **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets
#' (potentially through comparisons using binding frameworks such as RStata)*
#'
#' @srrstats {G5.4a, G5.4b, G5.4c}
#'
#'  *For new methods, it can be difficult to separate out correctness of the method
#' from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented
#' against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
#' @srrstats {G5.4b} *For new implementations of existing methods, correctness tests should include tests against previous implementations.
#' Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, o
#' r use stored outputs from those where that is not possible.*
#' @srrstats {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original
#' implementations is not available*
#' @srrstats {G5.7} **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*

set.seed(123)
seeds <- sample(1:1000, size = 5)

test_that("parameters for the linear regression are recovered as with lm", {
  # Need to skip on CRAN as the compilation can take too much time
  skip_on_cran()
  set.seed(seed[1])
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
  set.seed(seed[2])
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
  set.seed(seed[3])
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
  set.seed(seed[4])
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

  set.seed(seed[5])

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

    list(data = data.frame(
      y = c(y), x = c(x), z = c(z), id = seq_len(N),
      time = rep(seq_len(T_), each = N)),
      true_values = c(delta, tau, beta, sigma))
  }
  d <- create_data()
  dformula <- obs(y ~ -1 + z + varying(~ x), family = "gaussian") +
    splines(df = 50)
  # compile model only once
  code <- get_code(dformula,
    data = d$data,
    group = "id",
    time = "time")
  model <- rstan::stan_model(model_code = code)

  n <- 10
  diffs <- matrix(NA, length(d$true_values), n)
  pars <- c("alpha_y",  "delta_y", "tau_alpha_y", "tau_y", "beta_y", "sigma_y")
  for (i in seq_len(n)) {
    data <- get_data(dformula, group = "id", time = "time", data = d$data)
    diffs[, i] <- rstan::get_posterior_mean(
      rstan::sampling(model, data = data,
      refresh = 0, chains = 1, iter = 2000,
      pars = pars), pars = pars) - d$true_values
    d <- create_data()
    print(i)
  }

  expect_lt(mean(diffs^2), 0.005)

  d <- create_data(T_ = 1000, N = 1000, D = 500)
  data <- get_data(dformula, group = "id", time = "time", data = d$data)
  estimates <- c(rstan::get_posterior_mean(
    rstan::sampling(model, data = data,
      refresh = 0, chains = 1, iter = 2000,
      pars = pars), pars = pars))
  expect_equal(c(estimates), d$true_values)

  p <- get_priors(dformula,
    data = d$data,
    group = "id",
    time = "time")
  p$priors[] #TODO WIP
  d$data$y[] <- NA
  fit_prior <- dynamite(dformula,
    data = d$data,
    group = "id",
    time = "time",
    priors = p,
    iter = 2000,
    warmup = 1000,
    chains = 1,
    cores = 1,
    refresh = 0,
    save_warmup = FALSE
  )
  # TODO WIP
})

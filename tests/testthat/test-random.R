#' @srrstats {G5.10} Extended tests can be switched on via setting the
#'   environment variable DYNAMITE_EXTENDED_TESTS to "true".

run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

test_that("centered and noncentered parameterization for random effects work", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  n <- 40
  k <- 10
  x <- rnorm(n * k)
  u1 <- rep(rnorm(k, sd = 0.2), each = n)
  u2 <- rep(rnorm(k, sd = 0.1), each = n)
  mu <- exp(4 +x + u1 + u2 * x)
  phi <- 50
  y <- rnbinom(n * k, mu = mu, size = phi)
  hist(y)
  d <- data.frame(year = 1:n, person = rep(1:k, each = n), y = y, x = x)

  fit_centered <- dynamite(
    obs(y ~ x + random(~ 1 + x), family = "negbin") +
      random_spec(noncentered = FALSE, correlated = TRUE),
    data = d, time = "year", group = "person",
    chains = 2, iter = 4000, refresh = 0, seed = 1
  )
  fit_noncentered <- dynamite(
    obs(y ~ x + random(~ 1 + x), family = "negbin") +
      random_spec(noncentered = TRUE, correlated = TRUE),
    data = d, time = "year", group = "person",
    chains = 2, iter = 4000, refresh = 0, seed = 1
  )

  expect_equal(
    summary(fit_centered,
      types = c("alpha", "beta", "corr_nu", "sigma_nu", "nu"))$mean,
    summary(fit_noncentered,
      types = c("alpha", "beta", "corr_nu", "sigma_nu", "nu"))$mean,
    tolerance = 0.1
  )
  expect_equal(
    summary(fit_centered, parameter = "phi_y")$mean,
    summary(fit_noncentered, parameter = "phi_y")$mean,
    tolerance = 0.2
  )

  fit_centered_nocorr <- dynamite(
    obs(y ~ x + random(~ 1 + x), family = "negbin") +
      random_spec(noncentered = FALSE, correlated = FALSE),
    data = d, time = "year", group = "person",
    chains = 2, iter = 4000, refresh = 0, seed = 1
  )
  fit_noncentered_nocorr <- dynamite(
    obs(y ~ x + random(~ 1 + x), family = "negbin") +
      random_spec(noncentered = TRUE, correlated = FALSE),
    data = d, time = "year", group = "person",
    chains = 2, iter = 4000, refresh = 0, seed = 1
  )

  expect_equal(
    summary(fit_centered_nocorr,
      types = c("alpha", "beta", "sigma_nu", "nu"))$mean,
    summary(fit_noncentered_nocorr,
      types = c("alpha", "beta","sigma_nu", "nu"))$mean,
    tolerance = 0.1
  )
  expect_equal(
    summary(fit_noncentered,
      types = c("alpha", "beta", "sigma_nu", "nu"))$mean,
    summary(fit_noncentered_nocorr,
      types = c("alpha", "beta","sigma_nu", "nu"))$mean,
    tolerance = 0.1
  )
})

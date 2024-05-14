data.table::setDTthreads(1) # For CRAN

test_that("priors can be extracted", {
  expect_error(
    get_priors(gaussian_example_fit),
    NA
  )
})

f <- obs(y ~ -1 + random(~1) + z + varying(~ x + lag(y)), family = "gaussian") +
  random_spec() + splines(df = 20)
p <- get_priors(gaussian_example_fit)

test_that("manual prior setting works", {
  expect_error(
    fit <- dynamite(
      f,
      data = gaussian_example, time = "time", group = "id",
      priors = p, debug = list(no_compile = TRUE)
    ),
    NA
  )
})

test_that("extracted priors match initial priors", {
  fit <- dynamite(
    f,
    data = gaussian_example, time = "time", group = "id",
    priors = p, debug = list(no_compile = TRUE)
  )
  p <- get_priors(gaussian_example_fit)
  expect_identical(get_priors(fit), p)
})

test_that("inserting a valid prior works", {
  p$prior[2] <- "cauchy(0, 2)"
  p$prior[5:6] <- "std_normal()"
  expect_error(
    dynamite(f,
      data = gaussian_example, time = "time", group = "id",
      priors = p, debug = list(no_compile = TRUE)
    ),
    NA
  )
})

test_that("manual prior setting works", {
  testdata <- data.frame(
    y = c(0, rexp(9, 1)),
    x = c(NA, rbeta(9, 2, 2)),
    z = c(0, rnbinom(9, 5, 0.5)),
    w = c(0, 3 + 2 * rt(9, 3)),
    t = 1:10
  )
  f <-  obs(y ~ x, "gamma") +
    obs(x ~ z, "beta") +
    obs(z ~ 1, "negbin") +
    obs(w ~ 1, "student")
  expect_error(
    p <- get_priors(f, data = testdata, time = "t"),
    NA
  )
  expect_identical(
    p$parameter,
    c("alpha_y", "beta_y_x", "phi_y", "alpha_x", "beta_x_z", "phi_x",
      "alpha_z", "phi_z", "alpha_w", "sigma_w", "phi_w")
  )
  expect_error(
    fit <- dynamite(
      f, data = testdata, time = "t", priors = p,
      debug = list(no_compile = TRUE)
    ),
    NA
  )
  expect_identical(get_priors(fit), p)
})

test_that("manual priors for multivariate gaussian channel works", {
  y <- rnorm(10)
  x <- rexp(10)
  testdata <- data.frame(
    x = x,
    y1 = y + 0.5 * x,
    y2 = 0.25 * y + 1.5 * x + rnorm(10),
    t = 1:10
  )
  f <- obs(c(y1, y2) ~ x, family = "mvgaussian")
  expect_error(
    p <- get_priors(f, data = testdata, time = "t"),
    NA
  )
  expect_identical(
    p$parameter,
    c(
      "alpha_y1", "beta_y1_x", "sigma_y1", "alpha_y2", "beta_y2_x", "sigma_y2",
      "L_y1_y2"
    )
  )
  expect_error(
    fit <- dynamite(
      f, data = testdata, time = "t", priors = p,
      debug = list(no_compile = TRUE)
    ),
    NA
  )
  expect_identical(get_priors(fit), p)
})

test_that("manual priors for multinomial channel works", {
  x <- rnorm(10)
  y1 <- sample(5, 10, replace = TRUE)
  y2 <- sample(6, 10, replace = TRUE)
  n <- y1 + y2
  testdata <- data.frame(y1 = y1, y2 = y2, x = x, n = n, t = 1:10)
  f <- obs(c(y1, y2) ~ -1 + varying(~ x) + trials(n), family = "multinomial") +
    splines(df = 10)
  expect_error(
    p <- get_priors(f, data = testdata, time = "t"),
    NA
  )
  expect_identical(
    p$parameter,
    c("alpha_y2", "tau_alpha_y2", "delta_y2_x", "tau_y2_x")
  )
  expect_error(
    fit <- dynamite(
      f, data = testdata, time = "t", priors = p,
      debug = list(no_compile = TRUE)
    ),
    NA
  )
  expect_identical(get_priors(fit), p)
})

test_that("manual priors for cumulative channel works", {
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
    y[, i] <- apply(p, 1, sample, x = letters[1:4], size = 1, replace = FALSE)
  }

  d <- data.frame(
    y = factor(c(y)), x = c(x),
    time = rep(seq_len(t), each = n),
    id = rep(seq_len(n), t)
  )
  f <- obs(y ~ x, family = "cumulative", link = "logit")

  expect_error(
    p <- get_priors(
      f,
      data = d,
      time = "time",
      group = "id"
    ),
    NA
  )
  expect_identical(
    p$parameter,
    c("cutpoint_y_1", "cutpoint_y_2", "cutpoint_y_3", "beta_y_x")
  )
  expect_error(
    fit <- dynamite(
      f,
      data = d,
      time = "time",
      group = "id",
      priors = p,
      debug = list(no_compile = TRUE)
    ),
    NA
  )
  expect_identical(get_priors(fit), p)

  f <- obs(y ~ -1 + x + varying(~ 1), family = "cumulative", link = "probit") +
    splines()

  expect_error(
    p <- get_priors(
      f,
      data = d,
      time = "time",
      group = "id"
    ),
    NA
  )
  expect_identical(
    p$parameter,
    c("alpha_y_1", "alpha_y_2", "alpha_y_3", "tau_alpha_y_1", "tau_alpha_y_2",
      "tau_alpha_y_3", "beta_y_x")
  )
  expect_error(
    fit <- dynamite(
      f,
      data = d,
      time = "time",
      group = "id",
      priors = p,
      debug = list(no_compile = TRUE)
    ),
    NA
  )
  expect_identical(get_priors(fit), p)
})

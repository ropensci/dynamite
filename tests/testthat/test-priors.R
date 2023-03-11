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

test_that("prior order does not matter", {
  p2 <- p[seq.int(nrow(p), 1L), ]
  fit2 <- dynamite(f,
    data = gaussian_example, time = "time", group = "id",
    priors = p2, debug = list(no_compile = TRUE)
  )
  expect_identical(get_priors(fit2), p)
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
    c("alpha_y1_y2", "tau_alpha_y1_y2", "delta_y1_y2_x", "tau_y1_y2_x")
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

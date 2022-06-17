test_that("priors can be extracted", {
  expect_error(
    p <- get_priors(gaussian_example_fit),
    NA
  )

  f <- obs(y ~ -1 + z + varying(~ x + lag(y)), family = gaussian(),
           random_intercept = TRUE) + splines(df = 20)

  expect_error(
    fit <- dynamite(
      f, data = gaussian_example, time = "time", group = "id",
      priors = p, debug = list(no_compile = TRUE)),
    NA
  )

  expect_identical(get_priors(fit), p)

  p2 <- p[nrow(p):1,]
  fit2 <- dynamite(f, data = gaussian_example, time = "time", group = "id",
                   priors = p2, debug = list(no_compile = TRUE))

  expect_identical(get_priors(fit2), p)

  p2$prior[2] <- "cauchy(0, 2)"
  expect_error(
    dynamite(f, data = gaussian_example, time = "time", group = "id",
             priors = p2, debug = list(no_compile = TRUE)),
    NA
  )

  p2$prior[5] <- "aaa"
  expect_error(
    dynamite(f, data = gaussian_example, time = "time", group = "id",
             priors = p2, debug = list(no_compile = TRUE))
  )

  p2$prior[5] <- "gamma(2,1)"
  expect_error(
    dynamite(f, data = gaussian_example, time = "time", group = "id",
             priors = p2, debug = list(no_compile = TRUE))
  )
})

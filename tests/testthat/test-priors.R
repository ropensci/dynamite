test_that("priors can be extracted", {
  expect_error(
    get_priors(gaussian_example_fit),
    NA
  )
})

f <- obs(y ~ -1 + z + varying(~ x + lag(y)), family = "gaussian") +
  random() + splines(df = 20)
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

fit <- dynamite(
  f,
  data = gaussian_example, time = "time", group = "id",
  priors = p, debug = list(no_compile = TRUE)
)

test_that("extracted priors match initial priors", {
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
  expect_error(
    dynamite(f,
      data = gaussian_example, time = "time", group = "id",
      priors = p, debug = list(no_compile = TRUE)
    ),
    NA
  )
})

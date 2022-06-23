
# Data warnings -----------------------------------------------------------

test_that("ordered factor conversion to factor warns", {
  test_data <- data.frame(y = factor(c(1, 2), ordered = TRUE),
                          x = c(1, 1), z = c(1, 2))
  expect_warning(
    dynamite(dformula = obs(y ~ x, family = "categorical"),
             data = test_data, group = "x", time = "z",
             debug = list(no_compile = TRUE)),
    "Response variable `y` is of class <ordered factor> whose channel is categorical:\ni `y` will be converted to an unordered factor\\."
  )
})

test_that("perfect collinearity warns", {
  f <- obs(y ~ x + z, family = "gaussian")
  test_data1 <- data.frame(y = rnorm(10), x = rep(1, 10), z = rep(2, 10))
  test_data2 <- data.frame(y = rep(1, 10), x = rep(1, 10), z = rnorm(10))
  expect_warning(
    full_model.matrix(f, test_data1),
    "Perfect collinearity found between predictor variables of channel `y`\\."
  )
  expect_warning(
    full_model.matrix(f, test_data2),
    "Perfect collinearity found between response and predictor variable:\ni Response variable `y` is perfectly collinear with predictor variable `x`\\."
  )
})

# Specials warnings -------------------------------------------------------

test_that("multiple intercept warns", {
  expect_warning(
    obs(y ~ 1 + varying(~1), family = "gaussian"),
    "Both time-independent and time-varying intercept specified:\ni Defaulting to time-varying intercept\\."
  )
})

test_that("deterministic fixed warns", {
  expect_warning(
    aux(numeric(y) ~ fixed(~ x)),
    "fixed\\(\\) definitions of a determinstic channel `numeric\\(y\\)` will be ignored\\."
  )
})

test_that("deterministic varying warns", {
  expect_warning(
    aux(numeric(y) ~ varying(~ x)),
    "varying\\(\\) definitions of a determinstic channel `numeric\\(y\\)` will be ignored\\."
  )
})

test_that("untyped deterministic warns", {
  expect_warning(
    aux(y ~ 1 + x),
    "No type specified for deterministic channel `y`:\ni Assuming type is <numeric>\\."
  )
})


# Predict warnings --------------------------------------------------------

test_that("Too many n_draws warns", {
  expect_warning(
    predict(gaussian_example_fit, n_draws = 1e6),
    "You've supplied `n_draws` = 1000000 but there are only 400 samples available:\ni The available samples will be used for prediction\\."
  )
})

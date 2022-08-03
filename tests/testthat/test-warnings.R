# Data warnings -----------------------------------------------------------

test_that("ordered factor conversion to factor warns", {
  test_data <- data.frame(y = factor(c(1, 2, 2), ordered = TRUE),
                          x = c(1, 1, 2), z = c(1, 2, 3))
  expect_warning(
    dynamite(dformula = obs(y ~ x, family = "categorical"),
             data = test_data, group = "x", time = "z",
             debug = list(no_compile = TRUE)),
    paste0(
      "Response variable `y` is of class <ordered factor> ",
      "whose channel is categorical:\n",
      "i `y` will be converted to an unordered factor\\."
    )
  )
})

test_that("perfect collinearity warns", {
  f1 <- obs(y ~ -1 + x + z, family = "gaussian")
  f2 <- obs(y ~ z, family = "gaussian")
  test_data1 <- data.frame(y = rnorm(10), x = rep(1, 10), z = rep(2, 10))
  test_data2 <- data.frame(y = rep(1, 10), x = rep(1, 10), z = rnorm(10))
  expect_warning(
    full_model.matrix(f1, test_data1, TRUE),
    "Perfect collinearity found between predictor variables of channel `y`\\."
  )
  expect_warning(
    full_model.matrix(f2, test_data2, TRUE),
    paste0(
      "Perfect collinearity found between response and predictor variable:\n",
      "i Response variable `y` is perfectly collinear ",
      "with predictor variable `\\(Intercept\\)`\\."
    )
  )
  expect_warning(
    full_model.matrix(f1, test_data2, TRUE),
    paste0(
      "Perfect collinearity found between response and predictor variable:\n",
      "i Response variable `y` is perfectly collinear ",
      "with predictor variable `x`\\."
    )
  )
})

test_that("too few observations warns", {
  f <- obs(y ~ x + z + w, family = "gaussian")
  test_data <- data.frame(y = rnorm(3), x = rnorm(3), z = rnorm(3),
    w = rnorm(3))
  expect_warning(
    full_model.matrix(f, test_data, TRUE),
    paste0(
      "Number of non-missing observations 3 in channel `y` ",
      "is less than 4, the number of predictors \\(including possible ",
      "intercept\\)\\."
    )
  )
})

test_that("zero predictor warns", {
  f <- obs(y ~ -1 + x + z, family = "gaussian")
  test_data <- data.frame(
    y = rnorm(6),
    x = c(NA, rnorm(2), NA, rnorm(2)),
    z = factor(1:3))
  expect_warning(
    full_model.matrix(f, test_data, TRUE),
    paste0(
      "Predictor `z1` contains only zeros in the complete case rows of the ",
      "design matrix for the channel `y`\\."
    )
  )
})

# Specials warnings -------------------------------------------------------

test_that("multiple intercept warns", {
  expect_warning(
    obs(y ~ 1 + varying(~1), family = "gaussian"),
    paste0(
      "Both time-independent and time-varying intercept specified:\n",
      "i Defaulting to time-varying intercept\\."
    )
  )
})

test_that("untyped deterministic warns", {
  expect_warning(
    aux(y ~ 1 + x),
    paste0(
      "No type specified for deterministic channel `y`:\n",
      "i Assuming type is <numeric>\\."
    )
  )
})

# Predict warnings --------------------------------------------------------

test_that("too large n_draws warns", {
  expect_warning(
    predict(gaussian_example_fit, n_draws = 500),
    paste0(
      "You've supplied `n_draws` = 500 but ",
      "there are only 400 samples available:\n",
      "i The available samples will be used for prediction\\."
    )
  )
})

test_that("gaps in newdata with exogenous predictors and no impute warns", {
  newdata <- multichannel_example |>
    dplyr::mutate(b = ifelse(time > 5, NA, b)) |>
    dplyr::filter(time < 3 | time > 10)
  expect_warning(
    predict(multichannel_example_fit, newdata = newdata, n_draws = 4),
    paste0("Time index variable `time` of `newdata` has gaps:\n",
     "i Filling the `newdata` to regular time points\\. This will lead to ",
     "propagation of NA values if the model contains exogenous predictors ",
     "and `impute` is \"none\"\\."
    )
  )
})

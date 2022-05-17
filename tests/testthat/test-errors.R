obs_test <- obs(y ~ x, family = gaussian())


# Formula errors ----------------------------------------------------------

test_that("nonformula to dynamiteformula fails", {
  expect_error(
    obs(formula = numeric()),
    "Argument 'formula' is not a formula object"
  )
})

test_that("unsupported family fails", {
  expect_error(
    obs(y ~ x, family = "unknown_distr"),
    "Family 'unknown_distr' is not supported"
  )
})

test_that("unrecognized family call fails", {
  myfamily <- function() invisible(NULL)
  expect_error(
    obs(y ~ x, family = myfamily()),
    "Unsupported family call 'myfamily\\(\\)'"
  )
})

test_that("duplicate response definition fails", {
  expect_error(
    obs_test + obs_test,
    "Multiple definitions for response variable\\(s\\): y"
  )
})

test_that("duplicate spline definition fails", {
  obs_lhs <- obs_test + splines()
  obs_rhs <- obs(z ~ x, family = gaussian()) + splines()
  expect_error(
    obs_lhs + obs_rhs,
    "Multiple definitions for splines"
  )
})

test_that("duplicate lags definition fails", {
  obs_lhs <- obs_test + lags(k = 1)
  obs_rhs <- obs(z ~ x, family = gaussian()) + lags(k = 2)
  expect_error(
    obs_lhs + obs_rhs,
    "Multiple definitions for lags"
  )
})

test_that("simultaneity fails", {
  expect_error(
    obs_test + obs(x ~ y, family = "gaussian"),
    "Simultaneous regression is not supported, response variables 'y' appear in the formulas of 'x'"
  )
})

# Data errors -------------------------------------------------------------

test_that("data is not data.frame fails", {
  expect_error(
    dynamite(dformula = obs_test, data = list()),
    "Argument 'data' is not a data.frame object")
})

test_that("group variable not in data fails", {
  expect_error(
    dynamite(dformula = obs_test, data = data.frame(y = 1, x = 1), group = "z"),
    "Grouping variable 'z' is not present in the data"
  )
})

test_that("missing time variable fails", {
  expect_error(
    dynamite(dformula = obs_test, data = data.frame(z = 1), group = "z"),
    "Argument 'time' is missing"
  )
})

test_that("time variable not in data fails", {
  expect_error(
    dynamite(dformula = obs_test, data = data.frame(y = 1, x = 1), time = "z"),
    "Time index variable 'z' is not present in the data"
  )
})



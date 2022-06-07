obs_test <- obs(y ~ x + w, family = gaussian())

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

test_that("as-is use fails", {
  expect_error(
    obs(y ~ I(x), family = gaussian()),
    "The use of I\\(.\\) is not supported by dynamiteformula"
  )
})

test_that("duplicate response definition fails", {
  expect_error(
    obs_test + obs_test,
    "Multiple definitions for response variable\\(s\\): y"
  )
})

test_that("duplicate spline definition fails", {
  expect_error(
    obs_test + splines() + splines(),
    "Multiple definitions for splines"
  )
})

test_that("duplicate algs definition fails", {
  expect_error(
    obs_test + lags() + lags(),
    "Multiple definitions for lags"
  )
})

test_that("attempting to add dynamiteformulas with lag definitions fails", {
  obs_lhs <- obs_test + lags(k = 1)
  obs_rhs <- obs(z ~ x, family = gaussian()) + lags(k = 2)
  expect_error(
    obs_lhs + obs_rhs,
    "Both dynamiteformulas contain a lags definition"
  )
})

test_that("attempting to add dynamiteformulas with splines definitions fails", {
  obs_lhs <- obs_test + splines()
  obs_rhs <- obs(z ~ x, family = gaussian()) + splines()
  expect_error(
    obs_lhs + obs_rhs,
    "Both dynamiteformulas contain a splines definition"
  )
})

test_that("simultaneity fails", {
  obs_lhs <-
    obs(q ~ w + e + r + lag(i), family = gaussian()) +
    obs(t ~ y + u, family = gaussian()) +
    obs(i ~ o + p + a + lag(f), family = gaussian())
  obs_rhs <-
    obs(f ~ h + l + lag(x), family = gaussian()) +
    obs(x ~ q + z, family = gaussian())
  expect_error(
    obs_rhs + obs_lhs,
    "Simultaneous regression is not supported, response variable 'q' appears in the formula of 'x'"
  )
  # should fail for deterministic as well
  expect_error(
    obs(y ~ x, family = gaussian()) + aux(integer(x) ~ y),
    "Simultaneous regression is not supported, response variable 'x' appears in the formula of 'y'"
  )
})

test_that("adding nondynamiteformula to dynamiteformula fails", {
  expect_error(
    obs_test + 1.0,
    "Unable to add an object of class 'numeric' to an object of class 'dynamiteformula'"
  )
})

test_that("plus method fails for nondynamiteformula", {
  expect_error(
    `+.dynamiteformula`(data.frame(), numeric()),
    "Method '\\+\\.dynamiteformula' is not supported for 'data.frame' objects"
  )
})

test_that("categorical random intercept fails", {
  expect_error(
    obs(y ~ x, family = categorical(), random_intercept = TRUE),
    "Random intercepts are not yet supported for the categorical family"
  )
})

test_that("negative lb_tau fails", {
  expect_error(
    obs_test + splines(lb_tau = -1.0),
    "Lower bound for 'tau' should be non-negative"
  )
})

# Formula specials errors -------------------------------------------------

test_that("Specification as both fixed and varying fails", {
  expect_error(
    obs(y ~ x + varying(~ x), family = gaussian()),
    "Variables 'x' specified as both time-constant and time-varying"
  )
})

test_that("No intercept or predictors fails", {
  expect_error(
    obs(y ~ -1, family = gaussian()),
    "Invalid formula for response variable 'y', there are no predictors nor an intercept"
  )
})

test_that("past in the middle of formula fails", {
  expect_error(
    aux(y ~ x + past(0) + z),
    "Past values term must be the last term of the formula"
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

test_that("single time point fails", {
  expect_error(
    dynamite(dformula = obs_test, data = data.frame(y = 1, x = 1, z = 1),
             group = "x", time = "z"),
    "There must be at least two time points in the data"
  )
})

test_that("negative lag fails", {
  expect_error(
    dynamite(dformula = obs(y ~ lag(y, -1), family = gaussian()),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x", time = "z"),
    "Only positive shift values are allowed in lag()"
  )
})

test_that("missing lag variable fails", {
  expect_error(
    dynamite(dformula = obs(y ~ lag(d, 1), family = gaussian()),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x", time = "z"),
    "Unable to construct lagged values of 'd', no such variables are present in the data"
  )
})

test_that("missing predictor fails", {
  expect_error(
    dynamite(dformula = obs(y ~ w, family = gaussian()),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x", time = "z"),
    "Variables 'w' used in the formula are not present in the data"
  )
})

test_that("invalid deterministic channel definition fails", {
  expect_error(
    dynamite(dformula = aux(integer(d) ~ 1 + w),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x", time = "z"),
    "Unable to evaluate the definitions of deterministic channels.+"
  )
})

test_that("irregular time intervals fails", {
  data_irreg <- data.frame(
    y = c(1, 2, 3, 4, 5),
    x = c(1, 1, 1, 2, 2),
    t = c(2, 5, 7, 3.5, 5.75)
  )
  expect_error(
    dynamite(obs_test, data = data_irreg, group = "x", time = "t"),
    "Observations must occur at regular time intervals"
  )
})

test_that("deterministic insufficient initial values fails", {
  expect_error(
    dynamite(dformula = aux(numeric(d) ~ lag(d, 1)),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x",
             time = "z"),
    "Deterministic channel 'd' requires 1 initial values, but only 0 values have been specified"
  )
})

# Data type errors --------------------------------------------------------

test_that("factor types for non-categorical families fails", {
  test_data <- data.frame(y = factor(c(0, 1)), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = gaussian()),
             data = test_data, group = "x", time = "z"),
    "Response variable 'y' is invalid: gaussian family is not supported for factors"
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = bernoulli()),
             data = test_data, group = "x", time = "z"),
    "Response variable 'y' is invalid: bernoulli family is not supported for factors"
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = binomial()),
             data = test_data, group = "x", time = "z"),
    "Response variable 'y' is invalid: binomial family is not supported for factors"
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = poisson()),
             data = test_data, group = "x", time = "z"),
    "Response variable 'y' is invalid: Poisson family is not supported for factors"
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = negbin()),
             data = test_data, group = "x", time = "z"),
    "Response variable 'y' is invalid: negative binomial family is not supported for factors"
  )
})

test_that("negative values for binomial, negbin and poisson fails", {
  test_data <- data.frame(y = c(-1, -2), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = binomial()),
             data = test_data, group = "x", time = "z"),
    "Response variable 'y' is invalid: binomial family supports only non-negative integers"
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = poisson()),
             data = test_data, group = "x", time = "z"),
    "Response variable 'y' is invalid: Poisson family supports only non-negative integers"
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = negbin()),
             data = test_data, group = "x", time = "z"),
    "Response variable 'y' is invalid: negative binomial family supports only non-negative integers"
  )
})

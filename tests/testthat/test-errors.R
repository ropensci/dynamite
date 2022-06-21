obs_test <- obs(y ~ x + w, family = gaussian())

# Formula errors ----------------------------------------------------------

test_that("nonformula to dynamiteformula fails", {
  expect_error(
    obs(formula = numeric()),
    "Argument `formula` is not a <formula> object\\."
  )
})

test_that("unsupported family fails", {
  expect_error(
    obs(y ~ x, family = "unknown_distr"),
    'Family "unknown_distr" is not supported\\.'
  )
})

test_that("unrecognized family call fails", {
  myfamily <- function() invisible(NULL)
  expect_error(
    obs(y ~ x, family = myfamily()),
    "Unsupported family call `myfamily\\(\\)`\\."
  )
})

test_that("as-is use fails", {
  expect_error(
    obs(y ~ I(x), family = gaussian()),
    "The use of `I\\(\\.\\)` is not supported by `dynamiteformula\\(\\)`\\."
  )
})

test_that("duplicate response definition fails", {
  expect_error(
    obs_test + obs_test,
    "Multiple definitions for response variable `y`\\."
  )
})

test_that("duplicate spline definition fails", {
  expect_error(
    obs_test + splines() + splines(),
    "Multiple definitions for splines\\."
  )
})

test_that("duplicate algs definition fails", {
  expect_error(
    obs_test + lags() + lags(),
    "Multiple definitions for lags\\."
  )
})

test_that("attempting to add dynamiteformulas with lag definitions fails", {
  obs_lhs <- obs_test + lags(k = 1)
  obs_rhs <- obs(z ~ x, family = gaussian()) + lags(k = 2)
  expect_error(
    obs_lhs + obs_rhs,
    "Both dynamiteformulas contain a lags definition\\."
  )
})

test_that("attempting to add dynamiteformulas with splines definitions fails", {
  obs_lhs <- obs_test + splines()
  obs_rhs <- obs(z ~ x, family = gaussian()) + splines()
  expect_error(
    obs_lhs + obs_rhs,
    "Both dynamiteformulas contain a splines definition\\."
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
    "Simultaneous regression is not supported:\nx Response variable `q` appears in the formula of `x`\\."
  )
  # should fail for deterministic as well
  expect_error(
    obs(y ~ x, family = gaussian()) + aux(integer(x) ~ y),
    "Simultaneous regression is not supported:\nx Response variable `x` appears in the formula of `y`\\."
  )
})

test_that("adding nondynamiteformula to dynamiteformula fails", {
  expect_error(
    obs_test + 1.0,
    "Unable to add an object of class <numeric> to an object of class <dynamiteformula>\\."
  )
})

test_that("plus method fails for nondynamiteformula", {
  expect_error(
    `+.dynamiteformula`(data.frame(), numeric()),
    "Method `\\+\\.dynamiteformula\\(\\)` is not supported for <data.frame> objects\\."
  )
})

test_that("categorical random intercept fails", {
  expect_error(
    obs(y ~ x, family = categorical(), random_intercept = TRUE),
    "Random intercepts are not yet supported for the categorical family\\."
  )
})

test_that("negative lb_tau fails", {
  expect_error(
    obs_test + splines(lb_tau = -1.0),
    "Lower bound for `tau` must be non-negative\\."
  )
})

# Formula specials errors -------------------------------------------------

test_that("no intercept or predictors fails", {
  expect_error(
    obs(y ~ -1, family = gaussian()),
    "Invalid formula for response variable `y`:\nx There are no predictors nor an intercept term\\."
  )
})

test_that("past in the middle of formula fails", {
  expect_error(
    aux(y ~ x + past(0) + z),
    "Past values term must be the last term of the formula\\."
  )
})

# Data errors -------------------------------------------------------------

test_that("data is not data.frame fails", {
  expect_error(
    dynamite(dformula = obs_test, data = list()),
    "Argument `data` is not a <data.frame> object\\.")
})

test_that("group variable not in data fails", {
  expect_error(
    dynamite(dformula = obs_test, data = data.frame(y = 1, x = 1), group = "z"),
    "Can't find grouping variable `z` in `data`\\."
  )
})

test_that("missing time variable fails", {
  expect_error(
    dynamite(dformula = obs_test, data = data.frame(z = 1), group = "z"),
    "Argument `time` is missing\\."
  )
})

test_that("time variable not in data fails", {
  expect_error(
    dynamite(dformula = obs_test, data = data.frame(y = 1, x = 1), time = "z"),
    "Can't find time index variable `z` in `data`\\."
  )
})

test_that("single time point fails", {
  expect_error(
    dynamite(dformula = obs_test, data = data.frame(y = 1, x = 1, z = 1),
             group = "x", time = "z"),
    "There must be at least two time points in the data."
  )
})

test_that("negative lag fails", {
  expect_error(
    dynamite(dformula = obs(y ~ lag(y, -1), family = gaussian()),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x", time = "z"),
    "Shift values must be positive in `lag\\(\\)`\\."
  )
})

test_that("missing lag variable fails", {
  expect_error(
    dynamite(dformula = obs(y ~ lag(d, 1), family = gaussian()),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x", time = "z"),
    "Unable to construct lagged values of `d`:\nx Can't find such variables in `data`\\."
  )
})

test_that("missing predictor fails", {
  expect_error(
    dynamite(dformula = obs(y ~ w, family = gaussian()),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x", time = "z"),
    "Can't find variable `w` in `data`\\."
  )
})

test_that("invalid deterministic channel definition fails", {
  expect_error(
    dynamite(dformula = aux(integer(d) ~ 1 + w),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x", time = "z"),
    "Unable to evaluate definitions of deterministic channels:\ni Some variables are possibly missing or incorrect\\."
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
    "Observations must occur at regular time intervals\\."
  )
})

test_that("deterministic insufficient initial values fails", {
  expect_error(
    dynamite(dformula = aux(numeric(d) ~ lag(d, 1)),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x",
             time = "z"),
    "Deterministic channel `d` requires 1 initial value:\nx You've supplied no initial values\\."
  )
})

# Data type errors --------------------------------------------------------

#' @srrstats {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
#' @srrstats {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.*
test_that("invalid column types fail", {
  test_data <- data.frame(y = c(1i, 2i), x = c(1, 1), z = c(1, 2))
  test_data$w <- c(list(a = 1), list(b = 2))
  test_data$d <- as.raw(c(40, 20))
  expect_error(
    dynamite(dformula = obs(y ~ x, family = gaussian()),
             data = test_data, group = "x", time = "z"),
    "Columns `y`, `w`, and `d` of `data` are invalid:\nx Column types <complex/list/raw> are not supported\\."
  )
})

test_that("non-finite values in data fail", {
  test_data <- data.frame(y = c(1, Inf), x = c(1, 1),
                          z = c(1, 2), w = c(-Inf, 2), u = c(1, Inf))
  expect_error(
    dynamite(dformula = obs(y ~ x, family = gaussian()),
             data = test_data, group = "x", time = "z"),
    "Non-finite values in variables `y`, `w`, and `u` of `data`\\."
  )
})

test_that("factor types for non-categorical families fails", {
  test_data <- data.frame(y = factor(c(0, 1)), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = gaussian()),
             data = test_data, group = "x", time = "z"),
    "Response variable `y` is invalid:\nx .+ family is not supported for <factor> variables\\."
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = bernoulli()),
             data = test_data, group = "x", time = "z"),
    "Response variable `y` is invalid:\nx .+ family is not supported for <factor> variables\\."
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = binomial()),
             data = test_data, group = "x", time = "z"),
    "Response variable `y` is invalid:\nx .+ family is not supported for <factor> variables\\."
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = poisson()),
             data = test_data, group = "x", time = "z"),
    "Response variable `y` is invalid:\nx .+ family is not supported for <factor> variables\\."
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = negbin()),
             data = test_data, group = "x", time = "z"),
    "Response variable `y` is invalid:\nx .+ family is not supported for <factor> variables\\."
  )
})

test_that("negative values for binomial, negbin and poisson fails", {
  test_data <- data.frame(y = c(-1, -2), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = binomial()),
             data = test_data, group = "x", time = "z"),
    "Response variable `y` is invalid:\nx .+ family supports only non-negative integers\\."
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = poisson()),
             data = test_data, group = "x", time = "z"),
    "Response variable `y` is invalid:\nx .+ family supports only non-negative integers\\."
  )
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = negbin()),
             data = test_data, group = "x", time = "z"),
    "Response variable `y` is invalid:\nx .+ family supports only non-negative integers\\."
  )
})

# Predict errors ----------------------------------------------------------

gaussian_example_small <- gaussian_example |> dplyr::filter(.data$time < 6)

test_that("newdata without group variable fails when there are groups", {
  gaussian_example_nogroup <- gaussian_example_small |> dplyr::select(!.data$id)
  expect_error(
    predict(gaussian_example_fit, newdata = gaussian_example_nogroup),
    "Can't find grouping variable `id` in `newdata`\\."
  )
})

test_that("newdata with new groups fails when there are groups", {
  gaussian_example_newgroup <- rbind(
    gaussian_example_small,
    data.frame(y = 1, x = 1, z = 0, id = 101, time = 1)
  )
  expect_error(
    predict(gaussian_example_fit, newdata = gaussian_example_newgroup),
    'Grouping variable `id` contains unknown levels:\nx Level "101" is not present in the original data\\.'
  )
})

test_that("newdata without time variable fails", {
  gaussian_example_notime <- gaussian_example_small |> dplyr::select(!.data$time)
  expect_error(
    predict(gaussian_example_fit, newdata = gaussian_example_notime),
    "Can't find time index variable `time` in `newdata`\\."
  )
})

test_that("newdata with new time points fails", {
  gaussian_example_newtime <- rbind(
    gaussian_example_small,
    data.frame(y = 1, x = 1, z = 0, id = 1, time = 31)
  )
  expect_error(
    predict(gaussian_example_fit, newdata = gaussian_example_newtime),
    'Time index variable `time` contains unknown time points:\nx Time point "31" is not present in the original data\\.'
  )
})

test_that("newdata with missing response fails", {
  gaussian_example_misresp <- gaussian_example_small |> dplyr::select(!.data$y)
  expect_error(
    predict(gaussian_example_fit, newdata = gaussian_example_misresp),
    "Can't find response variable `y` in `newdata`."
  )
})

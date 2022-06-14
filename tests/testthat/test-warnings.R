
# Data warnings -----------------------------------------------------------

test_that("ordered factor conversion to factor warns", {
  test_data <- data.frame(y = factor(c(1, 2), ordered = TRUE),
                          x = c(1, 1), z = c(1, 2))
  expect_warning(
    dynamite(dformula = obs(y ~ x, family = categorical()),
             data = test_data, group = "x", time = "z",
             debug = list(no_compile = TRUE)),
    "Response variable `y` is of class <ordered factor> whose channel is categorical:\ni `y` will be converted to an unordered factor\\."
  )
})

# Specials warnings -------------------------------------------------------

test_that("multiple intercept warns", {
  expect_warning(
    obs(y ~ 1 + varying(~1), family = gaussian()),
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

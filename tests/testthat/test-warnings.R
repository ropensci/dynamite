# Specials warnings -------------------------------------------------------

test_that("multiple intercept warns", {
  expect_warning(
    obs(y ~ 1 + varying(~1), family = gaussian()),
    "Both time-independent and time-varying intercept specified\\. Defaulting to time-varying intercept\\."
  )
})

test_that("deterministic fixed warns", {
  expect_warning(
    aux(numeric(y) ~ fixed(~ x)),
    "fixed\\(\\) definitions of a determinstic channel 'numeric\\(y\\)' will be ignored"
  )
})

test_that("deterministic varying warns", {
  expect_warning(
    aux(numeric(y) ~ varying(~ x)),
    "varying\\(\\) definitions of a determinstic channel 'numeric\\(y\\)' will be ignored"
  )
})

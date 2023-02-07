#' @srrstats {G5.10} Extended tests can be switched on via setting the
#'   environment variable DYNAMITE_EXTENDED_TESTS to "true".

run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

test_that("cmdstanr backend works for categorical model", {
  skip_if_not(run_extended_tests)
  set.seed(1)

  fit_dynamite <- update(categorical_example_fit,
    backend = "cmdstanr"
  )
  expect_equal(coef(fit_dynamite)$mean[1], -0.5,
    tolerance = 0.1,
    ignore_attr = TRUE
  )
  expect_error(get_code(fit_dynamite), NA)
})

test_that("stanc_options argument works", {
  skip_if_not(run_extended_tests)
  set.seed(1)

  fit <- dynamite(
    dformula = obs(y ~ -1 + varying(~x), family = "gaussian") +
      lags(type = "varying") +
      splines(df = 20),
    gaussian_example,
    "time",
    "id",
    chains = 1,
    refresh = 0,
    backend = "cmdstanr",
    stanc_options = list("O0")
  )
  expect_equal(summary(fit, parameter = "alpha_y")$mean[2], 1.5,
    tolerance = 0.1,
    ignore_attr = TRUE
  )
})

test_that("LOO and LFO works for AR(1) model estimated with cmdstanr", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  fit <- dynamite(obs(LakeHuron ~ 1, "gaussian") + lags(),
    data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
    time = "time",
    group = "id",
    chains = 1,
    iter = 2000,
    refresh = 0,
    backend = "cmdstanr"
  )
  l <- loo(fit)
  expect_equal(l$estimates,
    structure(
      c(
        -107.877842970846, 2.86041434691809, 215.755685941693,
        7.36848739076899, 0.561813071004331, 14.736974781538
      ),
      dim = 3:2,
      dimnames = list(c("elpd_loo", "p_loo", "looic"), c("Estimate", "SE"))
    ),
    tolerance = 1
  )
  expect_error(plot(l), NA)

  l <- lfo(fit, L = 20)
  expect_equal(l$ELPD, -90.4188604974201, tolerance = 1)
  expect_equal(l$ELPD_SE, 7.58842574523583, tolerance = 1)
  expect_error(plot(l), NA)
})

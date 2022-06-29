# Capture both message and output types
capture_all_output <- function(x) {
  utils::capture.output(
    utils::capture.output(x, type = "message"),
    type = "output"
  )
}

test_that("conversion to data.frame works", {
  expect_error(
    as.data.frame(gaussian_example_fit),
    NA
  )
  expect_error(
    as.data.frame(
      gaussian_example_fit,
      types = c(
        "alpha", "beta", "delta", "tau", "tau_alpha",
        "sigma_nu", "sigma", "phi", "nu", "omega", "omega_alpha"
       )
    ),
    NA
  )
})

test_that("coefficients can be extracted", {
  expect_error(
    coef(gaussian_example_fit, type = "beta"),
    NA
  )
  expect_error(
    coef(gaussian_example_fit, type = "delta", include_alpha = FALSE),
    NA
  )
  expect_error(
    coef(gaussian_example_fit, type = "delta", include_alpha = TRUE),
    NA
  )
})

test_that("fit object can be printed", {
  expect_error(
    capture_all_output(print(gaussian_example_fit)),
    NA
  )
  gaussian_example_fit_null <- gaussian_example_fit
  gaussian_example_fit_null$stanfit <- NULL
  expect_message(
    print(gaussian_example_fit_null),
    "No Stan model fit is available\\."
  )
})

test_that("default plot informs", {
  expect_message(
    plot(gaussian_example_fit),
    paste0(
      "Please use `plot_deltas\\(\\)`, `plot_betas\\(\\)`, or ",
      "`plot_nus\\(\\)` to produce plots of a <dynamitefit> object\\."
    )
  )
})

test_that("betas can be plotted", {
  expect_error(
    plot_betas(gaussian_example_fit),
    NA
  )
})

test_that("deltas can be plotted", {
  expect_error(
    plot_deltas(gaussian_example_fit),
    NA
  )
})

test_that("nus can be plotted", {
  expect_error(
    plot_nus(gaussian_example_fit),
    NA
  )
})

test_that("formula can be extracted", {
  expect_error(
    formula(gaussian_example_fit),
    NA
  )
  f <- obs(y ~ 1 + varying(~-1 + x), family = "gaussian") +
    obs(w ~ x + lag(y), family = "exponential") +
    aux(numeric(f) ~ w - 1) +
    aux(numeric(d) ~ y + 1) +
    lags(k = 1) +
    splines(df = 5)
  set.seed(0)
  d <- data.frame(
    y = rnorm(5),
    w = rexp(5),
    x = rnorm(5),
    z = seq.int(5)
  )
  fit <- dynamite(f, d, time = "z", debug = list(no_compile = TRUE))
  expect_error(
    formula(fit),
    NA
  )
})

test_that("MCMC diagnostics can be computed", {
  expect_error(
    capture_all_output(mcmc_diagnostics(gaussian_example_fit)),
    NA
  )
  gaussian_example_fit_null <- gaussian_example_fit
  gaussian_example_fit_null$stanfit <- NULL
  expect_message(
    mcmc_diagnostics(gaussian_example_fit_null),
    "No Stan model fit is available\\."
  )
})

test_that("gets can be got", {
  expect_error(
    get_code(obs(y ~ -1 + z + varying(~ x + lag(y)), family = "gaussian",
             random_intercept = TRUE) + splines(df = 20),
             gaussian_example, "id", "time"),
    NA
  )
  expect_error(
    get_code(gaussian_example_fit),
    NA
  )
  expect_error(
    get_data(obs(y ~ -1 + z + varying(~ x + lag(y)), family = "gaussian",
             random_intercept = TRUE) + splines(df = 20),
             gaussian_example, "id", "time"),
    NA
  )
  expect_error(
    get_data(gaussian_example_fit),
    NA
  )
})

test_that("credible intervals can be computed", {
  expect_error(
    confint(gaussian_example_fit),
    NA
  )
})

test_that("number of observations can be extracted", {
  expect_error(
    nobs(gaussian_example_fit),
    NA
  )
  expect_error(
    nobs(gaussian_example_single_fit),
    NA
  )
})

test_that("summary can be extracted", {
  expect_error(
    summary(gaussian_example_fit),
    NA
  )
})

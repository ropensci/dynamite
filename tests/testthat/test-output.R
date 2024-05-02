data.table::setDTthreads(1) # For CRAN

# Capture both message and output types
capture_all_output <- function(x) {
  utils::capture.output(
    utils::capture.output(x, type = "message"),
    type = "output"
  )
}

test_that("dynamiteformula can be printed", {
  f <- obs(y ~ x + random(~1), family = "gaussian") +
    lags(k = c(1, 3)) +
    random_spec()
  expect_output(
    print(f),
    paste0(
      "  Family   Formula           \n",
      "y gaussian y ~ x \\+ random\\(~1\\)\n",
      "\n",
      "Lagged responses added as fixed predictors with: k = 1, 3\n",
      "Correlated random effects added for response\\(s\\): y"
    )
  )
})

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
  expect_error(
    capture_all_output(print(gaussian_example_fit, full_diagnostics = TRUE)),
    NA
  )
  gaussian_example_fit_null <- gaussian_example_fit
  gaussian_example_fit_null$stanfit <- NULL
  expect_output(
    print(gaussian_example_fit_null),
    "No Stan model fit is available"
  )
})

test_that("default plot works", {
  expect_error(
    plot(gaussian_example_fit, type = "beta"),
    NA
  )
})

test_that("default formula plot works", {
  f <- obs(y ~ x + lag(logz) + lag(y, 2) + lag(w), family = "gaussian") +
    obs(w ~ lag(x) + z, family = "gaussian") +
    aux(numeric(logz) ~ log(z))
  expect_error(
    plot(f),
    NA
  )
  expect_error(
    plot(f, show_auxiliary = TRUE),
    NA
  )
  expect_error(
    plot(f, show_covariates = TRUE),
    NA
  )
  expect_error(
    plot(f, show_auxiliary = TRUE, show_covariates = TRUE),
    NA
  )
  multichannel_formula <- obs(g ~ lag(g) + lag(logp), family = "gaussian") +
    obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
    obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
    aux(numeric(logp) ~ log(p + 1))
  expect_error(
    plot(multichannel_formula),
    NA
  )
  expect_error(
    plot(multichannel_formula, show_auxiliary = TRUE),
    NA
  )
  expect_error(
    plot(multichannel_formula, show_covariates = TRUE),
    NA
  )
  expect_error(
    plot(multichannel_formula, show_auxiliary = TRUE, show_covariates = TRUE),
    NA
  )
})

test_that("tikz formula plot works", {
  f <- obs(y ~ x + lag(logz) + lag(y, 2) + lag(w), family = "gaussian") +
    obs(w ~ lag(x) + z, family = "gaussian") +
    aux(numeric(logz) ~ log(z))
  expect_error(
    plot(f, tikz = TRUE),
    NA
  )
  expect_error(
    plot(f, show_auxiliary = TRUE, tikz = TRUE),
    NA
  )
  expect_error(
    plot(f, show_covariates = TRUE, tikz = TRUE),
    NA
  )
  expect_error(
    plot(f, show_auxiliary = TRUE, show_covariates = TRUE, tikz = TRUE),
    NA
  )
  multichannel_formula <- obs(g ~ lag(g) + lag(logp), family = "gaussian") +
    obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
    obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
    aux(numeric(logp) ~ log(p + 1))
  expect_error(
    plot(multichannel_formula, tikz = TRUE),
    NA
  )
  expect_error(
    plot(multichannel_formula, show_auxiliary = TRUE, tikz = TRUE),
    NA
  )
  expect_error(
    plot(multichannel_formula, show_covariates = TRUE, tikz = TRUE),
    NA
  )
  expect_error(
    plot(
      multichannel_formula,
      show_auxiliary = TRUE,
      show_covariates = TRUE,
      tikz = TRUE
    ),
    NA
  )
})

test_that("betas can be plotted", {
  expect_error(
    plot(gaussian_example_fit, types = "beta"),
    NA
  )
  expect_error(
    plot(categorical_example_fit, types = "beta"),
    NA
  )
})

test_that("deltas can be plotted", {
  expect_error(
    plot(gaussian_example_fit, types = "delta"),
    NA
  )
})

test_that("nus can be plotted", {
  expect_error(
    plot(gaussian_example_fit, types = "nu"),
    NA
  )
})

test_that("formula can be extracted", {
  expect_error(
    formula(gaussian_example_fit),
    NA
  )
  f <- obs(y ~ 1 + varying(~ -1 + x), family = "gaussian") +
    obs(w ~ x + lag(y), family = "exponential") +
    aux(numeric(f) ~ w - 1) +
    aux(numeric(d) ~ y + 1) +
    lags(k = 1) +
    splines(df = 5)

  set.seed(0)
  d <- data.frame(
    y = rnorm(6),
    w = rexp(6),
    x = rnorm(6),
    z = seq.int(6)
  )
  fit <- dynamite(f, d, time = "z", debug = list(no_compile = TRUE))
  expect_error(
    formula(fit),
    NA
  )
  expect_error(
    formula(multichannel_example_fit),
    NA
  )
})

test_that("formula extraction is correct", {
  expect_identical(
    deparse1(formula(gaussian_example_fit)),
    paste0(
      "obs(y ~ -1 + z + varying(~x + lag(y)) + ",
      "random(~1), family = \"gaussian\") + ",
      "splines(df = 20, degree = 3, lb_tau = 0, ",
      "noncentered = FALSE, override = FALSE) + ",
      "random_spec(correlated = FALSE, noncentered = TRUE)"
    )
  )
  expect_identical(
    deparse1(formula(multichannel_example_fit)),
    paste0(
      "obs(g ~ lag(g) + lag(logp), family = \"gaussian\") + ",
      "obs(p ~ lag(g) + lag(logp) + lag(b), family = \"poisson\") + ",
      "obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = \"bernoulli\") + ",
      "aux(numeric(logp) ~ log(p + 1))"
    )
  )
  expect_identical(
    deparse1(formula(categorical_example_fit)),
    paste0(
      "obs(x ~ z + lag(x) + lag(y), family = \"categorical\") + ",
      "obs(y ~ z + lag(x) + lag(y), family = \"categorical\")"
    )
  )

  f <- obs(z ~ w1, family = "gaussian") +
    obs(c(w1, w2, w3) ~ 1 | y | varying(~-1 + x), family = "mvgaussian") +
    lags(k = 1, type = "varying") +
    splines()
  d <- data.frame(
    x = rnorm(10),
    y = rnorm(10),
    z = rnorm(10),
    w1 = rnorm(10),
    w2 = rnorm(10),
    w3 = rnorm(10),
    time = 1:10,
    id = rep(1, 10)
  )
  fit <- dynamite(f, d, "time", "id", debug = list(no_compile = TRUE))
  expect_identical(
    deparse1(formula(fit)),
    paste0(
      "obs(z ~ w1, family = \"gaussian\") + ",
      "obs(c(w1, w2, w3) ~ 1 | y | varying(~-1 + x), family = \"mvgaussian\") + ",
      "lags(k = 1, type = \"varying\") + ",
      "splines(df = NULL, degree = 3, lb_tau = c(0, 0, 0, 0), ",
      "noncentered = c(FALSE, FALSE, FALSE, FALSE), override = FALSE)"
    )
  )

  f <- obs(c(a, b, c) ~ d + trials(n), family = "multinomial") +
    lags(k = 1, type = "fixed")
  d <- data.frame(
    a = sample.int(5, 10, replace = TRUE),
    b = sample.int(10, 10, replace = TRUE),
    c = sample.int(12, 10, replace = TRUE),
    d = rnorm(10),
    time = 1:10,
    id = rep(1, 10)
  )
  d$n <- d$a + d$b + d$c
  fit <- dynamite(f, d, "time", "id", debug = list(no_compile = TRUE))
  expect_identical(
    deparse1(formula(fit)),
    paste0(
      "obs(c(a, b, c) ~ d + trials(n), family = \"multinomial\") + ",
      "lags(k = 1, type = \"fixed\")"
    )
  )
})

test_that("MCMC diagnostics can be computed", {
  expect_error(
    capture_all_output(mcmc_diagnostics(gaussian_example_fit)),
    NA
  )
  gaussian_example_fit_null <- gaussian_example_fit
  gaussian_example_fit_null$stanfit <- NULL
  expect_output(
    mcmc_diagnostics(gaussian_example_fit_null),
    "No Stan model fit is available\\."
  )
})

test_that("gets can be got", {
  expect_error(
    get_code(
      obs(y ~ -1 + z + varying(~ x + lag(y)) +
        random(~1), family = "gaussian") + random_spec() + splines(df = 20),
      gaussian_example, time =  "time", group = "id"
    ),
    NA
  )
  expect_error(
    code1 <- get_code(gaussian_example_fit),
    NA
  )
  gaussian_example_fit_null <- gaussian_example_fit
  gaussian_example_fit_null$stanfit <- NULL
  expect_error(
    code2 <- get_code(gaussian_example_fit_null),
    NA
  )
  expect_error(
    get_code(gaussian_example_fit, blocks = c("parameters", "model")),
    NA
  )
  expect_error(
    get_data(
      obs(y ~ -1 + z + varying(~ x + lag(y)) + random(~1),
        family = "gaussian"
      ) + random_spec() + splines(df = 20),
      gaussian_example, time = "time", group = "id"
    ),
    NA
  )
  expect_error(
    get_data(gaussian_example_fit),
    NA
  )
  expect_equal(
    get_parameter_names(categorical_example_fit),
    c("alpha_x", "beta_x_z", "beta_x_x_lag1B", "beta_x_x_lag1C",
      "beta_x_y_lag1b", "beta_x_y_lag1c", "alpha_y", "beta_y_z",
      "beta_y_x_lag1B", "beta_y_x_lag1C", "beta_y_y_lag1b", "beta_y_y_lag1c")
  )
  expect_equal(
    get_parameter_dims(categorical_example_fit),
    list(
      beta_x_B = 5L,
      a_x_B = 1L,
      beta_x_C = 5L,
      a_x_C = 1L,
      beta_y_b = 5L,
      a_y_b = 1L,
      beta_y_c = 5L,
      a_y_c = 1L
    )
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
    n <- nobs(gaussian_example_fit),
    NA
  )
  expect_equal(n, 1450L)
  expect_error(
    n <- nobs(gaussian_example_single_fit),
    NA
  )
  expect_equal(n, 29L)
  expect_error(
    n <- nobs(multichannel_example_fit),
    NA
  )
  expect_equal(n, 2850L)
  expect_error(
    n <- nobs(categorical_example_fit),
    NA
  )
  expect_equal(n, 3800L)
})

test_that("summary can be extracted", {
  expect_error(
    summary(gaussian_example_fit),
    NA
  )
})

test_that("number of draws can be extraced", {
  expect_error(
    ndraws(gaussian_example_fit),
    NA
  )
})

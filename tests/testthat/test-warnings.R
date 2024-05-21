data.table::setDTthreads(1) # For CRAN

# Data warnings -----------------------------------------------------------

test_that("factor time conversion warns", {
  test_data <- data.frame(
    y = c(1, 2, 3),
    x = c(1, 1, 2),
    z = factor(c(1, 2, 3))
  )
  expect_warning(
    dynamite(
      dformula = obs(y ~ x, family = "negbin"),
      data = test_data, group = "x", time = "z",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Time index variable `z` is a <factor>:\n",
      "i Converting the variable to <integer> based on its levels\\."
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
  test_data <- data.frame(
    y = rnorm(3),
    x = rnorm(3),
    z = rnorm(3),
    w = rnorm(3)
  )
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
    z = factor(1:3)
  )
  expect_warning(
    full_model.matrix(f, test_data, TRUE),
    paste0(
      "Predictor `z1` contains only zeros in the complete case rows of the ",
      "design matrix for the channel `y`\\."
    )
  )
})

test_that("deterministic channel insufficient initial values warns", {
  expect_warning(
    dynamite(
      dformula = obs(y ~ x, family = "gaussian") + aux(numeric(d) ~ lag(d, 1)),
      data = data.frame(y = c(1, 2, 1), x = c(1, 2, 3), z = c(1, 2, 3)),
      group = NULL,
      time = "z",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Deterministic channel `d` has a maximum lag of 1 but ",
      "you've supplied no initial values:\n",
      "i This may result in NA values for `d`\\."
    )
  )
})

# Specials warnings -------------------------------------------------------

test_that("multiple intercept warns", {
  expect_warning(
    obs(y ~ 1 + varying(~1), family = "gaussian"),
    paste0(
      "Both time-constant and time-varying intercept specified:\n",
      "i Defaulting to time-varying intercept\\."
    )
  )
})

test_data <- data.frame(
  y = rnorm(10),
  x = rnorm(10),
  time = 1:5,
  id = rep(1:2, each = 5)
)

debug <- list(no_compile = TRUE)

test_that("time-varying intercept is removed", {
  expect_warning(
    dynamite(
      obs(y ~ -1 + x + varying(~1), family = "gaussian") +
        lfactor() +
        splines(4),
      test_data,
      "time",
      "id",
      debug = debug
    ),
    paste0(
      "The common time-varying intercept term of channel `y` was removed ",
      "as channel predictors contain latent factor specified with ",
      "`nonzero_kappa` as TRUE\\."
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
      "You've supplied `n_draws` = 500 but there are only ",
      ndraws(gaussian_example_fit),
      " samples available:\n",
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
    paste0(
      "Time index variable `time` of `newdata` has gaps:\n",
      "i Filling the `newdata` to regular time points\\. This will lead to ",
      "propagation of NA values if the model contains exogenous predictors ",
      "and `impute` is \"none\"\\."
    )
  )
  newdata <- gaussian_example |>
    dplyr::filter(id == 1) |>
    dplyr::mutate(y = ifelse(time > 5, NA, y)) |>
    dplyr::filter(time < 3 | time > 10)
  # capture due to multiple warnings
  w <- capture_warnings(
    predict(gaussian_example_single_fit, newdata = newdata, ndraws = 1)
  )
  expect_match(
    w[1L],
    paste0(
      "Time index variable `time` of `newdata` has gaps:\n",
      "i Filling the `newdata` to regular time points\\. This will lead to ",
      "propagation of NA values if the model contains exogenous predictors ",
      "and `impute` is \"none\"\\.|NAs produced"
    )
  )
})

# Stan warnings -----------------------------------------------------------

test_that("unrecognized arguments warns", {
  expect_warning(
    dynamite(
      obs(y ~ x, family = "gaussian") +
        splines(4),
      test_data,
      "time",
      "id",
      debug = debug,
      strange_arg1 = 1L,
      strange_arg2 = 1L,
    ),
    paste0(
      "Arguments `strange_arg1` and `strange_arg2` passed to rstan sampling ",
      "function are not recognized and will be ignored\\."
    )
  )
})

test_that("categorical non-glm availability warns", {
  expect_warning(
    mockthat::with_mock(
      dynamite_model = function(...) NULL,
      dynamite_sampling = function(...) NULL,
      stan_supports_categorical_logit_glm = function(...) FALSE,
      dynamite(
        dformula = obs(y ~ 1, family = "categorical"),
        data = data.frame(y = c("A", "B"), time = c(1, 2)),
        time = "time"
      )
    ),
    paste0(
      "Efficient GLM variant of the categorical likelihood is not available ",
      "in this version of rstan\\.\n",
      "i For more efficient sampling, please install ",
      "a newer version of rstan\\."
    )
  )
})

test_that("windows and old rstan warns on attach", {
  out <- mockthat::with_mock(
    stan_version = function(...) "2.23",
    is_windows = function(...) TRUE,
    R_version = function(...) "4.2.0",
     capture.output(startup(), type = "message")
  )
  expect_match(
    out[1L],
    paste0(
      "Please update your `rstan` and `StanHeaders` installations before ",
      "using `dynamite` with the `rstan` backend by running:"
    )
  )
})

# Plot warnings -----------------------------------------------------------

test_that("too many parameters warns in plot", {
  expect_warning(
    plot(gaussian_example_fit, types = "nu"),
    paste0(
      "Number of parameters to be plotted \\(50\\) exceeds the maximum ",
      "number of parameters \\(20\\) for parameters of type `nu`\\. ",
      "The remaining parameters of this type will not be plotted\\.\n",
      "i Please increase `n_params` to plot more parameters\\."
    )
  )
  expect_warning(
    plot(gaussian_example_fit, types = "nu", plot_type = "trace"),
    paste0(
      "Number of parameters to be plotted \\(50\\) exceeds the maximum ",
      "number of parameters \\(5\\)\\. ",
      "The remaining parameters will not be plotted\\.\n",
      "i Please increase `n_params` to plot more parameters\\."
    )
  )
})

# Deprecated --------------------------------------------------------------

test_that("deprecated functions warn", {
  expect_warning(
    plot_betas(gaussian_example_fit),
    "'plot_betas' is deprecated"
  )
  expect_warning(
    plot_deltas(gaussian_example_fit),
    "'plot_deltas' is deprecated"
  )
  expect_warning(
    plot_nus(gaussian_example_fit, n_params = 10),
    "'plot_nus' is deprecated"
  )
  expect_warning(
    try(plot_lambdas(gaussian_example_fit), silent = TRUE),
    "'plot_lambdas' is deprecated"
  )
  expect_warning(
    try(plot_psis(gaussian_example_fit), silent = TRUE),
    "'plot_psis' is deprecated"
  )
})

test_that("deprecated cmdstanr arguments warn", {
  dots <- list(seed = 0, cores = 4, num_sampling = 1000)
  expect_warning(
    check_stan_args(dots, verbose = TRUE, backend = "cmdstanr"),
  )
})

obs_test <- obs(y ~ x + w, family = "gaussian")

set.seed(0)
timepoints <- 10
individuals <- 5
total_obs <- timepoints * individuals

test_data <- data.frame(
  time = 1:timepoints,
  group = gl(individuals, timepoints),
  offset = sample(50:100, size = total_obs, replace = TRUE),
  trials = sample(50:100, size = total_obs, replace = TRUE)
) |>
  dplyr::mutate(
    y1 = as.factor(sample(5, size = total_obs, replace = TRUE)),
    y2 = rnorm(n = total_obs, mean = 1, sd = 2),
    y3 = rbinom(n = total_obs, size = trials, prob = 0.75),
    y4 = rbinom(n = total_obs, size = 1, prob = 0.66),
    y5 = rnbinom(n = total_obs, size = 100, prob = 0.33),
    y6 = rpois(n = total_obs, lambda = log(offset) + 1),
    y7 = rexp(n = total_obs, rate = 0.1),
    y8 = rgamma(n = total_obs, shape = 2, rate = 2 * exp(-5)),
    y9 = rbeta(n = total_obs, 6, 4),
    x1 = sample(letters[1:4], size = total_obs, replace = TRUE),
    x2 = rnorm(total_obs),
    x3 = as.factor(sample(4, size = total_obs, replace = TRUE))
  )

# Formula errors ----------------------------------------------------------

test_that("nonformula to dynamiteformula fails", {
  expect_error(
    obs(formula = numeric()),
    "Argument `formula` must be a <formula> object\\."
  )
})

test_that("noncharacter family fails", {
  expect_error(
    obs(y ~ x, family = data.frame()),
    "Argument `family` must be a single <character> string\\."
  )
})

test_that("unsupported family fails", {
  expect_error(
    obs(y ~ x, family = "unknown_distr"),
    'Family "unknown_distr" is not supported\\.'
  )
})

test_that("as-is use fails", {
  expect_error(
    obs(y ~ I(x), family = "gaussian"),
    "`I\\(\\.\\)` is not supported by `dynamiteformula\\(\\)`\\."
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

test_that("duplicate lags definition fails", {
  expect_error(
    obs_test + lags() + lags(),
    "Multiple definitions for lags\\."
  )
})

test_that("adding dynamiteformulas with existing lag definitions fails", {
  obs_lhs <- obs_test + lags(k = 1)
  obs_rhs <- obs(z ~ x, family = "gaussian") + lags(k = 2)
  expect_error(
    obs_lhs + obs_rhs,
    "Both dynamiteformulas contain a lags definition\\."
  )
})

test_that("adding dynamiteformulas with existing splines definitions fails", {
  obs_lhs <- obs_test + splines()
  obs_rhs <- obs(z ~ x, family = "gaussian") + splines()
  expect_error(
    obs_lhs + obs_rhs,
    "Both dynamiteformulas contain a splines definition\\."
  )
})

test_that("simultaneity fails", {
  obs_lhs <-
    obs(q ~ w + e + r + lag(i), family = "gaussian") +
    obs(t ~ y + u, family = "gaussian") +
    obs(i ~ o + p + a + lag(f), family = "gaussian")
  obs_rhs <-
    obs(f ~ h + l + lag(x), family = "gaussian") +
    obs(x ~ q + z, family = "gaussian")
  expect_error(
    obs_rhs + obs_lhs,
    paste0(
      "Simultaneous regression is not supported:\n",
      "x Response variable `q` appears in the formula of `x`\\."
    )
  )
  # should fail for deterministic as well
  expect_error(
    obs(y ~ x, family = "gaussian") + aux(integer(x) ~ y),
    paste0(
      "Simultaneous regression is not supported:\n",
      "x Response variable `x` appears in the formula of `y`\\."
    )
  )
})

test_that("adding nondynamiteformula to dynamiteformula fails", {
  expect_error(
    obs_test + 1.0,
    paste(
      "Unable to add an object of class <numeric>",
      "to an object of class <dynamiteformula>\\."
    )
  )
})

test_that("plus method fails for nondynamiteformula", {
  expect_error(
    `+.dynamiteformula`(data.frame(), numeric()),
    paste(
      "Method `\\+\\.dynamiteformula\\(\\)` is not supported",
      "for <data.frame> objects\\."
    )
  )
})

test_that("categorical random intercept fails", {
  expect_error(
    dynamite(obs(y ~ x, family = "categorical") + random(),
      data = data.frame(y=factor(1:4), x= runif(4), id = 1, time = 1:4),
      "id", "time",
      debug = list(no_compile = TRUE)),
    paste0(
      "No valid responses for random intercept component:\n",
      "x Random intercepts are not supported for the categorical family\\."
    )
  )
})

test_that("negative lb_tau fails", {
  expect_error(
    obs_test + splines(lb_tau = -1.0),
    "Argument `lb_tau` must be a <numeric> vector of non-negative values\\."
  )
})

test_that("time-varying definitions without splines fails", {
  obs_varying <- obs(y ~ 1 + varying(~-1 + x), family = "gaussian")
  test_data <- data.frame(y = c(1, 2, 3),
                          x = c(0.5, -1, 0.25),
                          z = c(1, 2, 3))
  expect_error(
    dynamite(obs_varying, test_data, time = "z"),
    paste(
      "Model for response variable `y` contains time-varying definitions",
      "but splines have not been defined\\."
    )
  )
})


test_that("noncentered definition throws error if not of correct length", {
  expect_error(
    obs_all_alpha <- obs(y1 ~ -1 + varying(~ x1), family = "categorical") +
      obs(x3 ~ varying(~ -1 + x1), family = "categorical") +
      obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ lag(x3) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~-1 + x2), family = "bernoulli") +
      obs(y9 ~ -1 + x1 + varying(~x2), family = "beta") +
      splines(df = 5, noncentered = rep(TRUE, 3)),
    NA
  )
  expect_error(
    dynamite(obs_all_alpha, test_data,
      "group", "time", debug = debug),
    paste("Length of the `noncentered` argument of `splines\\(\\)` function",
      "is not equal to 1 or 6, the number of the channels\\."
    )
  )
})

test_that("lb_tau definition throws error if not of correct length", {
  expect_error(
    obs_all_alpha <- obs(y1 ~ -1 + varying(~ x1), family = "categorical") +
      obs(x3 ~ varying(~ -1 + x1), family = "categorical") +
      obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ lag(x3) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~-1 + x2), family = "bernoulli") +
      obs(y9 ~ -1 + x1 + varying(~x2), family = "beta") +
      splines(df = 5, lb_tau = rep(1, 3)),
    NA
  )
  expect_error(
    dynamite(obs_all_alpha, test_data,
      "group", "time", debug = debug),
    paste("Length of the `lb_tau` argument of `splines\\(\\)` function is not",
      "equal to 1 or 6, the number of the channels\\."
    )
  )
})
# Formula specials errors -------------------------------------------------

test_that("no intercept or predictors fails", {
  expect_error(
    obs(y ~ -1, family = "gaussian"),
    paste0(
      "Invalid formula for response variable `y`:\n",
      "x There are no predictors nor an intercept term\\."
    )
  )
})

test_that("past in the middle of formula fails", {
  expect_error(
    aux(y ~ x + past(0) + z),
    "Past values term must be the last term of the formula\\."
  )
})

test_that("binomial channel without a trials term fails", {
  expect_error(
    obs(y ~ x, family = "binomial"),
    "Formula for a binomial channel must include a trials term\\."
  )
})

# Data errors -------------------------------------------------------------

test_that("data is not data.frame fails", {
  expect_error(
    dynamite(dformula = obs_test, data = list()),
    "Argument `data` must be a <data.frame> object\\.")
})

test_that("group variable not in data fails", {
  expect_error(
    dynamite(dformula = obs_test,
             data = data.frame(y = 1, x = 1), group = "z"),
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
    dynamite(dformula = obs_test,
             data = data.frame(y = 1, x = 1), time = "z"),
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

test_that("duplicated time points fail", {
  # groups
  expect_error(
    dynamite(
      dformula = obs(y ~ x, family = "gaussian"),
      data = data.frame(
       y = rep(1, 9),
       x = gl(3, 3),
       z = c(1, 2, 2, 1, 2, 3, 1, 3, 3)
      ),
      group = "x",
      time = "z"
    ),
    paste0(
      "Each time index must correspond to a single observation per group:\n",
      "x Groups `1` and `3` of `x` have duplicate observations\\."
    )
  )
  # no groups
  expect_error(
    dynamite(
      dformula = obs(y ~ x, family = "gaussian"),
      data = data.frame(
        y = rep(1, 3),
        z = c(1, 2, 2)
      ),
      time = "z"
    ),
    "Each time index must correspond to a single observation\\."
  )
})

test_that("missing lag variable fails", {
  expect_error(
    dynamite(dformula = obs(y ~ lag(d, 1), family = "gaussian"),
             data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
             group = "x", time = "z"),
    paste0(
      "Unable to construct lagged values of `d`:\n",
      "x Can't find such variables in `data`\\."
    )
  )
})

test_that("missing predictor fails", {
  expect_error(
    dynamite(dformula = obs(y ~ w, family = "gaussian"),
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
    paste0(
      "Unable to evaluate definitions of deterministic channels:\n",
      "i Some variables are possibly missing or incorrect\\."
    )
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
    paste0(
      "Deterministic channel `d` requires 1 initial value:\n",
       "x You've supplied no initial values\\."
    )
  )
})

# Data type errors --------------------------------------------------------

#' @srrstats {G2.11, G2.12} Tests for unsupported column types.
test_that("invalid column types fail", {
  test_data <- data.frame(y = c(1i, 2i), x = c(1, 1), z = c(1, 2))
  test_data$w <- c(list(a = 1), list(b = 2))
  test_data$d <- as.raw(c(40, 20))
  expect_error(
    dynamite(dformula = obs(y ~ x, family = "gaussian"),
             data = test_data, group = "x", time = "z"),
    paste0(
      "Columns `y`, `w`, and `d` of `data` are invalid:\n",
      "x Column types <complex/list/raw> are not supported\\."
    )
  )
})

test_that("non-finite values in data fail", {
  test_data <- data.frame(y = c(1, Inf), x = c(1, 1),
                          z = c(1, 2), w = c(-Inf, 2), u = c(1, Inf))
  expect_error(
    dynamite(dformula = obs(y ~ x, family = "gaussian"),
             data = test_data, group = "x", time = "z"),
    "Non-finite values in variables `y`, `w`, and `u` of `data`\\."
  )
})

test_that("non-factor categorical response fails", {
  test_data <- data.frame(y = c(0, 1), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = "categorical"),
             data = test_data, group = "x", time = "z"),
    paste0(
      "Response variable `y` is invalid:\n",
      "x Categorical family supports only <factor> variables\\."
    )
  )
})

test_that("factor types for non-categorical families fails", {
  test_data <- data.frame(y = factor(c(0, 1)), x = c(1, 1), z = c(1, 2))
  families <- c("gaussian", "exponential", "gamma", "beta",
                "bernoulli", "binomial", "poisson", "negbin")
  for (f in families) {
    form <- ifelse_(identical(f, "binomial"), y ~ 1 + trials(x), y ~ 1)
    expect_error(
      dynamite(dformula = obs(form, family = f),
               data = test_data, group = "x", time = "z"),
      paste0(
        "Response variable `y` is invalid:\n",
        "x .+ family is not supported for <factor> variables\\."
      )
    )
  }
})

test_that("negative values for distributions with positive support fails", {
  test_data <- data.frame(y = c(-1, -2), x = c(1, 1), z = c(1, 2))
  families <- c("exponential", "gamma", "binomial", "negbin", "poisson")
  for (f in families) {
    form <- ifelse_(identical(f, "binomial"), y ~ 1 + trials(x), y ~ 1)
    expect_error(
      dynamite(dformula = obs(form, family = f),
               data = test_data, group = "x", time = "z"),
      paste0(
        "Response variable `y` is invalid:\n",
        "x .+ family supports only non-negative .+\\."
      )
    )
  }
})

test_that("bernoulli without 0/1 values fails", {
  test_data <- data.frame(y = c(2, 3), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = "bernoulli"),
             data = test_data, group = "x", time = "z"),
    paste0(
      "Response variable `y` is invalid:\n",
      "x Bernoulli family supports only 0/1 integers\\."
    )
  )
})

test_that("beta without (0, 1) values fails", {
  test_data <- data.frame(y = c(2, 3), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(dformula = obs(y ~ 1, family = "beta"),
      data = test_data, group = "x", time = "z"),
    paste0(
      "Response variable `y` is invalid:\n",
      "x Beta family supports only values on open interval \\(0, 1\\)\\."
    )
  )
})
# Lag errors --------------------------------------------------------------

test_that("invalid lagged value definition fails", {
  expect_error(
    complete_lags(quote(lag(y, a:b))),
    "Invalid shift value expression `a:b`\\."
  )
})

test_that("non coerceable shift value fails", {
  expect_error(
    complete_lags(quote(lag(y, 'a'))),
    'Unable to coerce shift value to <integer> in `lag\\(y, "a"\\)`\\.'
  )
})

test_that("multiple shift values fail", {
  expect_error(
    complete_lags(quote(lag(y, 1:2))),
    paste0(
      "Shift value must be a single <integer> in `lag\\(\\)`:\n",
      "x Multiple shift values were found in `lag\\(y, 1:2\\)`\\."
    )
  )
})

test_that("negative lag shift value fails", {
  expect_error(
    complete_lags(quote(lag(y, -1))),
    paste0(
      "Shift value must be positive in `lag\\(\\)`:\n",
      "x Nonpositive shift value was found in `lag\\(y, -1\\)`\\."
    )
  )
})

test_that("too many arguments to lag fails", {
  expect_error(
    complete_lags(quote(lag(y, 1, 2))),
    paste0(
      "Invalid lag definition `lag\\(y, 1, 2\\)`:\n",
      "x Too many arguments supplied to `lag\\(\\)`\\."
    )
  )
})

# Output errors -----------------------------------------------------------

test_that("output for non dynamitefit objects fails", {
  expect_error(
    as.data.frame.dynamitefit(x = 1L),
    "Argument `x` must be a <dynamitefit> object\\."
  )
  expect_error(
    as_draws_df.dynamitefit(x = 1L),
    "Argument `x` must be a <dynamitefit> object\\."
  )
  expect_error(
    plot.dynamitefit(x = 1L),
    "Argument `x` must be a <dynamitefit> object\\."
  )
  expect_error(
    plot_deltas(x = 1L),
    "Argument `x` must be a <dynamitefit> object\\."
  )
  expect_error(
    plot_betas(x = 1L),
    "Argument `x` must be a <dynamitefit> object\\."
  )
  expect_error(
    plot_nus(x = 1L),
    "Argument `x` must be a <dynamitefit> object\\."
  )
  expect_error(
    print.dynamitefit(x = 1L),
    "Argument `x` must be a <dynamitefit> object\\."
  )
  expect_error(
    mcmc_diagnostics(x = 1L),
    "Argument `x` must be a <dynamitefit> object\\."
  )
  expect_error(
    nobs.dynamitefit(object = 1L),
    "Argument `object` must be a <dynamitefit> object\\."
  )
})

test_that("Invalid responses fail", {
  expect_error(
    as.data.frame.dynamitefit(
      gaussian_example_fit,
      responses = "resp"
    ),
    "Model does not contain response variable `resp`\\."
  )
})

test_that("Invalid confint level fails", {
  expect_error(
    confint.dynamitefit(gaussian_example_fit, level = -0.1),
    "Argument `level` must be a single <numeric> value between 0 and 1\\."
  )
})

# Predict errors ----------------------------------------------------------

gaussian_example_small <- gaussian_example |> dplyr::filter(.data$time < 6)

test_that("newdata without group variable fails when there are groups", {
  gaussian_example_nogroup <- gaussian_example_small |>
    dplyr::select(!.data$id)
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
    paste0(
      "Grouping variable `id` contains unknown levels:\n",
      "x Level \"101\" is not present in the original data\\."
    )
  )
})

test_that("newdata without time variable fails", {
  gaussian_example_notime <- gaussian_example_small |>
    dplyr::select(!.data$time)
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
    paste0(
      "Time index variable `time` contains unknown time points:\n",
      "x Time point \"31\" is not present in the original data\\."
    )
  )
})

test_that("newdata with duplicated time points fails", {
  # groups
  gaussian_example_duplicated <- rbind(
    gaussian_example_small,
    data.frame(y = 1, x = 1, z = 0, id = 1, time = 1)
  )
  expect_error(
    predict(gaussian_example_fit, newdata = gaussian_example_duplicated),
    paste0(
      "Each time index must correspond to a single observation per group:\n",
      "x Group `1` of `id` has duplicate observations\\."
    )
  )
  # no groups
  gaussian_example_duplicated <- rbind(
    gaussian_example_small |>
      dplyr::filter(.data$id == 1) |>
      dplyr::select(!.data$id),
    data.frame(y = 1, x = 1, z = 0, time = 1)
  )
  expect_error(
    predict(gaussian_example_single_fit, newdata = gaussian_example_duplicated),
    "Each time index must correspond to a single observation\\."
  )
})

test_that("newdata with unknown factor levels fails", {
  categorical_example_newlevel <- categorical_example |>
    dplyr::mutate(x = dplyr::recode(x, "C" = "D"))
  expect_error(
    predict(categorical_example_fit, newdata = categorical_example_newlevel),
    paste0(
      "<factor> variable `x` in `newdata` has new levels:\n",
      "x Level \"D\" is not present in the original data\\."
    )
  )
})

test_that("newdata with missing response fails", {
  gaussian_example_misresp <- gaussian_example_small |> dplyr::select(!.data$y)
  expect_error(
    predict(gaussian_example_fit, newdata = gaussian_example_misresp),
    "Can't find response variable `y` in `newdata`."
  )
})

test_that("newdata that is not a data.frame fails", {
  expect_error(
    predict(gaussian_example_fit, newdata = 1L),
    "Argument `newdata` must be a <data.frame> object\\."
  )
})

test_that("non-integer n_draws fails", {
  expect_error(
    predict(gaussian_example_fit, n_draws = data.frame()),
    "Argument `n_draws` must be a positive <integer>\\."
  )
})

test_that("negative n_draws fails", {
  expect_error(
    predict(gaussian_example_fit, n_draws = -1L),
    "Argument `n_draws` must be a positive <integer>\\."
  )
})

# Prior errors ------------------------------------------------------------

p <- get_priors(gaussian_example_fit)
f <- obs(y ~ -1 + z + varying(~ x + lag(y)), family = "gaussian") +
  random() + splines(df = 20)

test_that("incomplete priors fails", {
  p2 <- p[-1,]
  expect_error(
    dynamite(f, data = gaussian_example, time = "time", group = "id",
             priors = p2, debug = list(no_compile = TRUE)),
    paste0(
      "Argument `priors` must contain all relevant parameters:\n",
      "x Prior for parameter `sigma_nu_y` is not defined\\."
    )
  )
})

test_that("irrevelant parameters fails", {
  p2 <- rbind(p, data.frame(parameter = "extra",
                            response = "y",
                            prior = "normal(0, 1.0)",
                            type = "alpha",
                            category = ""))
  expect_error(
    dynamite(f, data = gaussian_example, time = "time", group = "id",
             priors = p2, debug = list(no_compile = TRUE)),
    paste0(
      "Argument `priors` must contain only relevant parameters:\n",
      "x Found a prior for parameter `extra` ",
      "but the model does not contain such a parameter\\."
    )
  )
})

test_that("unsupported prior distribution fails", {
  p$prior[5] <- "aaa"
  expect_error(
    dynamite(f, data = gaussian_example, time = "time", group = "id",
             priors = p, debug = list(no_compile = TRUE)),
    paste0(
      "Found an unsupported prior distribution in `priors`:\n",
      "x Distribution `aaa` is not available\\."
    )
  )
})

test_that("constrained prior for unconstrained parameter fails", {
  p$prior[5] <- "gamma(2, 1)"
  expect_error(
    dynamite(f, data = gaussian_example, time = "time", group = "id",
             priors = p, debug = list(no_compile = TRUE)),
    paste0(
      "Priors for parameters alpha, beta, and delta ",
      "should have unconstrained support:\n",
      "x Found an unconstrained distribution ",
      "`gamma` for parameter `delta_y_x`\\."
    )
  )
})

# Plot errors ----------------------------------------------------------

test_that("plot errors when the input is not a dynamitefit object", {
  expect_error(
    plot.dynamitefit(1),
    "Argument `x` must be a <dynamitefit> object."
  )
})
test_that("plot errors when no type is defined", {
  expect_error(
    plot(categorical_example_fit),
    paste0("Argument `type` is missing while it should be a single ",
    "<character> string.")
  )
})
test_that("plot errors when type is vector", {
  expect_error(
    plot(gaussian_example_fit, type = c("beta", "delta")),
    "Argument `type` must be a single <character> string."
  )
})

test_that("plot errors when no variable is found ", {
  expect_error(
    plot(categorical_example_fit, type = "delta"),
    paste0("No parameters of type `delta` found for any of the response ",
    "channels `x`, `y`.")
  )
})
test_that("plot_deltas errors when the model does not contain deltas", {
  expect_error(
    plot_deltas(categorical_example_fit),
    "The model does not contain varying coefficients delta."
  )
})
test_that("plot_nus errors when the model does not contain nus", {
  expect_error(
    plot_nus(categorical_example_fit),
    "The model does not contain random intercepts nu."
  )
})

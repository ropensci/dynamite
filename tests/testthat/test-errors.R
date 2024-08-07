data.table::setDTthreads(1) # For CRAN

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

test_that("missing formula fails", {
  expect_error(
    obs(),
    "Argument `formula` is missing\\."
  )
})

test_that("missing family fails", {
  expect_error(
    obs(y ~ x),
    "Argument `family` is missing\\."
  )
})

test_that("nonformula to dynamiteformula fails", {
  expect_error(
    obs(formula = numeric(), family = "gaussian"),
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

# test_that("simultaneity fails", {
#  obs_lhs <-
#    obs(q ~ w + e + r + lag(i), family = "gaussian") +
#    obs(t ~ y + u, family = "gaussian") +
#    obs(i ~ o + p + a + lag(f), family = "gaussian")
#  obs_rhs <-
#    obs(f ~ h + l + lag(x), family = "gaussian") +
#    obs(x ~ q + z, family = "gaussian")
#  expect_error(
#    obs_rhs + obs_lhs,
#    paste0(
#      "Simultaneous regression is not supported:\n",
#      "x Response variable `q` appears in the formula of `x`\\."
#    )
#  )
#  # should fail for deterministic as well
#  expect_error(
#    obs(y ~ x, family = "gaussian") + aux(integer(x) ~ y),
#    paste0(
#      "Simultaneous regression is not supported:\n",
#      "x Response variable `x` appears in the formula of `y`\\."
#    )
#  )
# })

test_that("cyclic dependency fails", {
  obs_lhs <- obs(y ~ x, family = "gaussian") +
    obs(z ~ y, family = "gaussian")
  obs_rhs <- aux(numeric(w) ~ z + 1) +
    obs(x ~ z, family = "gaussian")
  expect_error(
    obs_lhs + obs_rhs,
    "The model must be acyclic\\."
  )
})

test_that("contemporaneous self dependency within a channel fails", {
  expect_error(
    obs(y ~ y, family = "gaussian"),
    paste0(
      "Contemporaneous self-dependency found in model formula:\n",
      "x Variable `y` appears on both sides of the formula for \\(y\\)\\."
    )
  )
  expect_error(
    obs(c(y, x) ~ y | 1, family = "mvgaussian"),
    paste0(
      "Contemporaneous self-dependency found in model formula:\n",
      "x Variable `y` appears on both sides of the formula for \\(y, x\\)\\."
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

test_that("negative lb_tau fails", {
  expect_error(
    obs_test + splines(lb_tau = -1.0),
    "Argument `lb_tau` must be a <numeric> vector of non-negative values\\."
  )
})

test_that("time-varying definitions without splines fails", {
  obs_varying <- obs(y ~ 1 + varying(~ -1 + x), family = "gaussian")
  test_data <- data.frame(
    y = c(1, 2, 3),
    x = c(0.5, -1, 0.25),
    z = c(1, 2, 3)
  )
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
    obs_all_alpha <- obs(y1 ~ -1 + varying(~x1), family = "categorical") +
      obs(x3 ~ varying(~ -1 + x1), family = "categorical") +
      obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ lag(x3) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      obs(y9 ~ -1 + x1 + varying(~x2), family = "beta") +
      splines(df = 5, noncentered = rep(TRUE, 3)),
    NA
  )
  expect_error(
    dynamite(obs_all_alpha, test_data, "time", "group"),
    paste(
      "Length of the `noncentered` argument of `splines\\(\\)` function",
      "is not equal to 1 or 6, the number of the channels\\."
    )
  )
})

test_that("lb_tau definition throws error if not of correct length", {
  expect_error(
    obs_all_alpha <- obs(y1 ~ -1 + varying(~x1), family = "categorical") +
      obs(x3 ~ varying(~ -1 + x1), family = "categorical") +
      obs(y2 ~ -1 + x2 + varying(~1), family = "gaussian") +
      obs(y3 ~ lag(x3) + trials(trials), family = "binomial") +
      obs(y4 ~ x1 + varying(~ -1 + x2), family = "bernoulli") +
      obs(y9 ~ -1 + x1 + varying(~x2), family = "beta") +
      splines(df = 5, lb_tau = rep(1, 3)),
    NA
  )
  expect_error(
    dynamite(obs_all_alpha, test_data, "time", "group"),
    paste(
      "Length of the `lb_tau` argument of `splines\\(\\)` function is not",
      "equal to 1 or 6, the number of the channels\\."
    )
  )
})

test_that("pure deterministic formula to dynamite fails", {
  expect_error(
    dynamite(
      dformula = aux(numeric(d) ~ lag(d, 1)),
      data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
      group = "x",
      time = "z"
    ),
    "Argument `dformula` must contain at least one stochastic channel\\."
  )
})

test_that("latent factor errors with invalid responses", {
  expect_error(
    dynamite(
      obs(y ~ x, family = "gaussian") + lfactor(responses = 1),
      data = data.frame(y = rnorm(4), x = runif(4), id = 1, time = 1:4),
      time = "time",
      group = "id",
      debug = list(no_compile = TRUE)
    ),
    "Argument `responses` must be a <character> vector\\."
  )
  expect_error(
    dynamite(
      obs(y ~ x, family = "gaussian") + lfactor(responses = "x"),
      data = data.frame(y = rnorm(4), x = runif(4), id = 1, time = 1:4),
      time = "time",
      group = "id",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Argument `responses` of `lfactor\\(\\)` contains variable `x`:\n",
      "x No such response variables in the model\\."
    )
  )
})

test_that("latent factor errors with nonlogical value for nonzero_lambda", {
  expect_error(
    dynamite(
      obs(y ~ x, family = "gaussian") + lfactor(nonzero_lambda = 1),
      data = data.frame(y = rnorm(4), x = runif(4), id = 1, time = 1:4),
      time = "time",
      group = "id",
      debug = list(no_compile = TRUE)
    ),
    "Argument `nonzero_lambda` must be a <logical> vector\\."
  )
})

test_that("Random effect errors with single group", {
  expect_error(
    dynamite(
      obs(y ~ x + random(~1), family = "gaussian"),
      data = data.frame(y = rnorm(4), x = runif(4), id = 1, time = 1:4),
      time = "time",
      group = "id",
      debug = list(no_compile = TRUE)
    ),
    "Cannot estimate random effects using only one group\\."
  )
})

test_that("Latent factor errors with single group", {
  expect_error(
    dynamite(
      obs(y ~ x, family = "gaussian") + lfactor(),
      data = data.frame(y = rnorm(4), x = runif(4), id = 1, time = 1:4),
      time = "time",
      group = "id",
      debug = list(no_compile = TRUE)
    ),
    "Cannot estimate latent factors using only one group\\."
  )
})

test_that("latent factor fails with nonlogical value for noncentered_psi", {
  expect_error(
    dynamite(
      obs(y ~ x, family = "gaussian") + lfactor(noncentered_psi = 1),
      data = data.frame(y = rnorm(4), x = runif(4), id = 1, time = 1:4),
      time = "time",
      group = "id",
      debug = list(no_compile = TRUE)
    ),
    "Argument `noncentered_psi` must be a single <logical> value\\."
  )
})

test_that("latent factor fails with nonlogical value for correlated", {
  expect_error(
    dynamite(
      obs(y ~ x, family = "gaussian") + lfactor(correlated = 1),
      data = data.frame(y = rnorm(4), x = runif(4), id = 1, time = 1:4),
      time = "time",
      group = "id",
      debug = list(no_compile = TRUE)
    ),
    "Argument `correlated` must be a single <logical> value\\."
  )
})

test_that("update fails with incompatible formula", {
  expect_error(
    update(
      multichannel_example_fit,
      obs(y ~ x, family = "gaussian"),
      debug = list(no_compile = TRUE)
    ),
    "Can't find variable `x` in `data`\\."
  )
})

test_that("multivariate family fails with single response", {
  expect_error(
    obs(y1 ~ x, family = "mvgaussian"),
    "A multivariate channel must have more than one response variable\\."
  )
})

test_that("univariate family fails with multiple response variables", {
  expect_error(
    obs(c(y1, y2) ~ x, family = "gaussian"),
    "A univariate channel must have only one response variable\\."
  )
})

test_that("invalid number of formula components fails", {
  expect_error(
    obs(c(y1, y2) ~ x | x | x, family = "mvgaussian"),
    paste0(
      "Number of component formulas must be 1 ",
      "or the number of dimensions: 2\n",
      "x 3 formulas were provided\\."
    )
  )
  expect_error(
    obs(y1 ~ x | x, family = "gaussian"),
    "A univariate channel must have only one formula component\\."
  )
})

test_that("multinomial family fails with multiple formula components", {
  expect_error(
    obs(c(y1, y2, y3) ~ 1 + trials(n) | x | x, family = "multinomial"),
    "A multinomial channel must have only one formula component\\."
  )
})

test_that("cumulative channel fails without an intercept", {
  expect_error(
    obs(y ~ -1, family = "cumulative"),
    paste0(
      "A time-constant or a time-varying intercept must be specified ",
      "for a cumulative channel\\."
    )
  )
})

# Formula specials errors -------------------------------------------------

test_that("no intercept or predictors fails if no lfactor", {
  expect_error(
    dynamite(
      obs(y ~ -1, family = "gaussian"),
      data = gaussian_example,
      time = "time",
      group = "id"
    ),
    paste0(
      "Invalid formula for response variable `y`:\n",
      "x There are no predictors, intercept terms, or latent factors\\."
    )
  )
})

test_that("binomial channel without a trials term fails", {
  expect_error(
    obs(y ~ x, family = "binomial"),
    "Formula for a binomial channel must include a trials term\\."
  )
})

test_that("multinomial channel without a trials term fails", {
  expect_error(
    obs(c(y1, y2) ~ 1, family = "multinomial"),
    "Formula for a multinomial channel must include a trials term\\."
  )
})

test_that("deterministic fixed fails", {
  expect_error(
    aux(numeric(y) ~ fixed(~x)),
    paste0(
      "The use of `fixed\\(\\)` is not meaningful ",
      "for deterministic channels:\n",
      "x Time-invariant definition was found in ",
      "`numeric\\(y\\) ~ fixed\\(~x\\)`\\."
    )
  )
})

test_that("deterministic varying fails", {
  expect_error(
    aux(numeric(y) ~ varying(~x)),
    paste0(
      "The use of `varying\\(\\)` is not meaningful ",
      "for deterministic channels:\n",
      "x Time-varying definition was found in ",
      "`numeric\\(y\\) ~ varying\\(~x\\)`\\."
    )
  )
})

test_that("multiple special components fail", {
  expect_error(
    obs(y ~ fixed(~1) + fixed(~x), family = "gaussian"),
    "Multiple `fixed\\(\\)` terms are not supported\\."
  )
  expect_error(
    obs(y ~ varying(~1) + varying(~x), family = "gaussian"),
    "Multiple `varying\\(\\)` terms are not supported\\."
  )
  expect_error(
    obs(y ~ random(~1) + random(~x), family = "gaussian"),
    "Multiple `random\\(\\)` terms are not supported\\."
  )
})

test_that("specials with multiple arguments fail", {
  expect_error(
    obs(y ~ fixed(~1, 2), family = "gaussian"),
    "A `fixed\\(\\)` term must have a single formula argument\\."
  )
  expect_error(
    obs(y ~ varying(~1, 2), family = "gaussian"),
    "A `varying\\(\\)` term must have a single formula argument\\."
  )
  expect_error(
    obs(y ~ random(~1, 2), family = "gaussian"),
    "A `random\\(\\)` term must have a single formula argument\\."
  )
})

test_that("nested specials fail", {
  err <- paste0(
    "A model formula must not contain nested ",
    "`fixed\\(\\)`, `varying\\(\\)`, or `random\\(\\)` terms\\."
  )
  expect_error(
    obs(y ~ random(~1 + random(~1)), family = "gaussian"),
    err
  )
  expect_error(
    obs(y ~ varying(~1 + varying(~1)), family = "gaussian"),
    err
  )
  expect_error(
    obs(y ~ fixed(~1 + fixed(~1)), family = "gaussian"),
    err
  )
  expect_error(
    obs(y ~ random(~1 + varying(~1)), family = "gaussian"),
    err
  )
  expect_error(
    obs(y ~ varying(~1 + random(~1)), family = "gaussian"),
    err
  )
})

test_that("specials that cannot be evaluated fail", {
  expect_error(
    dynamite(
      obs(y ~ 1 + trials(log(-lag(y))), family = "binomial"),
      data = data.frame(y = 1:3, z = 1:3),
      time = "z"
    ),
    paste0(
      "Unable to evaluate `trials\\(\\)` for response variable `y`:\n",
      "x .+"
    )
  )
  expect_error(
    dynamite(
      obs(y ~ 1 + offset(log(-lag(x))), family = "poisson"),
      data = data.frame(y = 1:3, z = 1:3),
      time = "z"
    ),
    paste0(
      "Unable to evaluate `offset\\(\\)` for response variable `y`:\n",
      "x .+"
    )
  )
})

test_that("test that specials with invalid values fail", {
  expect_error(
    dynamite(
      obs(y ~ 1 + trials(n), family = "binomial"),
      data = data.frame(y = 1:3, z = 1:3, n = factor(1:3)),
      time = "z"
    ),
    paste0(
      "Invalid `trials\\(\\)` definition for response variable `y`:\n",
      "x Number of trials cannot be a <factor>\\."
    )
  )
  expect_error(
    dynamite(
      obs(y ~ 1 + trials(n), family = "binomial"),
      data = data.frame(y = 1:3, z = 1:3, n = -(1:3)),
      time = "z"
    ),
    paste0(
      "Invalid `trials\\(\\)` definition for response variable `y`:\n",
      "x Number of trials must contain only positive integers\\."
    )
  )
  expect_error(
    dynamite(
      obs(y ~ 1 + offset(n), family = "poisson"),
      data = data.frame(y = 1:3, z = 1:3, n = factor(1:3)),
      time = "z"
    ),
    paste0(
      "Invalid `offset\\(\\)` definition for response variable `y`:\n",
      "x Offset cannot be a <factor>\\."
    )
  )
})

# Data errors -------------------------------------------------------------

test_that("missing data object fails", {
  expect_error(
    dynamite(dformula = obs_test),
    "Argument `data` is missing\\."
  )
})

test_that("missing time variable fails", {
  expect_error(
    dynamite(dformula = obs_test, data = data.frame(z = 1), group = "z"),
    "Argument `time` is missing\\."
  )
})

test_that("data is not data.frame fails", {
  expect_error(
    dynamite(dformula = obs_test, data = list(), time = "z"),
    "Argument `data` must be a <data.frame> object\\."
  )
})

test_that("group variable not in data fails", {
  expect_error(
    dynamite(
      dformula = obs_test,
      data = data.frame(y = 1, x = 1), time = "x", group = "z"
    ),
    "Can't find grouping variable `z` in `data`\\."
  )
})

test_that("time variable not in data fails", {
  expect_error(
    dynamite(
      dformula = obs_test,
      data = data.frame(y = 1, x = 1),
      time = "z"
    ),
    "Can't find time index variable `z` in `data`\\."
  )
})

test_that("single time point fails", {
  expect_error(
    dynamite(
      dformula = obs_test,
      data = data.frame(y = 1, x = 1, z = 1),
      time = "z",
      group = "x"
    ),
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
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
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
      time = "z",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Each time index must correspond to a single observation per group:\n",
      "x Group `1` of `.group` has duplicate observations\\."
    )
  )
})

test_that("missing lag variable fails", {
  expect_error(
    dynamite(
      dformula = obs(y ~ lag(d, 1), family = "gaussian"),
      data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Unable to construct lagged values of `d`:\n",
      "x Can't find such variables in `data`\\."
    )
  )
})

test_that("missing predictor fails", {
  expect_error(
    dynamite(
      dformula = obs(y ~ w, family = "gaussian"),
      data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
    ),
    "Can't find variable `w` in `data`\\."
  )
})

test_that("invalid deterministic channel definition fails", {
  expect_error(
    dynamite(
      dformula = obs(y ~ x, family = "gaussian") +
        aux(integer(d) ~ 1 + w),
      data = data.frame(y = c(1, 1), x = c(1, 1), z = c(1, 2)),
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Unable to evaluate definitions of deterministic channels:\n",
      "x object 'w' not found"
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

# Data type errors --------------------------------------------------------

#' @srrstats {G2.11, G2.12} Tests for unsupported column types.
test_that("invalid column types fail", {
  test_data <- data.frame(y = c(1i, 2i), x = c(1, 1), z = c(1, 2))
  test_data$w <- c(list(a = 1), list(b = 2))
  test_data$d <- as.raw(c(40, 20))
  expect_error(
    dynamite(
      dformula = obs(y ~ x, family = "gaussian"),
      data = test_data,
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Columns `y`, `w`, and `d` of `data` are invalid:\n",
      "x Column types <complex/list/raw> are not supported\\."
    )
  )
})

test_that("non-finite values in data fail", {
  test_data <- data.frame(
    y = c(1, Inf), x = c(1, 1),
    z = c(1, 2), w = c(-Inf, 2), u = c(1, Inf)
  )
  expect_error(
    dynamite(
      dformula = obs(y ~ x, family = "gaussian"),
      data = test_data,
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
    ),
    "Non-finite values were found in variables `y`, `w`, and `u` of `data`\\."
  )
})

test_that("non-factor categorical response fails", {
  test_data <- data.frame(y = c(0, 1), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(
      dformula = obs(y ~ 1, family = "categorical"),
      data = test_data,
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Response variable `y` is invalid:\n",
      "x Categorical family supports only <factor> variables\\."
    )
  )
})

test_that("factor types for non-categorical families fails", {
  test_data <- data.frame(
    y = factor(c(0, 1)),
    w = c(1, 2),
    x = c(1, 1),
    z = c(1, 2)
  )
  families <- c(
    "gaussian",
    "exponential",
    "gamma",
    "beta",
    "bernoulli",
    "binomial",
    "poisson",
    "negbin",
    "student"
  )
  for (f in families) {
    form <- ifelse_(identical(f, "binomial"), y ~ 1 + trials(x), y ~ 1)
    expect_error(
      dynamite(
        dformula = obs(form, family = f),
        data = test_data,
        time = "z",
        group = "x",
        debug = list(no_compile = TRUE)
      ),
      paste0(
        "Response variable `y` is invalid:\n",
        "x .+ family is not supported for <factor> variables\\."
      )
    )
  }
  mvfamilies <- c(
    "mvgaussian",
    "multinomial"
  )
  for (f in mvfamilies) {
    form <- ifelse_(
      identical(f, "multinomial"),
      c(y, w) ~ 1 + trials(x),
      c(y, w) ~ 1
    )
    expect_error(
      dynamite(
        dformula = obs(form, family = f),
        data = test_data,
        time = "z",
        group = "x",
        debug = list(no_compile = TRUE)
      ),
      paste0(
        "Response variable `.+` is invalid:\n",
        "x .+ family is not supported for <factor> variables\\."
      )
    )
  }
})

test_that("negative values for distributions with positive support fails", {
  test_data <- data.frame(y = c(-1, -2), w = c(1, 2), x = c(1, 1), z = c(1, 2))
  families <- c(
    "exponential",
    "gamma",
    "binomial",
    "negbin",
    "poisson"
  )
  for (f in families) {
    form <- ifelse_(identical(f, "binomial"), y ~ 1 + trials(x), y ~ 1)
    expect_error(
      dynamite(
        dformula = obs(form, family = f),
        data = test_data,
        time = "z",
        group = "x",
        debug = list(no_compile = TRUE)
      ),
      paste0(
        "Response variable `y` is invalid:\n",
        "x .+ family supports only non-negative .+\\."
      )
    )
  }
  expect_error(
    dynamite(
      dformula = obs(c(y, w) ~ 1 + trials(c(2, 3)), family = "multinomial"),
      data = test_data,
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Response variable `y_w` is invalid:\n",
      "x Multinomial family supports only non-negative .+\\."
    )
  )
})

test_that("bernoulli without 0/1 values fails", {
  test_data <- data.frame(y = c(2, 3), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(
      dformula = obs(y ~ 1, family = "bernoulli"),
      data = test_data,
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Response variable `y` is invalid:\n",
      "x Bernoulli family supports only 0/1 integers\\."
    )
  )
})

test_that("beta without (0, 1) values fails", {
  test_data <- data.frame(y = c(2, 3), x = c(1, 1), z = c(1, 2))
  expect_error(
    dynamite(
      dformula = obs(y ~ 1, family = "beta"),
      data = test_data,
      time = "z",
      group = "x",
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Response variable `y` is invalid:\n",
      "x Beta family supports only values on the open interval \\(0, 1\\)\\."
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
    complete_lags(quote(lag(y, "a"))),
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

test_that("output for missing argument fails", {
  methods <- c(
    "as.data.frame",
    "as_draws_df",
    "confint",
    "coef",
    "fitted",
    "formula",
    "hmc_diagnostics",
    "lfo",
    "loo",
    "mcmc_diagnostics",
    "ndraws",
    "nobs",
    "plot",
    "predict",
    "print",
    "summary",
    "update"
  )
  for (m in methods) {
    call_fun <- paste0(m, ".dynamitefit")
    expect_error(
      do.call(call_fun, args = list()),
      "Argument `.+` is missing"
    )
  }
})

test_that("output for non dynamitefit objects fails", {
  methods <- c(
    "as.data.frame",
    "as_draws_df",
    "confint",
    "coef",
    "fitted",
    "formula",
    "hmc_diagnostics",
    "lfo",
    "loo",
    "mcmc_diagnostics",
    "ndraws",
    "nobs",
    "plot",
    "predict",
    "print",
    "summary",
    "update"
  )
  object_arg_methods <- c(
    "coef",
    "confint",
    "fitted",
    "nobs",
    "predict",
    "summary",
    "update"
  )
  for (m in methods) {
    args <- ifelse_(
      m %in% object_arg_methods,
      list(object = 1L),
      list(x = 1L)
    )
    call_fun <- paste0(m, ".dynamitefit")
    expect_error(
      do.call(call_fun, args = args),
      "Argument `.+` must be a <dynamitefit> object\\."
    )
  }
})

test_that("output without Stan fit fails", {
  methods <- c(
    "as.data.frame",
    "as_draws_df",
    "fitted",
    "lfo",
    "loo",
    "predict",
    "ndraws"
  )
  object_arg_methods <- c(
    "fitted",
    "predict"
  )
  fit <- gaussian_example_fit
  fit$stanfit <- NULL
  for (m in methods) {
    args <- ifelse_(
      m %in% object_arg_methods,
      list(object = fit),
      list(x = fit)
    )
    expect_error(
      do.call(paste0(!!m, ".dynamitefit"), args = args),
      "No Stan model fit is available\\."
    )
  }
})

test_that("invalid responses fail", {
  expect_error(
    as.data.table(gaussian_example_fit, responses = "resp"),
    paste0(
      "Argument `responses` contains invalid response variable names\\.\n",
      "x Response variable \"resp\" is not recognized\\.\n",
      "i The response variable of the model is \"y\"\\."
    )
  )
})

test_that("invalid parameters fail", {
  expect_error(
    as.data.table(gaussian_example_fit, parameter = "test"),
    paste0(
      "Argument `parameters` contains invalid parameter names\\.\n",
      "x Parameter \"test\" is not recognized\\.\n",
      "i Use `get_parameter_names\\(\\)` to check available parameters\\."
    )
  )
})

test_that("invalid types fail", {
  expect_error(
    as.data.table(gaussian_example_fit, types = c("aa", "bb")),
    paste0(
      "Argument `types` contains invalid types\\.\n",
      "x Types \"aa\" and \"bb\" are not recognized\\.\n",
      "i Use `get_parameter_types\\(\\)` to check available types\\."
    )
  )
})

test_that("not found parameters fail", {
  expect_error(
    as.data.table(categorical_example_fit, types = "delta"),
    paste0(
      "No parameters of type `delta` were found for any of the response ",
      "channels `x` and `y`\\."
    )
  )
})

test_that("invalid confint level fails", {
  expect_error(
    confint.dynamitefit(gaussian_example_fit, level = -0.1),
    "Argument `level` must be a single <numeric> value between 0 and 1\\."
  )
})

test_that("Invalid code blocks fail", {
  expect_error(
    get_code(gaussian_example_fit, blocks = mean),
    "Argument `blocks` must be a <character> vector or NULL\\."
  )
  expect_error(
    get_code(gaussian_example_fit, blocks = "block"),
    paste0(
      "Invalid Stan blocks provided: block\n",
      "i Argument `blocks` must be NULL or a subset of .*"
    )
  )
})

# Predict errors ----------------------------------------------------------

gaussian_example_small <- gaussian_example |> dplyr::filter(.data$time < 6)

# test_that("newdata without group variable fails when there are groups", {
#  gaussian_example_nogroup <- gaussian_example_small |>
#    dplyr::select(!"id")
#  expect_error(
#    predict(gaussian_example_fit, newdata = gaussian_example_nogroup),
#    "Can't find grouping variable `id` in `newdata`\\."
#  )
# })

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
    dplyr::select(!"time")
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
      dplyr::select(!"id"),
    data.frame(y = 1, x = 1, z = 0, time = 1)
  )
  expect_error(
    predict(gaussian_example_single_fit, newdata = gaussian_example_duplicated),
    paste0(
      "Each time index must correspond to a single observation per group:\n",
      "x Group `1` of `.group` has duplicate observations\\."
    )
  )
})

test_that("new group levels can't be included if new_levels is 'none'", {
  gaussian_example_new_levels <- rbind(
    gaussian_example,
    data.frame(
      y = c(0.5, rep(NA, 29L)),
      x = rnorm(30),
      z = rbinom(30, 1, 0.7),
      id = 226L, time = seq.int(1, 30)
    )
  )
  expect_error(
    predict(
      gaussian_example_fit,
      newdata = gaussian_example_new_levels,
      type = "response", n_draws = 2, new_levels = "none"
    ),
    paste(
      "Grouping variable `id` contains unknown levels:\nx Level \"226\"",
      "is not present in the original data\\.\ni Note: argument `new_levels`",
      "is \"none\" which disallows new levels\\."
    )
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
  gaussian_example_misresp <- gaussian_example_small |> dplyr::select(!"y")
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

test_that("non-logical expand fails", {
  expect_error(
    predict(gaussian_example_fit, expand = data.frame()),
    "Argument `expand` must be a single <logical> value\\."
  )
})

test_that("invalid funs fails", {
  expect_error(
    predict(gaussian_example_fit, funs = 1L),
    "Argument `funs` must be a <list>\\."
  )
  expect_error(
    predict(gaussian_example_fit, funs = list(1L)),
    "Argument `funs` must be named\\."
  )
  expect_error(
    predict(gaussian_example_fit, funs = list(w = 1L)),
    "The names of `funs` must be response variables of the model\\."
  )
  expect_error(
    predict(gaussian_example_fit, funs = list(y = 1L)),
    "Each element of `funs` must be a <list>\\."
  )
  expect_error(
    predict(gaussian_example_fit, funs = list(y = list(1L))),
    "Each element of `funs` must be named\\."
  )
  expect_error(
    predict(gaussian_example_fit, funs = list(y = list(fun = 1L))),
    "Each element of `funs` must contain only functions\\."
  )
})

# Prior errors ------------------------------------------------------------

p <- get_priors(gaussian_example_fit)
f <- obs(y ~ -1 + random(~1) + z + varying(~ x + lag(y)), family = "gaussian") +
  splines(df = 20)

test_that("incomplete priors fails", {
  p2 <- p[-1, ]
  expect_error(
    dynamite(
      f,
      data = gaussian_example,
      time = "time",
      group = "id",
      priors = p2,
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Argument `priors` must contain all relevant parameters:\n",
      "x Prior for parameter `sigma_nu_y_alpha` is not defined\\."
    )
  )

  expect_error(
    update(gaussian_example_fit,
      priors = p2,
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Argument `priors` must contain all relevant parameters:\n",
      "x Prior for parameter `sigma_nu_y_alpha` is not defined\\."
    )
  )
})

test_that("irrevelant parameters fails", {
  p2 <- rbind(p, data.frame(
    parameter = "extra",
    response = "y",
    prior = "normal(0, 1.0)",
    type = "alpha",
    category = ""
  ))
  expect_error(
    dynamite(
      f,
      data = gaussian_example,
      time = "time",
      group = "id",
      priors = p2,
      debug = list(no_compile = TRUE)
    ),
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
    dynamite(
      f,
      data = gaussian_example,
      time = "time",
      group = "id",
      priors = p,
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Found an unsupported prior distribution in `priors`:\n",
      "x Distribution `aaa` is not available\\."
    )
  )
})

test_that("constrained prior for unconstrained parameter fails", {
  p$prior[5] <- "gamma(2, 1)"
  expect_error(
    dynamite(
      f,
      data = gaussian_example,
      time = "time",
      group = "id",
      priors = p,
      debug = list(no_compile = TRUE)
    ),
    paste0(
      "Priors for parameters \"alpha\", \"beta\", and \"delta\" ",
      "should have unconstrained support:\n",
      "x Found an unconstrained distribution ",
      "`gamma` for parameter `delta_y_x`\\."
    )
  )
})

# Plot errors ----------------------------------------------------------

test_that("plot errors when the input is not a dynamitefit object", {
  expect_error(
    plot.dynamitefit(1, types = "beta"),
    "Argument `x` must be a <dynamitefit> object."
  )
})

# Model errors ------------------------------------------------------------

test_that("multinomial model fails if stan version < 2.24", {
  set.seed(1)
  n_id <- 10L
  n_time <- 5L
  d <- data.frame(
    y1 = sample(10, size = n_id * n_time, replace = TRUE),
    y2 = sample(15, size = n_id * n_time, replace = TRUE),
    y3 = sample(20, size = n_id * n_time, replace = TRUE),
    z = rnorm(n_id * n_time),
    time = seq_len(n_time),
    id = rep(seq_len(n_id), each = n_time)
  )
  d$n <- d$y1 + d$y2 + d$y3
  f <- obs(
    c(y1, y2, y3) ~ z + lag(y1) + lag(y2) + lag(y3) + trials(n),
    family = "multinomial"
  )
  expect_error(
    mockthat::with_mock(
      stan_version = function(...) "2.23",
      dynamite(
        dformula = f,
        data = d,
        time = "time",
        group = "id",
        backend = "rstan"
      )
    ),
    paste0(
      "Multinomial family is not supported for this version of rstan\\.\n",
      "i Please install a newer version of rstan\\."
    )
  )
})

# Stan errors -------------------------------------------------------------

test_that("Stan backend argument conversion duplicates fail", {
  dots <- list(iter = 1000, iter_sampling = 1000)
  expect_error(
    check_stan_args(dots, verbose = FALSE, backend = "rstan"),
    paste0(
      "Conflict in argument syntax conversion from cmdstanr to rstan\\.\n",
      "x Argument `iter` has been multiply specified\\."
    )
  )
  dots <- list(
    iter = 1000,
    iter_sampling = 1000,
    cores = 3,
    parallel_chains = 4
  )
  expect_error(
    check_stan_args(dots, verbose = FALSE, backend = "cmdstanr"),
    paste0(
      "Conflict in argument syntax conversion from rstan to cmdstanr\\.\n",
      "x Arguments `iter_sampling` and `parallel_chains` have been multiply ",
      "specified\\."
    )
  )
})

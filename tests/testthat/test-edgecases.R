#' @srrstats {G5.8} **Edge condition tests** *to test that these conditions produce expected behaviour such as clear warnings or errors when confronted with data with extreme properties including but not limited to:*
#' @srrstats {G5.8a} *Zero-length data*
#' @srrstats {G5.8b} *Data of unsupported types (e.g., character or complex numbers in for functions designed only for numeric data)*
#' @srrstats {G5.8c} *Data with all-`NA` fields or columns or all identical fields or columns*
#' @srrstats {G5.8d} *Data outside the scope of the algorithm (for example, data with more fields (columns) than observations (rows) for some regression algorithms)*

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
    x1 = sample(letters[1:4], size = total_obs, replace = TRUE),
    x2 = rnorm(total_obs),
    x3 = as.factor(sample(4, size = total_obs, replace = TRUE))
  )

debug <- list(no_compile = TRUE)

test_that("single channel models are valid", {
  expect_error(
    obs_categorical <- obs(y1 ~ x1, family = categorical()), NA
  )
  expect_error(
    obs_gaussian <- obs(y2 ~ x2, family = gaussian()), NA
  )
  expect_error(
    obs_binomial <- obs(y3 ~ x3 + trials(trials), family = binomial()), NA
  )
  expect_error(
    obs_bernoulli <- obs(y4 ~ x1, family = bernoulli()), NA
  )
  expect_error(
    obs_negbinom <- obs(y5 ~ x2, family = negbin()), NA
  )
  expect_error(
    obs_poisson <- obs(y6 ~ x3 + offset(offset), family = poisson()), NA
  )
  expect_error(
    obs_exponential <- obs(y7 ~ x3, family = exponential()), NA
  )
  expect_error(
    obs_gamma <- obs(y8 ~ x1, family = gamma()), NA
  )
  expect_error(
    dynamite(obs_categorical, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_gaussian, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_binomial, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_bernoulli, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_negbinom, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_poisson, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_exponential, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_gamma, test_data, "group", "time", debug = debug), NA
  )
})

test_that("multichannel models are valid", {
  expect_error(
    obs_all <- obs(y1 ~ x1, family = categorical()) +
      obs(y2 ~ x2, family = gaussian()) +
      obs(y3 ~ x3 + trials(trials), family = binomial()) +
      obs(y4 ~ x1, family = bernoulli()) +
      obs(y5 ~ x2, family = negbin()) +
      obs(y6 ~ x3 + offset(offset), family = poisson()) +
      obs(y7 ~ x3, family = exponential()) +
      obs(y8 ~ x1, family = gamma()),
    NA
  )
  expect_error(
    dynamite(obs_all, test_data, "group", "time", debug = debug), NA
  )
})

test_that("intercepts are handled correctly", {
  expect_error(
    obs_all_alpha <- obs(y1 ~ -1 + x1, family = categorical()) +
      obs(y2 ~ -1 + x2 + varying(~1), family = gaussian()) +
      obs(y3 ~ -1 + x3 + varying(~x1) + trials(trials), family = binomial()) +
      obs(y4 ~ x1 + varying(~-1 + x2), family = bernoulli()) +
      splines(df = 5),
    NA
  )
  expect_error(
    dynamite(obs_all_alpha, test_data, "group", "time", debug = debug), NA
  )
})

test_that("lags are parsed", {
  expect_error(
    obs_a <- obs(y1 ~ x1 + lag(y2, 1), family = categorical()) +
      obs(y2 ~ x2 + lag(y1, 1), family = gaussian()),
    NA
  )
  expect_error(
    obs_b <- obs(y1 ~ -1 + x1 + varying(~ lag(y2, 1)), family = categorical()) +
      obs(y2 ~ -1 + x2 + varying(~lag(y1, 1)), family = gaussian()) +
      splines(),
    NA
  )
  expect_error(
    obs_c <- obs(y1 ~ x1, family = categorical()) +
      obs(y2 ~ x2, family = gaussian()) +
      lags(k = 1, type = "fixed"),
    NA
  )
  expect_error(
    obs_d <- obs(y1 ~ x1, family = categorical()) +
      obs(y2 ~ x2, family = gaussian()) +
      lags(k = 1, type = "varying") +
      splines(),
    NA
  )
  expect_error(
    dynamite(obs_a, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_b, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_c, test_data, "group", "time", debug = debug), NA
  )
  expect_error(
    dynamite(obs_d, test_data, "group", "time", debug = debug), NA
  )
})

test_that("lags() and lag() give equal results", {
  f1 <- expect_error(
    obs(y1 ~ x1, family = categorical()) +
      obs(y2 ~ x2, family = gaussian()) +
      lags(k = 1, type = "fixed"),
    NA)
  f2 <- expect_error(
    obs(y1 ~ x1 + lag(y1) + lag(y2), family = categorical()) +
      obs(y2 ~ x2 + lag(y1) + lag(y2), family = gaussian()),
    NA)
  expect_identical(
    get_priors(f1, test_data, "group", "time"),
    get_priors(f2, test_data, "group", "time")
  )
})


test_that("higher order lags() and lag() give equal results", {
  f1 <- obs(y1 ~ x1, family = categorical()) +
    obs(y2 ~ x2, family = gaussian()) +
    lags(k = 1:2, type = "fixed")
  f2 <- obs(y1 ~ x1 + lag(y1, 1) + lag(y2, 1) +
              lag(y1, 2) + lag(y2, 2), family = categorical()) +
    obs(y2 ~ x2 + lag(y1, 1) + lag(y2, 1) +
          lag(y1, 2) + lag(y2, 2), family = gaussian())
  expect_identical(
    get_priors(f1, test_data, "group", "time"),
    get_priors(f2, test_data, "group", "time")
  )
})

test_that("vector of lags and lag of vector give equal resuls", {
  f1 <- obs(y1 ~ x1 + lag(y1, 1:2), family = categorical())
  f2 <- obs(y1 ~ x1 + lag(y1, c(1, 2)), family = categorical())
  f3 <- obs(y1 ~ x1 + lag(y1, seq(1, 2)), family = categorical())
  f4 <- obs(y1 ~ x1 + lag(y1, 1) + lag(y1, 2), family = categorical())
  expect_identical(
    get_priors(f1, test_data, "group", "time"),
    get_priors(f2, test_data, "group", "time")
  )
  expect_identical(
    get_priors(f1, test_data, "group", "time"),
    get_priors(f3, test_data, "group", "time")
  )
  expect_identical(
    get_priors(f1, test_data, "group", "time"),
    get_priors(f4, test_data, "group", "time")
  )
})

test_that("deterministic channels are parsed", {
  expect_error(
    obs_det <- obs(y5 ~ x1 + lag(d, 1) + lag(y5, 1) + lag(x1, 1), family = negbin()) +
      aux(numeric(d) ~ lag(d, 1) + lag(f, 2) + x2 + past(0)) +
      aux(numeric(f) ~ lag(y5, 1) + x2 * 3 + 1 + past(0, 1)),
    NA
  )
  expect_error(
    dynamite(obs_det, test_data, "group", "time", debug = debug), NA
  )
})

test_that("deterministic simultaneity is supported", {
  expect_error(
      obs(y5 ~ x1 + lag(d, 1) + lag(y5, 1) + lag(x1, 1), family = negbin()) +
      aux(numeric(d) ~ y5 + 3),
    NA
  )
})

test_that("deterministic types are supported", {
  expect_error(
    aux(factor(a) ~ factor(c(1,2,3), levels = c(1,2,3))) +
    aux(numeric(b) ~ log(1.0)) +
    aux(integer(c) ~ 1L) +
    aux(logical(d) ~ TRUE),
    NA
  )
})

test_that("manual fixed() terms work", {
  expect_error(
    obs_fixed <- obs(y1 ~ fixed(~x1 + lag(y2, 1)), family = categorical()),
    NA
  )
})

test_that("deterministic lags with zero observed lags is evaluated", {
  obs_zerolag <-
    obs(y2 ~ x1, family = gaussian()) +
    aux(numeric(d) ~ abs(y2) + lag(d) + past(0.5))
  expect_error(
    dynamite(obs_zerolag, test_data, "group", "time", debug = debug), NA
  )
})

test_that("data expansion to full time scale works", {
  set.seed(1)
  mis_rows <- sample(1:nrow(test_data), 10)
  test_data_mis <- test_data[-mis_rows,]
  fit <- dynamite(obs(y2 ~ x2, family = gaussian()),
                  test_data_mis, "group", "time", debug = debug)
  expected_data <- test_data
  expected_data[mis_rows,3:ncol(test_data)] <- NA
  expected_data$x1 <- factor(expected_data$x1)
  data.table::setDT(expected_data, key = c("group", "time"))
  expect_equal(fit$data, expected_data, ignore_attr = TRUE)
})

test_that("no groups data preparation works", {
  test_data_single <- test_data |>
    dplyr::filter(.data$group == 1) |>
    dplyr::select(!.data$group)
  obs_all <- obs(y1 ~ x1, family = categorical()) +
    obs(y2 ~ x2, family = gaussian()) +
    obs(y3 ~ x3 + trials(trials), family = binomial()) +
    obs(y4 ~ x1, family = bernoulli()) +
    obs(y5 ~ x2, family = negbin()) +
    obs(y6 ~ x3 + offset(offset), family = poisson()) +
    obs(y7 ~ x3, family = exponential()) +
    obs(y8 ~ x1, family = gamma())
  expect_error(
    dynamite(obs_all, test_data_single, time = "time", debug = debug),
    NA
  )
})

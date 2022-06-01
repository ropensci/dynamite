timepoints <- 10
individuals <- 5
total_obs <- timepoints * individuals

test_data <- data.frame(
  offset = sample(50:100, size = total_obs, replace = TRUE),
  trials = sample(50:100, size = total_obs, replace = TRUE),
  group = gl(individuals, timepoints),
  time = 1:timepoints
) |>
  dplyr::mutate(
    y1 = as.factor(sample(5, size = total_obs, replace = TRUE)),
    y2 = rnorm(n = total_obs, mean = 1, sd = 2),
    y3 = rbinom(n = total_obs, size = trials, prob = 0.75),
    y4 = rbinom(n = total_obs, size = 1, prob = 0.66),
    y5 = rnbinom(n = total_obs, size = 100, prob = 0.33),
    y6 = rpois(n = total_obs, lambda = log(offset) + 1),
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
})

test_that("Multichannel models are valid", {
  expect_error(
    obs_all <- obs(y1 ~ x1, family = categorical()) +
      obs(y2 ~ x2, family = gaussian()) +
      obs(y3 ~ x3 + trials(trials), family = binomial()) +
      obs(y4 ~ x1, family = bernoulli()) +
      obs(y5 ~ x2, family = negbin()) +
      obs(y6 ~ x3 + offset(offset), family = poisson()),
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

test_that("Lags are parsed", {
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

test_that("Deterministic channels are parsed", {
  expect_error(
    obs_det <- obs(y5 ~ x1 + lag(d, 1) + lag(y5, 1) + lag(x1, 1), family = negbin()) +
      aux(d ~ lag(d, 1) + lag(f, 2) + x2 + past(0, 0)) +
      aux(f ~ lag(y5, 1) + x2 * 3 + 1 + past(0, 1, 2)),
    NA
  )
  expect_error(
    dynamite(obs_det, test_data, "group", "time", debug = debug), NA
  )
})

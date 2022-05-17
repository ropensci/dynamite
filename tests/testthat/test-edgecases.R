testthat::test_that("single channel models are valid", {
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

  debug <- list(no_compile = TRUE)

  expect_error(
    dynamite(obs_categorical, test_data, group, time, debug = debug), NA
  )
  expect_error(
    dynamite(obs_gaussian, test_data, group, time, debug = debug), NA
  )
  expect_error(
    dynamite(obs_binomial, test_data, group, time, debug = debug), NA
  )
  expect_error(
    dynamite(obs_bernoulli, test_data, group, time, debug = debug), NA
  )
  expect_error(
    dynamite(obs_negbinom, test_data, group, time, debug = debug), NA
  )
  expect_error(
    dynamite(obs_poisson, test_data, group, time, debug = debug), NA
  )


})

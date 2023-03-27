run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

test_that("multivariate gaussian fit and predict work", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  N <- 100
  T_ <- 50
  S <- crossprod(matrix(rnorm(4), 2, 2))
  L <- t(chol(S))
  y1 <- matrix(0, N, T_)
  y2 <- matrix(0, N, T_)
  x <- matrix(0, N, T_)
  for (t in 2:T_) {
    for (i in 1:N){
      mu <- c(0.7 * y1[i, t-1], 0.4 * y2[i, t-1] - 0.2 * y1[i, t-1])
      y <- mu + L %*% rnorm(2)
      y1[i, t] <- y[1L]
      y2[i, t] <- y[2L]
      x[i, t] <- rnorm(1, c(0.5 * y1[i, t-1]), 0.5)
    }
  }
  d <- data.frame(
    y1 = c(y1),
    y2 = c(y2),
    x = c(x),
    t = rep(1:T_, each = N),
    id = 1:N
  )
  f <- obs(c(y1, y2) ~ -1 + lag(y1) | -1 + lag(y1) + lag(y2), "mvgaussian") +
    obs(x ~ -1 + lag(y1), "gaussian")
  fit <- dynamite(
    dformula = f,
    data = d,
    time = "t",
    group = "id",
    chains = 1,
    iter = 2000,
    refresh = 0
  )
  expect_error(sumr <- summary(fit, type = "corr"), NA)
  expect_equal(sumr$mean, cov2cor(S)[2, 1], tolerance = 0.1)
  expect_error(sumr <- summary(fit, type = "sigma"), NA)
  expect_equal(sumr$mean, c(0.5, sqrt(diag(S))), tolerance = 0.1)
  expect_error(sumr <- summary(fit, type = "beta"), NA)
  expect_equal(sumr$mean, c(0.5, 0.7, -0.2, 0.4), tolerance = 0.1)
  expect_error(predict(fit, n_draws = 5), NA)
})

test_that("multinomial fit and predict work", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  n_id <- 100L
  n_time <- 20L

  d <- data.frame(
    y1 = sample(10, size = n_id, replace = TRUE),
    y2 = sample(12, size = n_id, replace = TRUE),
    y3 = sample(15, size = n_id, replace = TRUE),
    time = 1,
    id = seq_len(n_id)
  )
  d$n <- d$y1 + d$y2 + d$y3

  d <- dplyr::right_join(
    d,
    data.frame(
      time = rep(seq_len(n_time), each = n_id),
      id = seq_len(n_id)
    )
  )

  d$z <- rnorm(nrow(d))
  d$n[is.na(d$n)] <- sample(37, size = sum(is.na(d$n)), replace = TRUE)

  f <- obs(
    c(y1, y2, y3) ~ z + lag(y1) + lag(y2) + lag(y3) + trials(n),
    family = "multinomial"
  )

  init <- list(
    beta_y2 = c(1.2, 0.8, 0.2, 0.1),
    beta_y3 = c(1, 0.5, 0.6, 0.3),
    a_y2 = -0.1,
    a_y3 = 0.2
  )

  expect_error(
    fit <- dynamite(
      dformula = f,
      data = d,
      time = "time",
      group = "id",
      chains = 1,
      iter = 1,
      algorithm = "Fixed_param",
      init = list(init),
    ),
    NA
  )
  expect_error(
    pred <- predict(fit, type = "response"),
    NA
  )
  expect_identical(
    pred$y1_new + pred$y2_new + pred$y3_new,
    pred$n
  )
})

test_that("non-glm categorical fit works", {
  skip_if_not(run_extended_tests)

  expect_error(
    mockthat::with_mock(
      stan_supports_categorical_logit_glm = function(...) FALSE,
      dynamite(
        dformula = obs(x ~ z + lag(x) + lag(y), family = "categorical") +
          obs(y ~ z + lag(x) + lag(y), family = "categorical"),
        data = categorical_example,
        time = "time",
        group = "id",
        chains = 1,
        refresh = 0,
        thin = 5,
        verbose = FALSE,
        save_warmup = FALSE
      )
    ),
    NA
  )
})

test_that("get_parameter_dims works for dynamiteformula", {
  skip_if_not(run_extended_tests)

  expect_error(
    dims <- get_parameter_dims(
      x = obs(x ~ z + lag(x) + lag(y), family = "categorical") +
        obs(y ~ z + lag(x) + lag(y), family = "categorical"),
      data = categorical_example,
      time = "time",
      group = "id"
    ),
    NA
  )
  expect_identical(
    dims,
    get_parameter_dims(categorical_example_fit)
  )
})

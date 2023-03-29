#' @srrstats {G5.10} Extended tests can be switched on via setting the
#'   environment variable DYNAMITE_EXTENDED_TESTS to "true".

run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

test_that("multiple random effects work in fit and predict", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  n <- 40
  k <- 10
  x <- rnorm(n * k)
  u1 <- rep(rnorm(k, sd = 0.2), each = n)
  u2 <- 0.5 * u1 + rep(rnorm(k, sd = 0.1), each = n)
  u3 <- 0.2 * u1 + rep(rnorm(k, sd = 0.3), each = n)
  y1 <- rbinom(n * k, size = 20, prob = plogis(x + u1 + u2 * x))
  y2 <- rnorm(n * k, u3 + 2 * x)
  d <- data.frame(
    year = 1:n, person = rep(1:k, each = n),
    y1 = y1, y2 = y2, x = x, tr = 20
  )

  expect_error(
    fit <- dynamite(
      obs(y1 ~ x + trials(tr) + random(~x), family = "binomial") +
        obs(y2 ~ x + random(~1), family = "gaussian") +
        random_spec(),
      data = d, time = "year", group = "person",
      chains = 1, iter = 2000, refresh = 0
    ),
    NA
  )

  expect_error(
    sumr <- summary(fit, type = "sigma_nu"),
    NA
  )
  expect_equal(sumr$mean, c(0.2338, 0.1331, 0.3339), tolerance = 0.2)

  expect_error(
    sumr <- summary(fit, type = "corr_nu"),
    NA
  )
  expect_equal(sumr$mean, c(0.336, -0.0522, -0.0712), tolerance = 0.2)

  expect_error(
    predict(fit, n_draws = 5),
    NA
  )

  newdata <- rbind(
    d[(n + 1):nrow(d), ], # remove one person and add two
    data.frame(
      y1 = c(4, rep(NA, n - 1), 4, rep(NA, n - 1)),
      y2 = c(0, rep(NA, n - 1), 0.5, rep(NA, n - 1)),
      x = rnorm(2 * n),
      person = rep(c(226L, 500L), each = n),
      year = seq.int(1, n),
      tr = rep(c(10, 25), each = n)
    )
  )
  expect_error(
    predict(
      fit,
      newdata = newdata,
      type = "response",
      n_draws = 5,
      new_levels = "bootstrap"
    ),
    NA
  )
  expect_error(
    predict(
      fit,
      newdata = newdata,
      type = "response",
      n_draws = 5,
      new_levels = "gaussian"
    ),
    NA
  )
  expect_error(
    predict(
      fit,
      newdata = newdata,
      type = "response",
      n_draws = 2,
      new_levels = "original"
    ),
    NA
  )
})

test_that("centered and noncentered parameterization for random effects work", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  n <- 40
  k <- 10
  x <- rnorm(n * k)
  u1 <- rep(rnorm(k, sd = 0.2), each = n)
  u2 <- rep(rnorm(k, sd = 0.1), each = n)
  mu <- exp(4 +x + u1 + u2 * x)
  phi <- 50
  y <- rnbinom(n * k, mu = mu, size = phi)
  hist(y)
  d <- data.frame(year = 1:n, person = rep(1:k, each = n), y = y, x = x)

  fit_centered <- dynamite(
    obs(y ~ x + random(~ 1 + x), family = "negbin") +
      random_spec(noncentered = FALSE, correlated = TRUE),
    data = d, time = "year", group = "person",
    chains = 2, iter = 4000, refresh = 0, seed = 1
  )
  fit_noncentered <- dynamite(
    obs(y ~ x + random(~ 1 + x), family = "negbin") +
      random_spec(noncentered = TRUE, correlated = TRUE),
    data = d, time = "year", group = "person",
    chains = 2, iter = 4000, refresh = 0, seed = 1
  )

  expect_equal(
    summary(
      fit_centered,
      types = c("alpha", "beta", "corr_nu", "sigma_nu", "nu")
    )$mean,
    summary(
      fit_noncentered,
      types = c("alpha", "beta", "corr_nu", "sigma_nu", "nu")
    )$mean,
    tolerance = 0.1
  )
  expect_equal(
    summary(fit_centered, parameter = "phi_y")$mean,
    summary(fit_noncentered, parameter = "phi_y")$mean,
    tolerance = 0.2
  )

  fit_centered_nocorr <- dynamite(
    obs(y ~ x + random(~ 1 + x), family = "negbin") +
      random_spec(noncentered = FALSE, correlated = FALSE),
    data = d, time = "year", group = "person",
    chains = 2, iter = 4000, refresh = 0, seed = 1
  )
  fit_noncentered_nocorr <- dynamite(
    obs(y ~ x + random(~ 1 + x), family = "negbin") +
      random_spec(noncentered = TRUE, correlated = FALSE),
    data = d, time = "year", group = "person",
    chains = 2, iter = 4000, refresh = 0, seed = 1
  )

  expect_equal(
    summary(
      fit_centered_nocorr,
      types = c("alpha", "beta", "sigma_nu", "nu")
    )$mean,
    summary(
      fit_noncentered_nocorr,
      types = c("alpha", "beta","sigma_nu", "nu")
    )$mean,
    tolerance = 0.1
  )
  expect_equal(
    summary(
      fit_noncentered,
      types = c("alpha", "beta", "sigma_nu", "nu")
    )$mean,
    summary(
      fit_noncentered_nocorr,
      types = c("alpha", "beta","sigma_nu", "nu")
    )$mean,
    tolerance = 0.1
  )
})

test_that("random effects for categorical distribution work", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  f <- obs(x ~  lag(x) + lag(y) + random(~1 + z), family = "categorical") +
    obs(y ~ -1 + lag(x) + lag(y) + random(~z), family = "categorical")
  fit <- try(
    suppressWarnings(
      dynamite(
        dformula = f,
        data = categorical_example,
        time = "time",
        group = "id",
        chains = 1,
        iter = 2000,
        refresh = 0
      )
    )
  )
  expect_false(inherits(fit, "try-error"))
  expect_error(
    pred <- predict(fit, n_draws = 5),
    NA
  )
})

test_that("random effects for multinomial distribution work", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  d <- as.data.frame(t(rmultinom(100, 150, prob = c(0.15, 0.55, 0.3))))
  names(d) <- c("y1", "y2", "y3")
  d$t <- rep(1:20, each = 5)
  d$id <- rep(1:5, 20)
  d$x <- rnorm(100)
  d$n <- d$y1 + d$y2 + d$y3

  fit <- try(
    suppressWarnings(
      dynamite(
        obs(c(y1, y2, y3) ~ 1 + random(~ -1 + x) + trials(n), "multinomial"),
        data = d,
        time = "t",
        group = "id",
        chains = 1,
        iter = 2000,
        refresh = 0
      )
    ),
    silent = TRUE
  )
  expect_false(inherits(fit, "try-error"))
  expect_error(
    pred <- predict(fit, n_draws = 5),
    NA
  )
})

test_that("random effects for multivariate gaussian distribution work", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  N <- 5
  T_ <- 20
  S <- crossprod(matrix(rnorm(4), 2, 2))
  L <- t(chol(S))
  y1 <- matrix(0, N, T_)
  y2 <- matrix(0, N, T_)
  x <- matrix(0, N, T_)
  for (t in 2:T_) {
    for (i in 1:N){
      mu <- c(0.7 * y1[i, t-1], 0.4 * y2[i, t-1] - 0.2 * y1[i, t-1])
      y <- mu + L %*% rnorm(2L)
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
  fit <- try(
    suppressWarnings(
      dynamite(
        dformula =
          obs(c(y1, y2) ~ -1 + random(~ 1) + varying(~ x) | x, "mvgaussian") +
          splines(df = 10) ,
        data = d,
        time = "t",
        group = "id",
        chains = 1,
        iter = 2000,
        refresh = 0
      )
    )
  )
  expect_false(inherits(fit, "try-error"))
  expect_error(
    pred <- predict(fit, ndraws = 5),
    NA
  )
})

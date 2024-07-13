run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

data.table::setDTthreads(1) # For CRAN

# Capture both message and output types
capture_all_output <- function(x) {
  utils::capture.output(
    utils::capture.output(x, type = "message"),
    type = "output"
  )
}

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
      mu <- c(0.7 * y1[i, t - 1], 0.4 * y2[i, t - 1] - 0.2 * y1[i, t - 1])
      y <- mu + L %*% rnorm(2)
      y1[i, t] <- y[1L]
      y2[i, t] <- y[2L]
      x[i, t] <- rnorm(1, c(0.5 * y1[i, t - 1]), 0.5)
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

  expect_error(
    sumr <- summary(fit, types = "corr"),
    NA
  )
  expect_equal(
    sumr$mean,
    cov2cor(S)[2, 1],
    tolerance = 0.1
  )
  expect_error(
    sumr <- summary(fit, types = "sigma"),
    NA
  )
  expect_equal(
    sumr$mean,
    c(0.5, sqrt(diag(S))),
    tolerance = 0.1
  )
  expect_error(
    sumr <- summary(fit, types = "beta"),
    NA
  )
  expect_equal(
    sumr$mean,
    c(0.5, 0.7, -0.2, 0.4),
    tolerance = 0.1
  )

  expect_error(
    pred <- predict(fit, n_draws = 5),
    NA
  )
  expect_true(
    all(!is.na(pred[, c("y1_new", "y2_new", "x_new")])),
    NA
  )

  expect_error(
    pred <- predict(fit, n_draws = 5, type = "mean"),
    NA
  )
  pred <- pred[pred$t > 1, ]
  expect_true(
    all(!is.na(pred[, c("y1_mean", "y2_mean", "x_mean")])),
    NA
  )

  expect_error(
    pred <- predict(fit, n_draws = 5, type = "link"),
    NA
  )
  pred <- pred[pred$t > 1, ]
  expect_true(
    all(!is.na(pred[, c("y1_link", "y2_link", "x_link")])),
    NA
  )

  expect_error(
    pred <- fitted(fit, n_draws = 5),
    NA
  )
  pred <- pred[pred$t > 1, ]
  expect_true(
    all(!is.na(pred[, c("y1_fitted", "y2_fitted", "x_fitted")])),
    NA
  )
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
  expect_error(
    predict(fit, type = "mean"),
    NA
  )
  expect_error(
    predict(fit, type = "link"),
    NA
  )
  expect_error(
    fitted(fit),
    NA
  )
})


test_that("time-invariant cumulative fit and predict work", {
  skip_if_not(run_extended_tests)

  set.seed(0)

  n <- 100
  t <- 30
  x <- matrix(0, n, t)
  y <- matrix(0, n, t)
  p <- matrix(0, n, 4)
  alpha <- c(-1, 0, 1)

  for (i in seq_len(t)) {
    x[, i] <- rnorm(n)
    eta <- 0.6 * x[, i]
    p[, 1] <- 1 - plogis(eta - alpha[1])
    p[, 2] <- plogis(eta - alpha[1]) - plogis(eta - alpha[2])
    p[, 3] <- plogis(eta - alpha[2]) - plogis(eta - alpha[3])
    p[, 4] <- plogis(eta - alpha[3])
    y[, i] <- apply(p, 1, sample, x = 1:4, size = 1, replace = FALSE)
  }

  d <- data.frame(
    y = factor(c(y), levels = 1:4), x = c(x),
    time = rep(seq_len(t), each = n),
    id = rep(seq_len(n), t)
  )

  expect_error(
    fit <- dynamite(
      dformula =
        obs(y ~ x, family = "cumulative", link = "logit"),
      data = d,
      time = "time",
      group = "id"
    ),
    NA
  )
  expect_error(
    plot(fit),
    NA
  )
  expect_error(
    as.data.table(fit, types = "cutpoint"),
    NA
  )
  expect_error(
    pred <- predict(fit, type = "response", n_draws = 10),
    NA
  )
  expect_true(
    all(pred$y_new %in% 1:4)
  )
  expect_error(
    pred <- predict(fit, type = "mean", n_draws = 10),
    NA
  )
  expect_true(
    all(!is.na(pred[, paste0("y_mean_", 1:4)]))
  )
  expect_error(
    pred <- predict(fit, type = "link", n_draws = 10),
    NA
  )
  expect_true(
    all(!is.na(pred$y_link))
  )
  expect_error(
    pred <- fitted(fit, n_draws = 10),
    NA
  )
  expect_true(
    all(!is.na(pred[, paste0("y_fitted_", 1:4)]))
  )
})

test_that("time-varying cutpoints for cumulative works", {
  skip_if_not(run_extended_tests)

  set.seed(34)

  n <- 50
  t <- 30
  x <- matrix(0, n, t)
  y <- matrix(0, n, t)
  alpha_spline <- t(replicate(3, cumsum(rnorm(t, sd = 0.5))))
  p <- array(0, c(n, 4, t))
  cutpoints <- matrix(0, 3, t)

  for (i in seq_len(t)) {
    tmp <- exp(c(0, alpha_spline[, i]))
    for (j in 1:3) {
      cutpoints[j, i] <- log(sum(tmp[1:j]) / sum(tmp[(j + 1):4]))
    }
    x[, i] <- rnorm(n)
    eta <- 0.6 * x[, i]
    p[, 1, i] <- 1 - plogis(eta - cutpoints[1, i])
    p[, 2, i] <- plogis(eta - cutpoints[1, i]) - plogis(eta - cutpoints[2, i])
    p[, 3, i] <- plogis(eta - cutpoints[2, i]) - plogis(eta - cutpoints[3, i])
    p[, 4, i] <- plogis(eta - cutpoints[3, i])
    y[, i] <- apply(p[, , i], 1, sample, x = 1:4, size = 1, replace = TRUE)
  }
  d <- data.frame(
    y = factor(y, levels = 1:4),
    x = c(x),
    time = rep(seq_len(t), each = n),
    id = rep(seq_len(n), t)
  )

  expect_error(
    fit <- dynamite(
      dformula =
        obs(y ~ -1 + x + varying(~1), family = "cumulative", link = "logit") +
        splines(10),
      data = d,
      time = "time",
      group = "id"
    ),
    NA
  )
  expect_error(
    plot(fit),
    NA
  )
  expect_error(
    as.data.table(fit, types = "cutpoint"),
    NA
  )
  expect_error(
    pred <- predict(fit, type = "response", n_draws = 10),
    NA
  )
  expect_true(
    all(pred$y_new %in% 1:4)
  )
  expect_error(
    pred <- predict(fit, type = "mean", n_draws = 10),
    NA
  )
  expect_true(
    all(!is.na(pred[, paste0("y_mean_", 1:4)]))
  )
  expect_error(
    pred <- predict(fit, type = "link", n_draws = 10),
    NA
  )
  expect_true(
    all(!is.na(pred$y_link))
  )
  expect_error(
    pred <- fitted(fit, n_draws = 10),
    NA
  )
  expect_true(
    all(!is.na(pred[, paste0("y_fitted_", 1:4)]))
  )
})

test_that("non-glm categorical fit works", {
  skip_if_not(run_extended_tests)
  skip_on_os("mac") # Seems to segfault on MacOS
  expect_error(
    mockthat::with_mock(
      stan_supports_categorical_logit_glm = function(...) FALSE,
      dynamite(
        dformula = obs(x ~ z + lag(x) + lag(y), family = "categorical") +
          obs(y ~ z + lag(x) + lag(y), family = "categorical"),
        data = categorical_example,
        time = "time",
        group = "id",
        iter = 2000,
        chains = 1,
        refresh = 0,
        verbose = FALSE,
        threads_per_chain = 2
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

test_that("predict with random variable trials works", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  N <- 100
  T_ <- 50
  y <- matrix(0, N, T_)
  n <- matrix(0, N, T_)
  for (t in 2:T_) {
    for (i in 1:N) {
      n[i, t] <- rpois(1, 5)
      lp <- -3.0 + 1.25 * y[i, t - 1]
      y[i, t] <- rbinom(1, 1 + n[i, t], plogis(lp))
    }
  }
  d <- data.frame(
    n = c(n),
    y = c(y),
    t = rep(1:T_, each = N),
    id = 1:N
  )
  f <- obs(n ~ 1, family = "poisson") +
    aux(integer(nn) ~ n + 1) +
    obs(y ~ lag(y) + trials(nn), family = "binomial")
  fit <- dynamite(
    dformula = f,
    data = d,
    time = "t",
    group = "id",
    chains = 1,
    iter = 2000,
    refresh = 0
  )
  expect_error(sumr <- summary(fit), NA)
  expect_equal(sumr$mean, c(log(5), -3.0, 1.25), tolerance = 0.1)
  expect_error(predict(fit, n_draws = 5), NA)
})

test_that("shrinkage for splines is functional", {
  skip("Shrinkage feature removed at least for now.")

  set.seed(1)
  expect_error(
    fit <- dynamite(
      dformula =
        obs(
          y ~ -1 + z + varying(~ x + lag(y)) + random(~1),
          family = "gaussian"
        ) +
        random_spec() +
        splines(df = 20, shrinkage = TRUE),
      data = gaussian_example,
      time = "time",
      group = "id",
      iter = 2000,
      warmup = 1000,
      chains = 1,
      refresh = 0
    ),
    NA
  )
  expect_error(
    summary(fit, types = "xi"),
    NA
  )
})

test_that("update without recompile works", {
  skip_if_not(run_extended_tests)

  set.seed(0)
  gaussian_fit <- dynamite(
    dformula =
      obs(
        y ~ -1 + z + varying(~ x + lag(y)) + random(~1),
        family = "gaussian"
      ) +
      random_spec() +
      splines(df = 20),
    data = gaussian_example,
    time = "time",
    group = "id",
    iter = 2000,
    warmup = 1000,
    thin = 1,
    chains = 2,
    cores = 2,
    refresh = 0,
    save_warmup = FALSE,
    pars = c(
      "omega_alpha_1_y", "omega_raw_alpha_y", "nu_raw", "nu", "L",
      "sigma_nu", "a_y"
    ),
    include = FALSE
  )
  expect_error(
    fit <- update(
      gaussian_fit,
      data = gaussian_example,
      warmup = 1000,
      iter = 2000,
      thin = 1
    ),
    NA
  )
  # Internal update_ function
  expect_error(
    lfo(gaussian_fit, L = 20, chains = 4, verbose_stan = FALSE),
    NA
  )
})

test_that("custom stan model works", {
  skip_if_not(run_extended_tests)

  # The same as gaussian_example_fit, but we need to refit because
  # the results may not be reproducible across different platforms
  set.seed(1)
  initial_fit <- dynamite(
    dformula = obs(y ~ -1 + z + varying(~ x + lag(y)) +
                     random(~1), family = "gaussian") +
      random_spec() +
      splines(df = 20),
    data = gaussian_example,
    time = "time",
    group = "id",
    iter = 2000,
    warmup = 1000,
    thin = 10,
    chains = 2,
    cores = 2,
    refresh = 0,
    save_warmup = FALSE,
    pars = c(
      "omega_alpha_1_y", "omega_raw_alpha_y", "nu_raw", "nu", "L",
      "sigma_nu", "a_y"
    ),
    include = FALSE
  )
  code <- get_code(initial_fit)
  set.seed(1)
  expect_error(
    custom_fit <- dynamite(
      dformula = obs(y ~ -1 + z + varying(~ x + lag(y)) +
                       random(~1), family = "gaussian") +
        random_spec() +
        splines(df = 20),
      data = gaussian_example,
      time = "time",
      group = "id",
      iter = 2000,
      warmup = 1000,
      thin = 10,
      chains = 2,
      cores = 2,
      refresh = 0,
      save_warmup = FALSE,
      pars = c(
        "omega_alpha_1_y", "omega_raw_alpha_y", "nu_raw", "nu", "L",
        "sigma_nu", "a_y"
      ),
      include = FALSE,
      custom_stan_model = code
    ),
    NA
  )
  expect_equal(
    rstan::extract(initial_fit$stanfit, permuted = FALSE),
    rstan::extract(custom_fit$stanfit, permuted = FALSE)
  )
})

test_that("dynamice works", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  n <- 50
  p <- 1000
  y <- replicate(p, stats::arima.sim(list(ar = 0.7), n, sd = 0.1))
  d <- data.frame(y = c(y), time = 1:n, id = rep(seq_len(p), each = n))
  dmiss <- d
  dmiss$y[sample(seq_len(nrow(d)), size = 0.2 * nrow(d))] <- NA

  # Long format imputation
  expect_error(
    fit_long <- dynamice(
      obs(y ~ lag(y), "gaussian"),
      time = "time",
      group = "id",
      data = dmiss,
      chains = 1,
      refresh = 0,
      backend = "rstan",
      impute_format = "long",
      keep_imputed = FALSE,
      mice_args = list(m = 3, print = FALSE)
    ),
    NA
  )

  # Wide format imputation
  expect_error(
    fit_wide <- dynamice(
      obs(y ~ lag(y), "gaussian"),
      time = "time",
      group = "id", data = dmiss,
      chains = 1,
      refresh = 0,
      backend = "rstan",
      impute_format = "wide",
      keep_imputed = FALSE,
      mice_args = list(m = 3, print = FALSE)
    ),
    NA
  )

  # single group
  set.seed(1)
  n <- 100
  y <- stats::arima.sim(list(ar = 0.7), n, sd = 0.1)
  d <- data.frame(y = c(y), time = 1:n)
  dmiss_single <- d
  dmiss_single$y[sample(seq_len(nrow(d)), size = 0.2 * nrow(d))] <- NA
  expect_error(
    fit_single <- dynamice(
      obs(y ~ lag(y), "gaussian"),
      time = "time",
      data = dmiss_single,
      chains = 1,
      refresh = 0,
      backend = "rstan",
      impute_format = "long",
      keep_imputed = FALSE,
      mice_args = list(m = 3, print = FALSE)
    ),
    NA
  )
})

test_that("information on >2 chains is summarized in print", {
  skip_if_not(run_extended_tests)

  set.seed(1)
  fit <- dynamite(
    dformula =
      obs(y ~ -1 + z +
            varying(~ x + lag(y)) + random(~1), family = "gaussian") +
      random_spec() +
      splines(df = 20),
    data = gaussian_example,
    time = "time",
    group = "id",
    iter = 2000,
    warmup = 1000,
    thin = 10,
    chains = 4,
    refresh = 0,
    save_warmup = FALSE,
    pars = c(
      "omega_alpha_1_y", "omega_raw_alpha_y", "nu_raw", "nu", "L",
      "sigma_nu", "a_y"
    ),
    include = FALSE
  )
  out <- capture_all_output(print(fit))
  expect_true(
    any(
      grepl(
        "Elapsed time (seconds) for fastest and slowest chains",
        out,
        fixed = TRUE
      )
    )
  )
})


test_that("latent factor models are identifiable", {
  skip_if_not(run_extended_tests)

  generate_data <- function(N, T_, D, alpha, mean_lambda, sd_lambda, sd_alpha) {
    y <- matrix(0, N, T_)
    x <- matrix(rnorm(N * T_), N, T_)
    psi <-
      c(splines::bs(seq_len(T_), df = D, intercept = TRUE) %*% cumsum(rnorm(D)))
    lambda <- rnorm(N, 0, sd_lambda)
    lambda <- mean_lambda + lambda - mean(lambda)
    a <- rnorm(N, alpha, sd_alpha)
    for (t in 1:T_) {
      y[, t] <- rnorm(N, a + lambda * psi[t] + x[, t])
    }
    list(data = data.frame(
      y = c(y), x = c(x),
      time = rep(seq_len(T_), each = N),
      id = rep(seq_len(N), times = T_)
    ),
    psi = psi, lambda = lambda,
    mean_lambda_psi = mean(lambda) * psi,
    alpha = alpha, beta = 1, kappa = sd_lambda / (1 + sd_lambda),
    sigma_lambda = sd_lambda, sigma_alpha = sd_alpha, sigma_y = 1, tau_psi = 1,
    zeta = 1 + sd_lambda)
  }

  set.seed(1)
  sim <- generate_data(
    N = 100,
    T_ = 50,
    D = 20,
    alpha = 1,
    mean_lambda = 1,
    sd_lambda = 1,
    sd_alpha = 1
  )
  dformula <- obs(y ~ x + random(~ 1), family = "gaussian") +
    lfactor() + splines(20)
  priors <- get_priors(dformula, data = sim$data, time = "time", group = "id")
  priors$prior[priors$parameter != "kappa_y"] <- "normal(0, 5)"
  fit1 <- dynamite(
    dformula, priors = priors,
    data = sim$data, time = "time", group = "id",
    backend = "cmdstanr", stanc_options = list("O1"),
    iter_sampling = 5000, iter_warmup = 5000,
    parallel_chains = 4, refresh = 0, seed = 1
  )
  sumr1 <- as_draws(fit1) |>
    posterior::summarise_draws()
  expect_true(all(sumr1$rhat < 1.1))
  expect_true(all(sumr1$ess_bulk > 500))
  expect_true(all(sumr1$ess_tail > 500))
  expect_equal(summary(fit1, type = "psi")$mean, sim$mean_lambda_psi,
               tolerance = 0.5)

  set.seed(2)
  sim <- generate_data(
    100, 50, 20,
    alpha = 0, mean_lambda = 1, sd_lambda = 0.1, sd_alpha = 2)
  dformula <- obs(y ~ -1 + x + random(~ 1), family = "gaussian") +
    lfactor() + splines(20)
  priors <- get_priors(dformula, data = sim$data, time = "time", group = "id")
  priors$prior[priors$parameter != "kappa_y"] <- "normal(0, 5)"
  fit2 <- dynamite(
    dformula, priors = priors,
    data = sim$data, time = "time", group = "id",
    backend = "cmdstanr", stanc_options = list("O1"),
    iter_sampling = 5000, iter_warmup = 5000,
    parallel_chains = 4, refresh = 0, seed = 1
  )
  sumr2 <- as_draws(fit2) |>
    posterior::summarise_draws()
  expect_true(all(sumr2$rhat < 1.1))
  expect_true(all(sumr2$ess_bulk > 500))
  expect_true(all(sumr2$ess_tail > 500))
  expect_equal(summary(fit2, type = "psi")$mean, sim$mean_lambda_psi,
               tolerance = 0.5)

  set.seed(3)
  sim <- generate_data(
    100, 50, 20,
    alpha = -1, mean_lambda = 0.5, sd_lambda = 2, sd_alpha = 0)
  dformula <- obs(y ~ x, family = "gaussian") +
    lfactor() + splines(20)
  priors <- get_priors(dformula, data = sim$data, time = "time", group = "id")
  priors$prior[priors$parameter != "kappa_y"] <- "normal(0, 5)"
  # need some informativeness on prior of kappa (or zeta?)
  priors$prior[priors$parameter == "kappa_y"] <- "beta(10, 10)"
  fit3 <- dynamite(
    dformula, priors = priors,
    data = sim$data, time = "time", group = "id",
    backend = "cmdstanr", stanc_options = list("O1"),
    iter_sampling = 5000, iter_warmup = 5000,
    parallel_chains = 4, refresh = 0, seed = 1
  )
  sumr3 <- as_draws(fit3) |>
    posterior::summarise_draws()

  expect_true(all(sumr3$rhat < 1.1))
  expect_true(all(sumr3$ess_bulk > 500))
  expect_true(all(sumr3$ess_tail > 500))
  expect_equal(summary(fit3, type = "psi")$mean, sim$mean_lambda_psi,
               tolerance = 0.5)

  set.seed(4)
  sim <- generate_data(
    100, 50, 20,
    alpha = 0, mean_lambda = -1, sd_lambda = 0.5, sd_alpha = 0.5)
  dformula <- obs(y ~ x + random(~ 1), family = "gaussian") +
    lfactor() + splines(20)
  priors <- get_priors(dformula, data = sim$data, time = "time", group = "id")
  priors$prior[priors$parameter != "kappa_y"] <- "normal(0, 5)"
  fit4 <- dynamite(
    dformula, priors = priors,
    data = sim$data, time = "time", group = "id",
    backend = "cmdstanr", stanc_options = list("O1"),
    iter_sampling = 5000, iter_warmup = 5000,
    parallel_chains = 4, refresh = 0, seed = 1
  )
  sumr4 <- as_draws(fit4) |>
    posterior::summarise_draws()

  expect_true(all(sumr4$rhat < 1.1))
  expect_true(all(sumr4$ess_bulk > 500))
  expect_true(all(sumr4$ess_tail > 500))
  expect_equal(summary(fit4, type = "psi")$mean, sim$mean_lambda_psi,
               tolerance = 0.5)

  set.seed(5)
  sim <- generate_data(
    100, 50, 20,
    alpha = 1, mean_lambda = 0, sd_lambda = 0.5, sd_alpha = 1)
  dformula <- obs(y ~ x + random(~ 1), family = "gaussian") +
    lfactor(nonzero_lambda = FALSE) + splines(20)
  priors <- get_priors(dformula, data = sim$data, time = "time", group = "id")
  priors$prior[] <- "normal(0, 5)"
  fit5 <- dynamite(
    dformula, priors = priors,
    data = sim$data, time = "time", group = "id",
    backend = "cmdstanr", stanc_options = list("O1"),
    iter_sampling = 5000, iter_warmup = 5000,
    parallel_chains = 4, refresh = 0, seed = 1
  )
  sumr5 <- as_draws(fit5) |>
    posterior::summarise_draws()

  expect_true(all(sumr5$rhat < 1.1))
  expect_true(all(sumr5$ess_bulk > 500))
  expect_true(all(sumr5$ess_tail > 500))
  expect_equal(summary(fit5, type = "psi")$mean, -sim$psi,
               tolerance = 0.5)

  set.seed(6)
  sim <- generate_data(
    100, 50, 20,
    alpha = 1, mean_lambda = 0, sd_lambda = 0.5, sd_alpha = 0)
  dformula <- obs(y ~ x, family = "gaussian") +
    lfactor(nonzero_lambda = FALSE) + splines(20)
  priors <- get_priors(dformula, data = sim$data, time = "time", group = "id")
  priors$prior[] <- "normal(0, 5)"
  fit6 <- dynamite(
    dformula, priors = priors,
    data = sim$data, time = "time", group = "id",
    backend = "cmdstanr", stanc_options = list("O1"),
    iter_sampling = 5000, iter_warmup = 5000,
    parallel_chains = 4, refresh = 0, seed = 1
  )
  sumr6 <- as_draws(fit6) |>
    posterior::summarise_draws()

  expect_true(all(sumr6$rhat < 1.1))
  expect_true(all(sumr6$ess_bulk > 500))
  expect_true(all(sumr6$ess_tail > 500))
  expect_equal(summary(fit6, type = "psi")$mean, -sim$psi,
               tolerance = 0.5)

  # Test bivariate case with nonzero_lambda
  set.seed(123)
  N <- 50
  T_ <- 50
  D <- 10
  x <- y <- matrix(0, N, T_)
  psi <- matrix(NA, 2, T_)
  lambda_y <- rnorm(N)
  lambda_y <- lambda_y - mean(lambda_y)
  lambda_x <- rnorm(N)
  lambda_x <- lambda_x - mean(lambda_x)
  L <- t(chol(matrix(c(1, 0.7, 0.7, 1), 2, 2)))
  B <- t(splines::bs(seq_len(T_), df = D, intercept = TRUE))
  omega <- matrix(NA, 2, D)
  omega[, 1] <- c(5, 2) + L %*% rnorm(2)
  for (i in 2:D) {
    omega[, i] <- omega[, i - 1] + L %*% rnorm(2)
  }
  psi[1, ] <- omega[1, ] %*% B
  psi[2, ] <- omega[2, ] %*% B
  for (t in 1:T_) {
    y[, t] <- rnorm(N, lambda_y * psi[1, t])
    x[, t] <- rnorm(N, lambda_x * psi[2, t])
  }
  d <- data.frame(
    y = c(y), x = c(x),
    time = rep(seq_len(T_), each = N),
    id = rep(seq_len(N), times = T_)
  )
  dformula <- obs(y ~ 1, family = "gaussian") +
    obs(x ~ 1, family = "gaussian") +
    lfactor(nonzero_lambda = FALSE, flip_sign = TRUE) + splines(10)
  fit <- dynamite(
    dformula,
    data = d, time = "time", group = "id",
    backend = "cmdstanr", stanc_options = list("O1"),
    iter_sampling = 5000, iter_warmup = 5000,
    parallel_chains = 4, refresh = 0, seed = 1
  )
  sumr <- as_draws(fit) |>
    posterior::summarise_draws()

  expect_true(all(sumr$rhat < 1.1))
  expect_true(all(sumr$ess_bulk > 500))
  expect_true(all(sumr$ess_tail > 500))
})

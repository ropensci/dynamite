#' @srrstats {G5.10} Extended tests can be switched on via setting the
#'   environment variable DYNAMITE_EXTENDED_TESTS to "true".
#' @srrstats {G5.5, G5.6b} Seeds are used appropriately in the tests.
#' @srrstats {G5.4, G5.4a, G5.4b, G5.4c, G5.6, G5.6a, BS7.0, BS7.1, BS7.2}
#'   Simple linear regression and GLM models are tested so that they match with
#'   lm and glm function outputs (within a tolerance due to MCMC, use of
#'   default priors, and discrepancy between ML estimate vs posterior mean).
#'   Further recovery and correctness tests are also implemented.
#' @srrstats {G5.7} Tested that the parameters of the true data generating
#'   process are recovered when increasing the data size.
run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

test_that("parameters for the linear regression are recovered as with lm", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  y <- 2 - 1 * x + rnorm(n, sd = 0.1)
  d <- data.frame(time = 1:n, y = y, x = x)

  fit_lm <- lm(y ~ x, data = d)
  priors <- get_priors(obs(y ~ x, family = "gaussian"),
    data = d, time = "time"
  )
  priors$prior <- c("normal(0, 5)", "std_normal()", "exponential(1)")
  fit_dynamite <- dynamite(obs(y ~ x, family = "gaussian"),
    data = d, time = "time", priors = priors,
    chains = 1, iter = 2000, refresh = 0
  )
  expect_equal(coef(fit_dynamite)$mean, coef(fit_lm),
    tolerance = 0.01,
    ignore_attr = TRUE
  )
})

test_that("parameters for the poisson glm are recovered as with glm", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  y <- rpois(n, exp(2 - 1 * x))
  d <- data.frame(time = 1:n, y = y, x = x)

  fit_glm <- glm(y ~ x, data = d, family = poisson)
  fit_dynamite <- dynamite(obs(y ~ x, family = "poisson"),
    data = d, time = "time", chains = 1, iter = 2000, refresh = 0
  )
  expect_equal(coef(fit_dynamite)$mean, coef(fit_glm),
    tolerance = 0.01,
    ignore_attr = TRUE
  )
})

test_that("parameters for the binomial glm are recovered as with glm", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  n <- 100
  u <- sample(1:10, n, TRUE)
  x <- rnorm(n)
  y <- rbinom(n, u, plogis(1 - x))
  d <- data.frame(time = 1:n, y = y, x = x, u = u)

  fit_glm <- glm(cbind(y, u - y) ~ x, data = d, family = binomial)
  fit_dynamite <- dynamite(obs(y ~ x + trials(u), family = "binomial"),
    data = d, time = "time", chains = 1, iter = 2000, refresh = 0
  )
  expect_equal(coef(fit_dynamite)$mean, coef(fit_glm),
    tolerance = 0.01,
    ignore_attr = TRUE
  )
})

test_that("parameters for the gamma glm are recovered as with glm", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  y <- rgamma(n, 2, 2 / exp(1 - 2 * x))
  d <- data.frame(time = 1:n, y = y, x = x)

  fit_glm <- glm(y ~ x, data = d, family = Gamma(link = "log"))
  fit_dynamite <- dynamite(obs(y ~ x, family = "gamma"),
    data = d, time = "time", chains = 1, iter = 2000, refresh = 0
  )
  expect_equal(coef(fit_dynamite)$mean[1:2], coef(fit_glm),
    tolerance = 0.01,
    ignore_attr = TRUE
  )
})

test_that("parameters for poisson mixed model are recovered", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  n <- 40
  k <- 10
  x <- rnorm(n * k)
  u1 <- rep(rnorm(k, sd = 0.2), each = n)
  u2 <- rep(rnorm(k, sd = 0.1), each = n)
  y <- rpois(n * k, exp(2 - x + u1 + u2 * x))
  d <- data.frame(year = 1:n, person = rep(1:k, each = n), y = y, x = x)

  p <- data.frame(
    parameter = c(
      "sigma_nu_y_alpha", "sigma_nu_y_x", "alpha_y", "beta_y_x",
      "L_nu"
    ),
    response = c(rep("y", 4), ""),
    prior = c(
      "std_normal()", "std_normal()",
      "student_t(3, 2, 2)", "normal(0, 10)", "lkj_corr_cholesky(1)"
    ),
    type = c("sigma_nu", "sigma_nu", "alpha", "beta", "L"),
    category = ""
  )
  fit_dynamite <- dynamite(
    obs(y ~ x + random(~ 1 + x), family = "poisson"),
    data = d, time = "year", group = "person", priors = p,
    init = 0, chains = 2, cores = 2, iter = 2000, refresh = 0, seed = 1
  )
  # "ground truth" obtained from one long dynamite run
  expect_equal(coef(fit_dynamite)$mean, c(2, -0.99),
    tolerance = 0.1
  )
  expect_equal(coef(fit_dynamite, type = "nu")$mean,
    c(
      0.17, 0.42, -0.09, -0.13, -0.07, -0.12, -0.2, -0.12, 0.28,
      -0.1, -0.03, 0, 0.1, -0.11, -0.03, 0.02, 0.04, -0.02, -0.14, 0.16
    ),
    tolerance = 0.1
  )
})

test_that("parameters for an AR(1) model are recovered as with arima", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  fit <- dynamite(obs(LakeHuron ~ 1, "gaussian") + lags(),
    data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
    time = "time",
    group = "id",
    chains = 1,
    iter = 2000,
    refresh = 0
  )
  fit_arima <- arima(LakeHuron, c(1, 0, 0))
  expect_equal(coef(fit)$mean[2], coef(fit_arima)[1L],
    tolerance = 0.01,
    ignore_attr = TRUE
  )
  expect_equal(
    coef(fit)$mean[1L],
    coef(fit_arima)[2L] * (1 - coef(fit_arima)[1L]),
    tolerance = 1, ignore_attr = TRUE
  )
})

test_that("LOO works for AR(1) model", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  fit <- dynamite(obs(LakeHuron ~ 1, "gaussian") + lags(),
    data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
    time = "time",
    group = "id",
    chains = 1,
    iter = 2000,
    refresh = 0
  )
  l <- loo(fit)
  expect_equal(
    l$estimates,
    structure(
      c(
        -107.877842970846, 2.86041434691809, 215.755685941693,
        7.36848739076899, 0.561813071004331, 14.736974781538
      ),
      dim = 3:2,
      dimnames = list(c("elpd_loo", "p_loo", "looic"), c("Estimate", "SE"))
    ),
    tolerance = 1
  )
  expect_error(plot(l), NA)
})

test_that("LOO works with separate channels", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  # Fit again so that recompile with update works with all platforms
  multichannel_fit <- dynamite(
    dformula = obs(g ~ lag(g) + lag(logp), family = "gaussian") +
      obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
      obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
      aux(numeric(logp) ~ log(p + 1)),
    data = multichannel_example,
    time = "time",
    group = "id",
    verbose = FALSE,
    chains = 1,
    cores = 1,
    iter = 2000,
    warmup = 1000,
    init = 0,
    refresh = 0,
    thin = 5,
    save_warmup = FALSE
  )
  expect_error(
    l <- loo(update(multichannel_fit, thin = 1), separate_channels = TRUE),
    NA
  )
  expect_equal(
    l$g_loglik$estimates,
    structure(
      c(
        127.7731689, 3.9598420, -255.5463377,
        21.1943047,  0.2433661,   42.3886094
      ),
      dim = 3:2,
      dimnames = list(c("elpd_loo", "p_loo", "looic"), c("Estimate", "SE"))
    ),
    tolerance = 1
  )
  expect_equal(
    l$p_loglik$estimates,
    structure(
      c(
        -2128.5452197, 4.5260226, 4257.0904393,
        26.5452884,    0.3107372, 53.0905768
      ),
      dim = 3:2,
      dimnames = list(c("elpd_loo", "p_loo", "looic"), c("Estimate", "SE"))
    ),
    tolerance = 1
  )
  expect_equal(
    l$b_loglik$estimates,
    structure(
      c(
        -583.3724555, 6.8573891, 1166.7449111,
        12.1459613,   0.3097697, 24.2919227
      ),
      dim = 3:2,
      dimnames = list(c("elpd_loo", "p_loo", "looic"), c("Estimate", "SE"))
    ),
    tolerance = 1
  )
})

test_that("LFO works for AR(1) model", {
  # This also implicitly tests update method
  skip_if_not(run_extended_tests)
  set.seed(1)
  fit <- dynamite(obs(LakeHuron ~ 1, "gaussian") + lags(),
    data = data.frame(LakeHuron, time = seq_len(length(LakeHuron)), id = 1),
    time = "time",
    group = "id",
    chains = 1,
    iter = 2000,
    refresh = 0
  )
  l <- lfo(fit, L = 20)
  expect_equal(l$ELPD, -90.4188604974201, tolerance = 1)
  expect_equal(l$ELPD_SE, 7.58842574523583, tolerance = 1)
  expect_error(plot(l), NA)
  expect_error(print(l), NA)
})

test_that("parameters of a time-varying gaussian model are recovered", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  create_data <- function(N = 10L, T_ = 100L, D = 50L) {
    K_fixed <- 1L
    K_varying <- 2L
    tau <- c(0.2, 0.4)
    sigma <- 0.1
    beta <- 2.0
    Bs <-
      t(splines::bs(seq.int(1L, T_), df = D, degree = 3L, intercept = TRUE))
    D <- nrow(Bs)
    a <- array(0.0, c(K_varying, D))
    delta <- array(NA, c(T_, K_varying))
    for (k in seq_len(K_varying)) {
      a[k, ] <- cumsum(rnorm(D, 0, tau[k]))
      for (t in seq.int(1L, T_)) {
        delta[t, k] <- a[k, ] %*% Bs[, t]
      }
    }
    x <- matrix(rnorm(T_ * N), N, T_)
    z <- matrix(rbinom(T_ * N, 1.0, 0.7), N, T_)
    y <- matrix(NA, N, T_)
    y[, 1L] <- rnorm(N)
    for (t in seq.int(1L, T_)) {
      m <- beta * z[, t] + delta[t, 1L] +
        delta[t, 2L] * x[, t]
      y[, t] <- rnorm(N, m, sigma)
    }

    list(
      data = data.frame(
        y = c(y), x = c(x), z = c(z), id = seq_len(N),
        time = rep(seq_len(T_), each = N)
      ),
      true_values = c(delta = delta, tau = tau, beta = beta, sigma = sigma)
    )
  }
  d <- create_data()
  dformula <- obs(y ~ -1 + z + varying(~x), family = "gaussian") +
    splines(df = 50)
  # compile model only once
  code <- get_code(dformula,
    data = d$data,
    time = "time",
    group = "id"
  )
  model <- rstan::stan_model(model_code = code)

  # simulate multiple datasets
  n <- 10
  diffs <- matrix(NA, length(d$true_values), n)
  pars <- c("alpha_y", "delta_y", "tau_alpha_y", "tau_y", "beta_y", "sigma_y")
  for (i in seq_len(n)) {
    data <- get_data(dformula, group = "id", time = "time", data = d$data)
    diffs[, i] <- rstan::get_posterior_mean(
      rstan::sampling(model,
        data = data,
        refresh = 0, chains = 1, iter = 2000,
        pars = pars
      ),
      pars = pars
    ) - d$true_values
    d <- create_data()
  }
  # small MSE
  expect_lt(mean(diffs^2), 0.005)

  # test with a single large dataset
  d <- create_data(T_ = 500, N = 500, D = 100)
  data <- get_data(obs(y ~ -1 + z + varying(~x), family = "gaussian") +
      splines(df = 100), time = "time", group = "id", data = d$data)
  fit_long <- rstan::sampling(model,
    data = data,
    refresh = 0, chains = 1, iter = 2000,
    pars = pars
  )
  estimates <- c(rstan::get_posterior_mean(fit_long, pars = pars))
  expect_equal(c(estimates), d$true_values,
    ignore_attr = TRUE, tolerance = 0.1
  )
})

test_that("prior parameters are recovered with zero observations", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  d <- data.frame(y = rep(NA, 10), x = rnorm(10), id = 1, time = 1:10)
  p <- get_priors(obs(y ~ x, "gaussian"), d, time = "time", group = "id")
  p$prior[] <- c("normal(2, 0.1)", "normal(5, 0.5)", "exponential(10)")
  fit_prior <- dynamite(obs(y ~ x, "gaussian"),
    data = d,
    time = "time",
    group = "id",
    priors = p,
    iter = 55000,
    warmup = 5000,
    chains = 1,
    cores = 1,
    refresh = 0,
    save_warmup = FALSE
  )
  sumr <- summary(fit_prior) |>
    dplyr::select(parameter, mean, sd, q5, q95) |>
    as.data.frame()

  sigma_y <- sumr |>
    dplyr::filter(parameter == "alpha_y") |>
    dplyr::select(mean, sd, q5, q95)

  m <- 2 - d$x[1L] * 5
  s <- sqrt(0.1^2 + d$x[1L]^2 * 0.5^2)
  expect_equal(
    unlist(sumr[1, 2:5]),
    c(m, s, qnorm(c(0.05, 0.95), m, s)),
    tolerance = 0.1, ignore_attr = TRUE
  )
  expect_equal(
    unlist(sumr[2, 2:5]),
    c(5, 0.5, qnorm(c(0.05, 0.95), 5, 0.5)),
    tolerance = 0.1, ignore_attr = TRUE
  )
  expect_equal(
    unlist(sumr[3, 2:5]),
    c(0.1, 0.1, qexp(c(0.05, 0.95), 10)),
    tolerance = 0.1, ignore_attr = TRUE
  )
})

test_that("predict recovers correct estimates", {
  skip_if_not(run_extended_tests)
  set.seed(1)
  N <- 20
  T_ <- 30
  y <- matrix(0, N, T_)
  nu <- rnorm(N)
  y[, 1] <- rbinom(N, size = 1, prob = 0.5)
  for (t in 2:T_) y[, t] <- rbinom(N, 1, plogis(nu + y[, t-1]))

  ## check these if tests fail ##
  # model <- rstan::stan_model("testmodel.stan")
  # fit <- rstan::sampling(model, data = list(N = N, T = T_, y = y), chains = 1,
  #   iter = 2e4, warmup = 1000)
  # rstan_obs_results_id1_time4 <- rstan::summary(fit, "y_rep[1, 4]",
  #   use_cache = FALSE)$summary[, 1:3]
  # rstan_obs_results_avg4 <- setNames(c(rstan::summary(fit,
  #   c("mean_y[4]", "sd_y[4]"),
  #   use_cache = FALSE)$summary[, 1:3]),
  # c("mean_m", "mean_s", "se_m", "se_s", "sd_m", "sd_s"))
  #
  # rstan_prob_results_id1_time4 <- rstan::summary(fit, "y_m[1, 4]",
  #   use_cache = FALSE)$summary[, 1:3]
  # rstan_prob_results_avg4 <- setNames(c(rstan::summary(fit,
  #   c("mean_y_m[4]", "sd_y_m[4]"),
  #   use_cache = FALSE)$summary[, 1:3]),
  #   c("mean_m", "mean_s", "se_m", "se_s", "sd_m", "sd_s"))

  rstan_obs_results_id1_time4 <- c(mean = 0.6098, se_mean = 0.0035, sd = 0.4878)
  rstan_obs_results_avg4 <- c(mean_m = 0.7136, mean_s = 0.4508,
    se_m = 7e-04, se_s = 4e-04, sd_m = 0.0939, sd_s = 0.0511)
  rstan_prob_results_id1_time4 <- c(mean = 0.6062, se_mean = 9e-04, sd = 0.1409)
  rstan_prob_results_avg4 <- c(mean_m = 0.7138, mean_s = 0.213, se_m = 2e-04,
    se_s = 2e-04, sd_m = 0.0264, sd_s = 0.025)

  d <- data.frame(y = c(y), time = rep(1:T_, each = N), id = 1:N)
  p <- get_priors(obs(y ~ lag(y) + random(~1), "bernoulli"),
    data = d, time = "time", group = "id")
  p$prior[] <- "std_normal()"
  fitd <- dynamite(obs(y ~ lag(y) + random(~1), "bernoulli"),
    data = d, time = "time", group = "id", priors = p,
    chains = 1, iter = 2e4, warmup = 1000, refresh = 0)

  pred <- predict(fitd)
  y_new <- pred$y_new[pred$time == 4 & pred$id == 1]
  expect_equal(
    c(
      mean = mean(y_new),
      se_mean = sd(y_new) / sqrt(length(y_new)),
      sd = sd(y_new)
    ),
    rstan_obs_results_id1_time4,
    tolerance = 0.05
  )

  res <- pred |>
    dplyr::filter(time == 4) |>
    dplyr::group_by(.draw) |>
    dplyr::summarise(m = mean(y_new), s = sd(y_new)) |>
    dplyr::summarise(
      mean_m = mean(m), mean_s = mean(s),
      se_m = sd(m) / sqrt(dplyr::n()), se_s = sd(s) / sqrt(dplyr::n()),
      sd_m = sd(m), sd_s = sd(s)
      ) |> unlist()

  expect_equal(res,
    rstan_obs_results_avg4,
    tolerance = 0.01
  )

  pred_m <- predict(fitd, type = "mean")
  y_mean <- pred_m$y_mean[pred_m$time == 4 & pred_m$id == 1]
  expect_equal(
    c(
      mean = mean(y_mean),
      se_mean = sd(y_mean) / sqrt(length(y_mean)),
      sd = sd(y_mean)
    ),
    rstan_prob_results_id1_time4,
    tolerance = 0.01
  )

  res <- pred_m |>
    dplyr::filter(time == 4) |>
    dplyr::group_by(.draw) |>
    dplyr::summarise(m = mean(y_mean), s = sd(y_mean)) |>
    dplyr::summarise(
      mean_m = mean(m), mean_s = mean(s),
      se_m = sd(m) / sqrt(dplyr::n()), se_s = sd(s) / sqrt(dplyr::n()),
      sd_m = sd(m), sd_s = sd(s)
    ) |> unlist()

  expect_equal(res,
    rstan_prob_results_avg4,
    tolerance = 0.01
  )
})

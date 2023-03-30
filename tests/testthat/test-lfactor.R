#' @srrstats {G5.10} Extended tests can be switched on via setting the
#'   environment variable DYNAMITE_EXTENDED_TESTS to "true".

run_extended_tests <- identical(Sys.getenv("DYNAMITE_EXTENDED_TESTS"), "true")

set.seed(1)
T_ <- 20
N <- 50
x <- matrix(rnorm(T_ * N, 2, 0.5), N, T_)
D <- 10
B <- t(splines::bs(1:T_, df = D, intercept = TRUE))
a1 <- cumsum(rnorm(D, 0, 0.1))
a2 <- cumsum(rnorm(D, 0, 0.5))
psi1 <- numeric(T_)
psi2 <- numeric(T_)
for (t in 1:T_) {
  psi1[t] <- B[, t] %*% a1
  psi2[t] <- B[, t] %*% a2
}
lambda1 <- rnorm(N, 0.4, 1)
lambda2 <- rnorm(N, 0, 0.5)
y1 <- matrix(0, N, T_)
y2 <- matrix(0, N, T_)
for (t in 1:T_) {
  y1[, t] <- rpois(N, exp(2 + x[, t] + lambda1 * psi1[t]))
  y2[, t] <- rnorm(N, 1 + x[, t] + lambda2 * psi2[t], 0.2)
}
d <- data.frame(
  y1 = c(y1),
  y2 = c(y2),
  x = c(x),
  id = seq_len(N),
  time = rep(seq_len(T_), each = N)
)

test_that("nonidentifiable lfactor specification gives warning", {
  expect_error(
    dynamite(
      obs(y1 ~ -1 + x, family = "poisson") +
        obs(y2 ~ x, family = "gaussian") +
        lfactor(responses = c("y1", "y2"),
          nonzero_lambda = TRUE,
          correlated = TRUE,
          noncentered_psi = TRUE
        ) + splines(30),
      data = d, time = "time", group = "id", debug = list(no_compile = TRUE)),
    NA
  )
  expect_warning(
    dynamite(
      obs(y1 ~ x, family = "poisson") +
        obs(y2 ~ x + random(~1), family = "gaussian") +
        lfactor(
          responses = c("y1", "y2"),
          nonzero_lambda = TRUE,
          correlated = TRUE,
          noncentered_psi = TRUE
        ) +
        splines(30),
      data = d,
      time = "time",
      group = "id",
      debug = list(no_compile = TRUE)
    ),
    paste0("The common time-invariant intercept term of channel `y2` was ",
      "removed as channel predictors contain random intercept and latent ",
      "factor specified with `nonzero_lambda` as TRUE\\.")
  )
})

# WIP
# test_that("latent factors work", {
#   skip_if_not(run_extended_tests)
#
#   fit1 <- dynamite(
#     obs(y1 ~ x, family = "poisson") + obs(y2 ~ x, family = "gaussian") +
#       lfactor(nonzero_lambda = c(TRUE, FALSE),
#         noncentered_psi = TRUE) +
#       splines(10),
#     data = d, time = "time", group = "id",
#     chains = 1, refresh = 0, seed = 1
#   )
#
#   fit2 <- dynamite(
#     obs(y1 ~ x, family = "poisson") + obs(y2 ~ x, family = "gaussian") +
#       lfactor(
#         nonzero_lambda = c(TRUE, FALSE),
#         noncentered_psi = FALSE) +
#       splines(10),
#     data = d, time = "time", group = "id",
#     chains = 1, refresh = 0, seed = 1
#   )
#
#   fit3 <- dynamite(
#     obs(y1 ~ x, family = "poisson") + obs(y2 ~ x, family = "gaussian") +
#       lfactor(
#         nonzero_lambda = c(TRUE, FALSE),
#         correlated = FALSE,
#         noncentered_psi = FALSE) +
#       splines(10),
#     data = d, time = "time", group = "id",
#     chains = 1, refresh = 0, seed = 1
#   )
#   as_draws(fit3,types=c("alpha","beta","sigma_lambda","tau_psi")) |> posterior:::summarise_draws()
#   expect_equal(
#     summary(fit1, types = c("alpha", "beta", "sigma"))$mean,
#     c(1, 1.05, 1, 0.53),
#     tolerance = 0.1
#   )
#
#   expect_equal(
#     summary(fit1)$mean,
#     summary(fit2)$mean,
#     tolerance = 0.1
#   )
#   expect_equal(
#     summary(fit1, types = c("alpha", "beta", "sigma"))$mean,
#     summary(fit3, types = c("alpha", "beta", "sigma"))$mean,
#     tolerance = 0.1
#   )
# })

# Tests involving `latent_factor_example` and `latent_factor_example_fit` -----

set.seed(123)
N <- 40L
T_ <- 20L
D <- 10
B <- t(splines::bs(1:T_, df = D, intercept = TRUE))
a <- cumsum(rnorm(D))
psi <- numeric(T_)
lambda_i <- rnorm(N, 1, 0.2)
for (t in 1:T_) {
  psi[t] <- B[, t] %*% a
}
y <- matrix(0, N, T_)
for (t in 1:T_) {
  y[, t] <- rnorm(N, lambda_i * psi[t], 0.2)
}
latent_factor_example <- data.frame(
  y = c(y),
  id = seq_len(N),
  time = rep(seq_len(T_), each = N)
)

set.seed(1)
latent_factor_example_fit <- onlyif(
  run_extended_tests,
  dynamite(
    dformula = obs(y ~ 1, family = "gaussian") +
      lfactor() +
      splines(df = 10),
    data = latent_factor_example,
    group = "id",
    time = "time",
    iter = 2000,
    warmup = 1000,
    thin = 10,
    chains = 2,
    cores = 2,
    refresh = 0,
    save_warmup = FALSE,
    pars = c(
      "omega_alpha_1_y", "omega_raw_alpha_y", "omega_raw_psi",
      "omega_raw_psi_1_y", "L_lf", "lambda_raw_y"
    ),
    include = FALSE
  )
)

test_that("latent factor related parameters can be got", {
  skip_if_not(run_extended_tests)
  expect_equal(
    get_parameter_types(latent_factor_example_fit),
    c("alpha", "sigma", "lambda", "sigma_lambda", "psi", "tau_psi", "omega_psi")
  )
})

test_that("lambdas can be plotted", {
  skip_if_not(run_extended_tests)
  expect_error(
    plot_lambdas(latent_factor_example_fit),
    NA
  )
})

test_that("psis can be plotted", {
  skip_if_not(run_extended_tests)
  expect_error(
    plot_psis(latent_factor_example_fit),
    NA
  )
})

test_that("new group levels can't be included if model has latent factor", {
  skip_if_not(run_extended_tests)
  nd <- latent_factor_example
  nd$id[nd$id == 1] <- 100
  expect_error(
    predict(
      latent_factor_example_fit,
      newdata = nd, n_draws = 2
    ),
    paste(
      "Grouping variable `id` contains unknown levels:\nx Level \"100\"",
      "is not present in the original data\\.\ni Models with latent",
      "factors do not support new levels because of identifiability",
      "constraints\\."
    )
  )
})

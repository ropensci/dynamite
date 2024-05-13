#' @srrstats {BS7.4}. Tests are made that the input and fitted values are on a
#'   same scale.
#' @srrstats {G5.4} Predict and fitted method produce results identical with
#'   "manual" computation based on the same posterior samples.

data.table::setDTthreads(1) # For CRAN

test_that("predictions are on the same scale as input data", {
  expect_error(
    pred <- predict(gaussian_example_fit, type = "response", n_draws = 1),
    NA
  )
  expect_equal(sd(pred$y), sd(pred$y_new), tolerance = 0.5)
  expect_equal(mean(pred$y), mean(pred$y_new), tolerance = 0.5)

  expect_error(
    fit <- fitted(gaussian_example_fit, n_draws = 1),
    NA
  )
  expect_equal(sd(fit$y), sd(fit$y_fitted, na.rm = TRUE), tolerance = 0.5)
  expect_equal(
    mean(fit$y),
    mean(fit$y_fitted, na.rm = TRUE),
    tolerance = 0.5
  )
  expect_error(
    pred <- predict(multichannel_example_fit, type = "response", n_draws = 1),
    NA
  )
  expect_equal(mean(pred$g), mean(pred$g_new), tolerance = 0.5)
  expect_equal(mean(pred$p), mean(pred$p_new), tolerance = 0.5)
  expect_equal(mean(pred$b), mean(pred$b_new), tolerance = 0.1)
  expect_equal(sd(pred$g), sd(pred$g_new), tolerance = 0.5)
  expect_equal(sd(pred$p), sd(pred$p_new), tolerance = 0.5)
  expect_equal(sd(pred$b), sd(pred$b_new), tolerance = 0.1)
})

test_that("prediction works", {
  expect_error(
    predict(gaussian_example_fit, type = "response", n_draws = 2),
    NA
  )
  expect_error(
    predict(gaussian_example_fit, type = "mean", n_draws = 2),
    NA
  )
  expect_error(
    predict(gaussian_example_fit, type = "link", n_draws = 2),
    NA
  )
  expect_error(
    predict(categorical_example_fit, type = "response", n_draws = 2),
    NA
  )
  expect_error(
    predict(categorical_example_fit, type = "mean", n_draws = 2),
    NA
  )
  expect_error(
    predict(categorical_example_fit, type = "link", n_draws = 2),
    NA
  )
})

test_that("prediction works when starting from an arbitrary time point", {
  newdata <- gaussian_example |>
    dplyr::mutate(y = ifelse(time > 8, NA, y))

  set.seed(1)
  expect_error(
    pred1 <- predict(gaussian_example_fit, newdata = newdata, n_draws = 4),
    NA
  )
  set.seed(1)
  expect_error(
    pred2 <- predict(gaussian_example_fit,
      newdata = newdata |> dplyr::filter(time > 5),
      n_draws = 4
    ),
    NA
  )
  expect_equal(pred1 |> dplyr::filter(time > 5), pred2)
  fit <- gaussian_example_fit
  fit$data <- fit$data |>
    dplyr::mutate(time = time + 10)

  newdata <- fit$data |>
    dplyr::mutate(y = ifelse(time > 20, NA, y))
  set.seed(1)
  expect_error(
    pred1 <- predict(fit, newdata = newdata, n_draws = 4),
    NA
  )
  set.seed(1)
  expect_error(
    pred2 <- predict(fit,
      newdata = newdata |> dplyr::filter(time > 15),
      n_draws = 4
    ),
    NA
  )
  expect_equal(pred1 |> dplyr::filter(time > 15), pred2)
})


gaussian_example_single_fit <- get0(
  x = "gaussian_example_single_fit",
  envir = asNamespace("dynamite")
)

test_that("no groups prediction works", {
  expect_error(
    predict(gaussian_example_single_fit, type = "response", n_draws = 2),
    NA
  )
  expect_error(
    predict(gaussian_example_single_fit, type = "mean", n_draws = 2),
    NA
  )
  expect_error(
    predict(gaussian_example_single_fit, type = "link", n_draws = 2),
    NA
  )
})

test_that("fitted works", {
  expect_error(fitg <- fitted(gaussian_example_fit, n_draws = 1), NA)

  n <- ndraws(gaussian_example_fit) %/%
    gaussian_example_fit$stanfit@sim$chains
  idx <- gaussian_example_fit$permutation[1L]
  iter <- idx %% n
  chain <- 1 + idx %/% n
  xzy <- gaussian_example_fit$data |>
    dplyr::filter(id == 5 & time == 20)
  manual <- as_draws(gaussian_example_fit) |>
    dplyr::filter(.iteration == iter & .chain == chain) |>
    dplyr::summarise(fit = `alpha_y[20]` + nu_y_alpha_id5 +
      `delta_y_x[20]` * xzy$x + beta_y_z * xzy$z +
      `delta_y_y_lag1[20]` * xzy$y_lag1) |>
    dplyr::pull(fit)
  automatic <- fitg |>
    dplyr::filter(id == 5 & time == 20) |>
    dplyr::pull(y_fitted)
  expect_equal(automatic, manual)

  expect_error(fitc <- fitted(categorical_example_fit, n_draws = 1), NA)

  n <- ndraws(categorical_example_fit) %/%
    categorical_example_fit$stanfit@sim$chains
  idx <- categorical_example_fit$permutation[1L]
  iter <- idx %% n
  chain <- 1 + idx %/% n
  xzy <- categorical_example_fit$data |>
    dplyr::filter(id == 5 & time == 20)
  manual <- as_draws(categorical_example_fit) |>
    dplyr::filter(.iteration == iter & .chain == chain) |>
    dplyr::summarise(
      x_A = 0,
      x_B = alpha_x_B + beta_x_z_B * xzy$z + beta_x_x_lag1C_B +
        beta_x_y_lag1b_B,
      x_C = alpha_x_C + beta_x_z_C * xzy$z + beta_x_x_lag1C_C +
        beta_x_y_lag1b_C
    )
  manual <- (exp(manual) / sum(exp(manual)))[1, "x_C"]
  automatic <- fitc |>
    dplyr::filter(id == 5 & time == 20) |>
    dplyr::pull(x_fitted_C)
  expect_equal(automatic, manual)
})

test_that("categorical predict with type = link works", {
  expect_error(
    fitc <- predict(categorical_example_fit, type = "link", n_draws = 1),
    NA
  )
  n <- ndraws(categorical_example_fit) %/%
    categorical_example_fit$stanfit@sim$chains
  idx <- categorical_example_fit$permutation[1L]
  iter <- idx %% n
  chain <- 1 + idx %/% n
  xzy <- categorical_example_fit$data |>
    dplyr::filter(id == 5 & time == 2)
  manual <- as_draws(categorical_example_fit) |>
    dplyr::filter(.iteration == iter & .chain == chain) |>
    dplyr::summarise(x = alpha_x_C + beta_x_z_C * xzy$z + beta_x_x_lag1C_C) |>
    dplyr::pull(x)
  automatic <- fitc |>
    dplyr::filter(id == 5 & time == 2) |>
    dplyr::pull(x_link_C)
  expect_equal(automatic, manual)
})

test_that("fitted and predict give equal results for the first time point", {
  expect_equal(
    predict(gaussian_example_fit, type = "mean", n_draws = 2) |>
      dplyr::filter(time == 2) |>
      dplyr::pull("y_mean"),
    fitted(gaussian_example_fit, n_draws = 2) |>
      dplyr::filter(time == 2) |>
      dplyr::pull("y_fitted")
  )
  expect_equal(
    predict(multichannel_example_fit, type = "mean", n_draws = 2) |>
      dplyr::filter(time == 2) |>
      dplyr::select("g_mean", "p_mean", "b_mean"),
    fitted(multichannel_example_fit, n_draws = 2) |>
      dplyr::filter(time == 2) |>
      dplyr::select("g_fitted", "p_fitted", "b_fitted"),
    ignore_attr = TRUE
  )
})

test_that("predict with NA-imputed newdata works as default NULL", {
  # gaussian example
  set.seed(1)
  pred1 <- predict(gaussian_example_fit, type = "mean", n_draws = 2)
  newdata <- gaussian_example_fit$data
  newdata$y[newdata$time > 1] <- NA
  set.seed(1)
  pred2 <- predict(
    gaussian_example_fit,
    type = "mean",
    n_draws = 2,
    newdata = newdata
  )
  expect_equal(
    pred1 |> dplyr::pull("y_mean"),
    pred2 |> dplyr::pull("y_mean")
  )
  # categorical example
  set.seed(1)
  pred1 <- predict(categorical_example_fit, type = "mean", n_draws = 2)
  newdata <- categorical_example_fit$data
  newdata$y[newdata$time > 1] <- NA
  newdata$x[newdata$time > 1] <- NA
  set.seed(1)
  pred2 <- predict(
    categorical_example_fit,
    type = "mean",
    n_draws = 2,
    newdata = newdata
  )
  expect_equal(
    pred1 |> dplyr::select(-c(y, x)),
    pred2 |> dplyr::select(-c(y, x))
  )
})

test_that("permuting newdata for predict does not alter results", {
  newdata <- gaussian_example_fit$data
  newdata$y[newdata$time > 1] <- NA
  set.seed(1)
  pred1 <- predict(
    gaussian_example_fit,
    type = "mean",
    n_draws = 2,
    newdata = newdata
  )
  newdata2 <- newdata[sample(seq_len(nrow(newdata))), ]
  set.seed(1)
  expect_error(
    pred2 <- predict(
      gaussian_example_fit,
      type = "mean",
      n_draws = 2,
      newdata = newdata2
    ),
    NA
  )
  expect_equal(pred1, pred2)
})

test_that("factor time and integer time for predict give equal results", {
  newdata <- gaussian_example_fit$data
  newdata$time <- factor(newdata$time)
  set.seed(1)
  pred1 <- predict(
    gaussian_example_fit,
    type = "mean",
    n_draws = 2,
    newdata = newdata
  )
  newdata2 <- gaussian_example_fit$data
  set.seed(1)
  expect_error(
    pred2 <- predict(
      gaussian_example_fit,
      type = "mean",
      n_draws = 2,
      newdata = newdata2
    ),
    NA
  )
  expect_equal(pred1$simulated, pred2$simulated)
})

test_that("no groups fitted works", {
  expect_error(fitted(gaussian_example_single_fit, n_draws = 2), NA)
})

test_that("missing factor levels are restored", {
  categorical_example_noC <- categorical_example |>
    dplyr::mutate(x = dplyr::recode(x, "C" = "B")) |>
    dplyr::filter(time < 5)
  pred <- predict(
    categorical_example_fit,
    newdata = categorical_example_noC,
    ndraws = 2
  )
  expect_equal(levels(pred$x), c("A", "B", "C"))
})

test_that("new group levels can be included in newdata", {
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
    predict(gaussian_example_fit,
      newdata = gaussian_example_new_levels,
      type = "response", n_draws = 2, new_levels = "bootstrap"
    ),
    NA
  )
  expect_error(
    predict(gaussian_example_fit,
      newdata = gaussian_example_new_levels,
      type = "response", n_draws = 2, new_levels = "gaussian"
    ),
    NA
  )
  expect_error(
    predict(gaussian_example_fit,
      newdata = gaussian_example_new_levels,
      type = "response", n_draws = 2, new_levels = "original"
    ),
    NA
  )
})

test_that("imputation works", {
  set.seed(0)
  mis_x <- sample(seq_len(nrow(gaussian_example)), 150, replace = FALSE)
  mis_z <- sample(seq_len(nrow(gaussian_example)), 150, replace = FALSE)
  gaussian_example_impute <- gaussian_example
  gaussian_example_impute[mis_x, "x"] <- NA
  gaussian_example_impute[mis_z, "z"] <- NA
  expect_error(
    predict(
      gaussian_example_fit,
      newdata = gaussian_example_impute,
      type = "response",
      n_draws = 2L,
      impute = "locf"
    ),
    NA
  )
  expect_error(
    predict(
      gaussian_example_fit,
      newdata = gaussian_example_impute,
      type = "response",
      n_draws = 2L,
      impute = "nocb"
    ),
    NA
  )
})

test_that("global_fixed options produce equal results with balanced data", {
  set.seed(3)
  pred1 <- predict(gaussian_example_fit, n_draws = 2L, global_fixed = TRUE)
  set.seed(3)
  pred2 <- predict(gaussian_example_fit, n_draws = 2L, global_fixed = FALSE)
  expect_equal(pred1, pred2)
})

test_that("summarising via funs is equivalent to manual summary", {
  # type = "response"
  set.seed(1)
  pred1 <- predict(
    gaussian_example_fit,
    funs = list(y = list(mean = mean, sd = sd)),
    n_draws = 2L
  )
  pred1 <- pred1$simulated |> dplyr::filter(time > 1)
  set.seed(1)
  pred2 <- predict(
    gaussian_example_fit, n_draws = 2L, expand = FALSE
  )
  pred2 <- pred2$simulated |>
    dplyr::group_by(time, .draw) |>
    dplyr::filter(time > 1) |>
    dplyr::summarise(mean_y = mean(y_new), sd_y = sd(y_new)) |>
    dplyr::arrange(.draw)
  expect_equal(pred1$mean_y, pred2$mean_y)
  expect_equal(pred1$sd_y, pred2$sd_y)
  # type = "mean"
  set.seed(1)
  pred3 <- predict(
    gaussian_example_fit,
    type = "mean",
    funs = list(y = list(mean = mean, sd = sd)),
    n_draws = 2L
  )
  pred3 <- pred3$simulated |> dplyr::filter(time > 1)
  set.seed(1)
  pred4 <- predict(
    gaussian_example_fit, type = "mean", n_draws = 2L, expand = FALSE
  )
  pred4 <- pred4$simulated |>
    dplyr::group_by(time, .draw) |>
    dplyr::filter(time > 1) |>
    dplyr::summarise(mean_y = mean(y_mean), sd_y = sd(y_mean)) |>
    dplyr::arrange(.draw)
  expect_equal(pred3$mean_y, pred4$mean_y)
  expect_equal(pred3$sd_y, pred4$sd_y)
})

test_that("predict with loglik works", {
  out <- expect_error(
    initialize_predict(
      gaussian_example_fit,
      newdata = NULL,
      type = "mean",
      eval_type = "loglik",
      funs = list(),
      impute = "none",
      new_levels = "none",
      global_fixed = FALSE,
      idx_draws = 1:ndraws(gaussian_example_fit),
      expand = FALSE,
      df = TRUE
    )$simulated,
    NA
  )

  # n <- ndraws(gaussian_example_fit) %/%
  #   gaussian_example_fit$stanfit@sim$chains
  # idx <- gaussian_example_fit$permutation[1L]
  # iter <- idx %% n
  # chain <- 1 + idx %/% n
  xzy <- gaussian_example_fit$data |>
    dplyr::filter(id == 5 & time == 20)
  manual <- as_draws(gaussian_example_fit) |>
    dplyr::filter(.iteration == 1 & .chain == 1) |>
    dplyr::summarise(
      loglik = dnorm(
        xzy$y, `alpha_y[20]` +
          nu_y_alpha_id5 + `delta_y_x[20]` * xzy$x +
          beta_y_z * xzy$z + `delta_y_y_lag1[20]` * xzy$y_lag1,
        sigma_y,
        log = TRUE
      )
    ) |>
    dplyr::pull(loglik)
  automatic <- out |>
    dplyr::filter(id == 5 & time == 20 & .draw == 1) |>
    dplyr::pull(y_loglik)
  expect_equal(manual, automatic)
})

test_that("thin works", {
  expect_error(
    pred <- predict(gaussian_example_fit, thin = 10),
    NA
  )
  expect_equal(
    unique(pred$.draw),
    1L:20L
  )
  expect_error(
    pred <- fitted(gaussian_example_fit, thin = 10),
    NA
  )
  expect_equal(
    unique(pred$.draw),
    1L:20L
  )
})

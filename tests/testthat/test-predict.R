#' @srrstats {BS7.4}. Tests are made that the input and fitted values are on a
#'   same scale.

test_that("scale of predictions are on the same scale as input data", {
  expect_error(pred <- predict(gaussian_example_fit,
    type = "response", n_draws = 1), NA)
  expect_equal(sd(pred$y), sd(pred$y_new), tolerance = 0.5)
  expect_equal(mean(pred$y), mean(pred$y_new), tolerance = 0.5)

  expect_error(pred <- fitted(gaussian_example_fit, n_draws = 1), NA)
  expect_equal(sd(pred$y), sd(pred$y_fitted, na.rm = TRUE), tolerance = 0.5)
  expect_equal(mean(pred$y), mean(pred$y_fitted, na.rm = TRUE),
    tolerance = 0.5)

  expect_error(pred <- predict(multichannel_example_fit,
    type = "response", n_draws = 1), NA)
  expect_equal(sd(pred$g), sd(pred$g_new), tolerance = 0.5)
  expect_equal(mean(pred$g), mean(pred$g_new), tolerance = 0.5)
  expect_equal(sd(pred$p), sd(pred$p_new), tolerance = 0.5)
  expect_equal(mean(pred$p), mean(pred$p_new), tolerance = 0.5)
  expect_equal(sd(pred$b), sd(pred$b_new), tolerance = 0.1)
  expect_equal(mean(pred$b), mean(pred$b_new), tolerance = 0.1)
})

test_that("prediction works", {
  expect_error(predict(gaussian_example_fit,
                       type = "response", n_draws = 2), NA)
  expect_error(predict(gaussian_example_fit,
                       type = "mean", n_draws = 2), NA)
  expect_error(predict(gaussian_example_fit,
                       type = "link", n_draws = 2), NA)
  expect_error(predict(categorical_example_fit,
                       type = "response", n_draws = 2), NA)
  expect_error(predict(categorical_example_fit,
                       type = "mean", n_draws = 2), NA)
  expect_error(predict(categorical_example_fit,
                       type = "link", n_draws = 2), NA)
})

gaussian_example_single_fit <- get0("gaussian_example_single_fit",
                                    envir = asNamespace("dynamite"))

test_that("no groups prediction works", {
  expect_error(predict(gaussian_example_single_fit,
                       type = "response", n_draws = 2), NA)
  expect_error(predict(gaussian_example_single_fit,
                       type = "mean", n_draws = 2), NA)
  expect_error(predict(gaussian_example_single_fit,
                       type = "link", n_draws = 2), NA)
})

test_that("fitted works", {
  expect_error(fitted(gaussian_example_fit, n_draws = 2), NA)
  expect_error(fitted(categorical_example_fit, n_draws = 2), NA)
})

test_that("no groups fitted works", {
  expect_error(fitted(gaussian_example_single_fit, n_draws = 2), NA)
})

test_that("new group levels can be included in newdata", {
  gaussian_example_new_levels <- rbind(
    gaussian_example,
    data.frame(y = c(0.5, rep(NA, 29L)),
               x = rnorm(30),
               z = rbinom(30, 1, 0.7),
               id = 226L, time = seq.int(1, 30))
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
  expect_error(
    predict(gaussian_example_fit,
            newdata = gaussian_example_new_levels,
            type = "response", n_draws = 2, new_levels = "bootstrap"),
    NA
  )
  expect_error(
    predict(gaussian_example_fit,
            newdata = gaussian_example_new_levels,
            type = "response", n_draws = 2, new_levels = "gaussian"),
    NA
  )
  expect_error(
    predict(gaussian_example_fit,
            newdata = gaussian_example_new_levels,
            type = "response", n_draws = 2, new_levels = "original"),
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
    predict(gaussian_example_fit, newdata = gaussian_example_impute,
            type = "response", n_draws = 2, impute = "locf"),
    NA
  )
})

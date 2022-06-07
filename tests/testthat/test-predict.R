test_that("prediction works", {
  expect_error(predict(gaussian_example_fit, type = "response", n_draws = 2), NA)
  expect_error(predict(gaussian_example_fit, type = "mean", n_draws = 2), NA)
  expect_error(predict(gaussian_example_fit, type = "link", n_draws = 2), NA)
})

test_that("fitted works", {
  expect_error(fitted(gaussian_example_fit, n_draws = 2), NA)
})

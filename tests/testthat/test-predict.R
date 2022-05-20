test_that("prediction works", {
  expect_error(predict(gaussian_example_fit, n_draws = 1), NA)
})

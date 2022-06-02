test_that("prediction works", {
  expect_error(predict(gaussian_example_fit, n_draws = 2), NA)
})

test_that("fitted works", {
  expect_error(fitted(gaussian_example_fit, n_draws = 2), NA)
})

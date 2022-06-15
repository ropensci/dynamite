# TODO put these back when fixed + 1 is used in Stan code generation
#test_that("prediction works", {
#  expect_error(predict(gaussian_example_fit,
#                       type = "response", n_draws = 2), NA)
#  expect_error(predict(gaussian_example_fit,
#                       type = "mean", n_draws = 2), NA)
#  expect_error(predict(gaussian_example_fit,
#                       type = "link", n_draws = 2), NA)
#})
#
#test_that("no groups prediction works", {
#  expect_error(predict(gaussian_example_single_fit,
#                       type = "response", n_draws = 2), NA)
#  expect_error(predict(gaussian_example_single_fit,
#                       type = "mean", n_draws = 2), NA)
#  expect_error(predict(gaussian_example_single_fit,
#                       type = "link", n_draws = 2), NA)
#})
#
#test_that("fitted works", {
#  expect_error(fitted(gaussian_example_fit, n_draws = 2), NA)
#})
#
#test_that("no groups fitted works", {
#  expect_error(fitted(gaussian_example_single_fit, n_draws = 2), NA)
#})
#

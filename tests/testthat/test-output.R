test_that("conversion to data.frame works", {
  expect_error(as.data.frame(gaussian_example_fit), NA)
})

test_that("coefficients can be extracted", {
  expect_error(coef(gaussian_example_fit), NA)
})

test_that("fit object can be printed", {
  expect_error(capture.output(print(gaussian_example_fit)), NA)
})

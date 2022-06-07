# TODO test custom priors
test_that("priors can be extracted", {
  expect_error(
    get_priors(gaussian_example_fit),
    NA
  )
})

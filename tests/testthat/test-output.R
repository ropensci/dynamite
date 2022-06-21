test_that("conversion to data.frame works", {
  expect_error(as.data.frame(gaussian_example_fit), NA)
})

test_that("coefficients can be extracted", {
  expect_error(
    coef(gaussian_example_fit, type = "beta"),
    NA
  )
  expect_error(
    coef(gaussian_example_fit, type = "delta", include_alpha = FALSE),
    NA
  )
  expect_error(
    coef(gaussian_example_fit, type = "delta", include_alpha = TRUE),
    NA
  )
})

test_that("fit object can be printed", {
  expect_error(
    capture.output(
      capture.output(
        print(gaussian_example_fit),
        type = "message"),
      type = "output"),
    NA)
})

test_that("betas can be plotted", {
  expect_error(
    plot_betas(gaussian_example_fit),
    NA
  )
})

test_that("deltas can be plotted", {
  expect_error(
    plot_deltas(gaussian_example_fit),
    NA
  )
})

test_that("nus can be plotted", {
  expect_error(
    plot_nus(gaussian_example_fit),
    NA
  )
})

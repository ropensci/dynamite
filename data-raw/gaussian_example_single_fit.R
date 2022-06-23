## code to create `gaussian_example_single_fit` object

set.seed(1)

gaussian_example_single_fit <- dynamite(
  obs(y ~ -1 + z + varying(~ x + lag(y)), family = "gaussian") +
    splines(df = 20),
  data = gaussian_example |>
    dplyr::filter(.data$id == 1) |>
    dplyr::select(!.data$id),
  time = "time", iter = 2000, chains = 2, cores = 2, refresh = 0
)

usethis::use_data(gaussian_example_single_fit, overwrite = TRUE, compress = "xz")

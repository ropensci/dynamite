# Code to create `gaussian_example_fit` object

library(dynamite)

# Note the very small number of post-warmup iterations due to the data size
# restrictions in CRAN.
set.seed(1)
gaussian_example_fit <- dynamite(
  dformula =
    obs(
      y ~ -1 + z + varying(~ x + lag(y)),
      family = "gaussian"
    ) +
      random() +
      splines(df = 20),
  data = gaussian_example,
  group = "id",
  time = "time",
  iter = 2000,
  warmup = 1000,
  thin = 5,
  chains = 2,
  cores = 2,
  refresh = 0,
  save_warmup = FALSE
)

usethis::use_data(gaussian_example_fit, overwrite = TRUE, compress = "xz")

# Code to create `gaussian_example_single_fit` object

# use only first id
d <- gaussian_example |> dplyr::filter(id == 1)

# convergence issues with the current setup but doesn't matter for tests
set.seed(1)
gaussian_example_single_fit <- dynamite(
  obs(y ~ -1 + z + varying(~ x + lag(y)), family = "gaussian") +
    splines(df = 20),
  data = d,
  time = "time",
  init = 0,
  iter = 1100,
  warmup = 1000,
  chains = 1,
  refresh = 0,
  save_warmup = FALSE
)

usethis::use_data(gaussian_example_single_fit,
  overwrite = TRUE,
  compress = "xz", internal = TRUE
)

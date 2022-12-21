# Code to create `latent_factor_example_fit` object

library(dynamite)

# Note the very small number of post-warmup iterations due to the data size
# restrictions in CRAN.
set.seed(1)
latent_factor_example_fit <- dynamite(
  dformula =
    obs(
      y ~ 1,
      family = "gaussian"
    ) +
      random() +
      lfactor() +
      splines(df = 10),
  data = latent_factor_example,
  group = "id",
  time = "time",
  iter = 2000,
  warmup = 1000,
  thin = 5,
  chains = 1,
  cores = 1,
  refresh = 0,
  save_warmup = FALSE
)

usethis::use_data(
  latent_factor_example_fit,
  overwrite = TRUE,
  compress = "xz"
)

# Code to create `latent_factor_example_fit` object

library(dynamite)

# Note the very small number of post-warmup iterations due to the data size
# restrictions in CRAN.
set.seed(1)
latent_factor_example_fit <- dynamite(
  dformula =
    obs(y ~ 1, family = "gaussian") +
      lfactor() +
      splines(df = 10),
  data = latent_factor_example,
  group = "id",
  time = "time",
  iter = 2000,
  warmup = 1000,
  thin = 10,
  chains = 2,
  cores = 2,
  refresh = 0,
  save_warmup = FALSE,
  pars = c(
    "omega_alpha_1_y", "omega_raw_alpha_y", "omega_raw_psi",
    "omega_raw_psi_1_y", "L_lf", "lambda_raw_y", "lambda_std_y"
  ),
  include = FALSE
)

usethis::use_data(
  latent_factor_example_fit,
  overwrite = TRUE,
  compress = "xz"
)

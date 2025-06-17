# Code to create `gaussian_example_fit` object

library("dynamite")

# Note the very small number of post-warmup iterations due to the data size
# restrictions in CRAN.
set.seed(1)
gaussian_example_fit <- dynamite(
  dformula =
    obs(y ~ -1 + z + varying(~ x + lag(y)) + random(~1), family = "gaussian") +
    random_spec() +
    splines(df = 20),
  data = gaussian_example,
  time = "time",
  group = "id",
  iter = 2000,
  warmup = 1000,
  thin = 10,
  chains = 2,
  cores = 2,
  refresh = 0,
  save_warmup = FALSE,
  pars = c(
    "omega_alpha_1_y", "omega_raw_alpha_y", "nu_raw", "nu", "L",
    "sigma_nu", "a_y"
  ),
  include = FALSE,
  backend = "rstan"
)

usethis::use_data(
  gaussian_example_fit,
  overwrite = TRUE,
  compress = "xz"
)

# Code to create `gaussian_example_single_fit` object

# use only first id
d <- gaussian_example |> dplyr::filter(id == 1)

# convergence issues with the current setup but doesn't matter for tests
set.seed(1)
gaussian_example_single_fit <- dynamite(
  dformula = obs(y ~ -1 + z + varying(~ x + lag(y)), family = "gaussian") +
    splines(df = 20),
  data = d,
  time = "time",
  init = 0,
  iter = 1100,
  warmup = 1000,
  chains = 1,
  refresh = 0,
  save_warmup = FALSE,
  backend = "rstan"
)

usethis::use_data(
  gaussian_example_single_fit,
  overwrite = TRUE,
  compress = "xz",
  internal = TRUE
)

# Code to create `multichannel_example_fit` object

library("dynamite")

# Note the very small number of post-warmup iterations due to the data size
# restrictions in CRAN.
set.seed(1)
multichannel_example_fit <- dynamite(
  dformula = obs(g ~ lag(g) + lag(logp), family = "gaussian") +
    obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
    obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
    aux(numeric(logp) ~ log(p + 1)),
  data = multichannel_example,
  time = "time",
  group = "id",
  verbose = FALSE,
  chains = 1,
  cores = 1,
  iter = 2000,
  warmup = 1000,
  init = 0,
  refresh = 0,
  thin = 5,
  save_warmup = FALSE
)

usethis::use_data(
  multichannel_example_fit,
  overwrite = TRUE,
  compress = "xz"
)

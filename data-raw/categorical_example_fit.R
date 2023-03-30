# Code to create `categorical_example_fit` object

library(dynamite)

set.seed(1)
categorical_example_fit <- dynamite(
  dformula = obs(x ~ z + lag(x) + lag(y), family = "categorical") +
    obs(y ~ z + lag(x) + lag(y), family = "categorical"),
  data = categorical_example,
  time = "time",
  group = "id",
  chains = 1,
  refresh = 0,
  thin = 5,
  save_warmup = FALSE
)

usethis::use_data(
  categorical_example_fit,
  overwrite = TRUE,
  compress = "xz"
)

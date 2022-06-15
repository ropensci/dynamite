set.seed(1)
library(dynamite)
f <- obs(x ~ z + lag(x) + lag(y), family = categorical()) +
  obs(y ~ z + lag(x) + lag(y), family = categorical())

categorical_example_fit <- dynamite(f, categorical_example, "id", "time",
                                    chains = 1, refresh = 0)

usethis::use_data(categorical_example_fit, overwrite = TRUE)

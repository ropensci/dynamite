## code to create `multichannel_example_fit` object

set.seed(1)
library(dynamite)
f <- obs(g ~ -1 + lag(g) + lag(logp) + varying(~x), family = gaussian()) +
  obs(p ~ x + lag(g) + lag(logp) + lag(b), family = poisson()) +
  obs(b ~ lag(b) * lag(logp) + lag(b) * x + lag(b) * lag(g),
      family = bernoulli()) +
  aux(numeric(logp) ~ log(p + 1)) + splines(df = D)


multichannel_example_fit <- dynamite(f, multichannel_example, "id", "time",
                 chains = 2, cores = 2, init = 0, refresh = 0)

usethis::use_data(multichannel_example_fit, overwrite = TRUE)





<!-- README.md is generated from README.Rmd. Please edit that file -->

# dynamite: Bayesian Modeling and Causal Inference for Multivariate Longitudinal Data <a href="https://docs.ropensci.org/dynamite/"><img src="man/figures/logo.png" align="right" height="139"/></a>

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/dynamite/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/dynamite/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/dynamite/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ropensci/dynamite?branch=main)
[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/554_status.svg)](https://github.com/ropensci/software-review/issues/554)
[![dynamite status
badge](https://ropensci.r-universe.dev/badges/dynamite)](https://ropensci.r-universe.dev)
[![dynamite CRAN
badge](https://www.r-pkg.org/badges/version/dynamite)](https://cran.r-project.org/package=dynamite)
<!-- badges: end -->

The `dynamite` [R](https://www.r-project.org/) package provides an
easy-to-use interface for Bayesian inference of complex panel (time
series) data comprising of multiple measurements per multiple
individuals measured in time via dynamic multivariate panel models
(DMPM). The main features distinguishing the package and the underlying
methodology from many other approaches are:

- Support for regular time-invariant effects, group-level random
  effects, and time-varying effects modeled via Bayesian P-splines.
- Joint modeling of multiple measurements per individual (multiple
  channels) based directly on the assumed data-generating process.
  Individual channels can be univariate or multivariate.
- Support for various distributions: Currently Gaussian, Multivariate
  Gaussian, Student t, Categorical, Ordered, Multinomial, Poisson,
  Bernoulli, Binomial, Negative Binomial, Gamma, Exponential, and Beta
  distributions are available, and these can be combined arbitrarily in
  multichannel models.
- Allows evaluating realistic long-term counterfactual predictions that
  take into account the dynamic structure of the model by efficient
  posterior predictive distribution simulation.
- Transparent quantification of parameter and predictive uncertainty due
  to a fully Bayesian approach.
- Various visualization methods including a method for drawing and
  producing a TikZ code of the directed acyclic graph (DAG) of the model
  structure.
- User-friendly and efficient R interface with state-of-the-art
  estimation via Stan. Both `rstan` and `cmdstanr` backends are
  supported, with both parallel chains and within-chain parallelization.

The `dynamite` package is developed with the support of the Research
Council of Finland grant 331817
([PREDLIFE](https://sites.utu.fi/predlife/en/)). For further information
on DMPMs and the `dynamite` package, see the related papers:

- Helske J. and Tikka S. (2024). Estimating Causal Effects from Panel
  Data with Dynamic Multivariate Panel Models. *Advances in Life Course
  Research*, 60, 100617. ([Journal
  version](https://doi.org/10.1016/j.alcr.2024.100617), [SocArXiv
  preprint](https://doi.org/10.31235/osf.io/mdwu5))
- Tikka S. and Helske J. (2025). `dynamite`: An R Package for Dynamic
  Multivariate Panel Models. *Journal of Statistical Software*, 115(5),
  1-42. ([Journal version](https://doi.org/10.18637/jss.v115.i05),
  [arXiv preprint](https://arxiv.org/abs/2302.01607))

## Installation

You can install the most recent stable version of `dynamite` from
[CRAN](https://cran.r-project.org/package=dynamite) or the development
version from [R-universe](https://r-universe.dev/search) by running one
the following lines:

``` r
install.packages("dynamite")
install.packages("dynamite", repos = "https://ropensci.r-universe.dev")
```

## Example

A single-channel model with time-invariant effect of `z`, time-varying
effect of `x`, lagged value of the response variable `y` and a
group-specific random intercepts:

``` r
set.seed(1)
library("dynamite")
gaussian_example_fit <- dynamite(
  obs(y ~ -1 + z + varying(~ x + lag(y)) + random(~1), family = "gaussian") +
    splines(df = 20),
  data = gaussian_example, time = "time", group = "id",
  iter = 2000, chains = 2, cores = 2, refresh = 0
)
```

Summary of the model:

``` r
print(gaussian_example_fit)
#> Model:
#>   Family   Formula                                       
#> y gaussian y ~ -1 + z + varying(~x + lag(y)) + random(~1)
#> 
#> Correlated random effects added for response(s): y
#> 
#> Data: gaussian_example (Number of observations: 1450)
#> Grouping variable: id (Number of groups: 50)
#> Time index variable: time (Number of time points: 30)
#> 
#> NUTS sampler diagnostics:
#> 
#> No divergences, saturated max treedepths or low E-BFMIs.
#> 
#> Smallest bulk-ESS: 651 (sigma_nu_y_alpha)
#> Smallest tail-ESS: 853 (sigma_nu_y_alpha)
#> Largest Rhat: 1.01 (delta_y_y_lag1[2])
#> 
#> Elapsed time (seconds):
#>         warmup sample
#> chain:1  7.646  5.113
#> chain:2  8.010  4.943
#> 
#> Summary statistics of the time- and group-invariant parameters:
#> # A tibble: 113 × 10
#>    variable      mean median     sd    mad      q5   q95  rhat ess_bulk ess_tail
#>    <chr>        <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl> <dbl>    <dbl>    <dbl>
#>  1 alpha_y[2]  0.0577 0.0579 0.0321 0.0327 0.00492 0.110 1.000    2093.    1577.
#>  2 alpha_y[3]  0.0900 0.0908 0.0437 0.0437 0.0180  0.161 0.999    2148.    1681.
#>  3 alpha_y[4]  0.165  0.164  0.0390 0.0392 0.102   0.229 1.000    2301.    1644.
#>  4 alpha_y[5]  0.263  0.263  0.0377 0.0385 0.200   0.325 1.00     2342.    1505.
#>  5 alpha_y[6]  0.303  0.302  0.0365 0.0366 0.244   0.364 0.999    2439.    1573.
#>  6 alpha_y[7]  0.332  0.333  0.0399 0.0406 0.268   0.398 1.000    2123.    1573.
#>  7 alpha_y[8]  0.421  0.421  0.0364 0.0382 0.361   0.480 1.00     2160.    1475.
#>  8 alpha_y[9]  0.457  0.456  0.0400 0.0395 0.391   0.522 1.00     2090.    1797.
#>  9 alpha_y[10] 0.412  0.412  0.0421 0.0410 0.340   0.482 1.000    2159.    1760.
#> 10 alpha_y[11] 0.404  0.404  0.0404 0.0412 0.337   0.468 0.999    2036.    1756.
#> # ℹ 103 more rows
```

Posterior estimates of time-varying effects:

``` r
plot(gaussian_example_fit, types = c("alpha", "delta"), scales = "free")
```

<img src="man/figures/README-posterior-1.png" alt="" style="display: block; margin: auto;" />

And group-specific intercepts (for first 10 groups):

``` r
plot(gaussian_example_fit, types = "nu", groups = 1:10)
```

<img src="man/figures/README-random-1.png" alt="" style="display: block; margin: auto;" />

Traceplots and density plots for time-invariant parameters:

``` r
plot(gaussian_example_fit, plot_type = "trace", types = "beta")
```

<img src="man/figures/README-trace-1.png" alt="" style="display: block; margin: auto;" />

Posterior predictive samples for the first 4 groups (using the samples
based on the posterior distribution of the model parameters and observed
data on the first time point):

``` r
library("ggplot2")
pred <- predict(gaussian_example_fit, n_draws = 100)
pred |>
  dplyr::filter(id < 5) |>
  ggplot(aes(time, y_new, group = .draw)) +
  geom_line(alpha = 0.25) +
  # observed values
  geom_line(aes(y = y), colour = "tomato") +
  facet_wrap(~id) +
  theme_bw()
```

<img src="man/figures/README-prediction-1.png" alt="" style="display: block; margin: auto;" />

Visualizing the model structure as a DAG (a snapshot at time `t`):

``` r
plot(gaussian_example_fit, plot_type = "dag", show_covariates = TRUE)
```

<img src="man/figures/README-dag-1.png" alt="" style="display: block; margin: auto;" />

For more examples, see the package vignettes and the [blog post about
dynamite](https://ropensci.org/blog/2023/01/31/dynamite-r-package/).

## Related packages

- The `dynamite` package uses Stan via
  [`rstan`](https://CRAN.R-project.org/package=rstan) and
  [`cmdstanr`](https://mc-stan.org/cmdstanr/) (see also
  <https://mc-stan.org>), which is a probabilistic programming language
  for general Bayesian modelling.
- The [`brms`](https://CRAN.R-project.org/package=brms) package also
  uses Stan, and can be used to fit various complex multilevel models.
- Regression modeling with time-varying coefficients based on kernel
  smoothing and least squares estimation is available in package
  [`tvReg`](https://CRAN.R-project.org/package=tvReg). The
  [`tvem`](https://CRAN.R-project.org/package=tvem) package provides
  similar functionality for gaussian, binomial and poisson responses
  with [`mgcv`](https://CRAN.R-project.org/package=mgcv) backend.
- [`plm`](https://CRAN.R-project.org/package=plm) contains various
  methods to estimate linear models for panel data, e.g., fixed effect
  models.
- [`lavaan`](https://CRAN.R-project.org/package=lavaan) provides tools
  for structural equation modeling, and as such can be used to model
  various panel data models as well.

## Contributing

Contributions are very welcome, see
[CONTRIBUTING.md](https://github.com/ropensci/dynamite/blob/main/.github/CONTRIBUTING.md)
for general guidelines.

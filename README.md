
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dynamite

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/santikka/dynamite/workflows/R-CMD-check/badge.svg)](https://github.com/santikka/dynamite/actions)
[![Codecov test
coverage](https://codecov.io/gh/santikka/dynamite/branch/main/graph/badge.svg)](https://app.codecov.io/gh/santikka/dynamite?branch=main)
<!-- badges: end -->

The `dynamite` package provides easy-to-use interface for Bayesian
inference of complex panel data. The main features distinguishing the
package and the underlying methodology from many other approaches are:

-   Support for both time-varying and time-invariant effects.
-   Joint modeling of multiple measurements per individual (multiple
    channels).
-   Support for non-gaussian observations.
-   Realistic counterfactual predictions which take into account the
    dynamic structure of the model.
-   Clear quantification of parameter and predictive uncertainty due to
    a Bayesian approach.
-   User-friendly and efficient R interface with state-of-the-art
    estimation via Stan.

The `dynamite` package is developed with the support of Academy of
Finland grant 331817.

## Installation

You can install the development version of dynamite from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("santikka/dynamite")
```

## Example

A single-channel model with time-invariant effect of `z`, time-varying
effect of `x`, lagged value of the response variable `y` and a
group-specific random intercepts:

``` r
set.seed(1)
library(dynamite)
gaussian_example_fit <- dynamite(
  obs(y ~ -1 + z + varying(~ x + lag(y)), family = "gaussian",
      random_intercept = TRUE) + splines(df = 20),
  data = gaussian_example, time = "time", group = "id",
  iter = 2000, warmup = 1000, thin = 5,
  chains = 2, cores = 2, refresh = 0, save_warmup = FALSE
)
```

Posterior estimates of the fixed effects:

``` r
plot_betas(gaussian_example_fit)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="50%" />

Posterior estimates of time-varying effects

``` r
plot_deltas(gaussian_example_fit, scales = "free")
#> Warning: Removed 1 row(s) containing missing values (geom_path).
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="50%" />

And group-specific intercepts:

``` r
plot_nus(gaussian_example_fit)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="50%" />

## Related packages

The `dynamite` package uses Stan via `rstan` (see
<https://mc-stan.org>).

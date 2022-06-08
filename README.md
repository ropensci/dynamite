
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dynamite

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/santikka/dynamite/workflows/R-CMD-check/badge.svg)](https://github.com/santikka/dynamite/actions)
[![Codecov test
coverage](https://codecov.io/gh/santikka/dynamite/branch/main/graph/badge.svg)](https://app.codecov.io/gh/santikka/dynamite?branch=main)
<!-- badges: end -->

The `dynamite` package provides easy-to-use interface for Bayesian
inference of complex panel data. The main features distinguishing the
package and the underlying methodology from many other approaches are

-   Support for both time-varying and time-invariant effects
-   Joint modelling of multiple measurements per individual (multiple
    channels)
-   Support for non-gaussian observations
-   Realistic counterfactual predictions which take account the dynamic
    structure of the model
-   Clear quantification of parameter and predictive uncertainty due to
    Bayesian approach
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

TODO

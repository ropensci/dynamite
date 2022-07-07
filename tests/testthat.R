#' @srrstats {G3.0} Floating point numbers are not compared for equality
#' @srrstats {G5.2, G5.2a, G5.2b, G5.3, RE7.2, RE7.3}
#'   Demonstrated by package tests.
#' @srrstats {G5.12} Instructions on how to run extended tests is provided in
#'   tests/README.md.
#' @srrstats {G5.5, G5.6b} Fixed and random seeds are used appropriately in
#'   tests.
#' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstats {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*

library(testthat)
library(dynamite)

test_check("dynamite")

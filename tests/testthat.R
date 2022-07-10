#' @srrstats {G3.0} Floating point numbers are not compared for equality
#' @srrstats {G5.2, G5.2a, G5.2b, G5.3, RE7.2, RE7.3}
#'   Demonstrated by package tests.
#' @srrstats {G5.12} Instructions on how to run extended tests is provided in
#'   tests/README.md.
#' @srrstats {G5.5, G5.6b} Fixed and random seeds are used appropriately in
#'   tests.
library(testthat)
library(dynamite)

test_check("dynamite")

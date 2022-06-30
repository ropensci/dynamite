#' @srrstats {G3.0} Floating point numbers are not compared for equality
#' @srrstats {G5.0} *Where applicable or practicable, tests should use standard data sets with known properties (for example, the NIST Standard Reference Datasets, or data sets provided by other widely-used R packages).*
#' @srrstats {G5.2, G5.2a, G5.2b, G5.3, RE7.2, RE7.3}
#'   Demonstrated by package tests.
#' @srrstats {G5.4} **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets (potentially through comparisons using binding frameworks such as RStata)*
#' @srrstats {G5.4a} *For new methods, it can be difficult to separate out correctness of the method from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
#' @srrstats {G5.4b} *For new implementations of existing methods, correctness tests should include tests against previous implementations. Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, or use stored outputs from those where that is not possible.*
#' @srrstats {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available*
#' @srrstats {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstats {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
#' @srrstats {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstats {G5.6b} *Parameter recovery tests should be run with multiple random seeds when either data simulation or the algorithm contains a random component. (When long-running, such tests may be part of an extended, rather than regular, test suite; see G4.10-4.12, below).*
#' @srrstats {G5.7} **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*
#' @srrstats {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstats {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstats {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*
#' @srrstats {G5.10} *Extended tests should included and run under a common framework with other tests but be switched on by flags such as as a MYPKG_EXTENDED_TESTS=1 environment variable.*
#' @srrstats {G5.11} *Where extended tests require large data sets or other assets, these should be provided for downloading and fetched as part of the testing workflow.*
#' @srrstats {G5.11a} *When any downloads of additional data necessary for extended tests fail, the tests themselves should not fail, rather be skipped and implicitly succeed with an appropriate diagnostic message.*
#' @srrstats {G5.12} *Any conditions necessary to run extended tests such as platform requirements, memory, expected runtime, and artefacts produced that may need manual inspection, should be described in developer documentation such as a CONTRIBUTING.md or tests/README.md file.*
#' @srrstats {BS7.0} *Software should demonstrate and confirm recovery of parametric estimates of a prior distribution*
#' @srrstats {BS7.1} *Software should demonstrate and confirm recovery of a prior distribution in the absence of any additional data or information*
#' @srrstats {BS7.2} *Software should demonstrate and confirm recovery of a expected posterior distribution given a specified prior and some input data*
#' @srrstats {BS7.3} *Bayesian software should include tests which demonstrate and confirm the scaling of algorithmic efficiency with sizes of input data.*
#' @srrstats {BS7.4} *Bayesian software should implement tests which confirm that predicted or fitted values are on (approximately) the same scale as input values.*
#' @srrstats {BS7.4a} *Bayesian software should implement tests which confirm that predicted or fitted values are on (approximately) the same scale as input values.*

library(testthat)
library(dynamite)

test_check("dynamite")

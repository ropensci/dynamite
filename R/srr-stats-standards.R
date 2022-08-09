#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#'
#' @noRd
NULL

#' NA_standards
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#' @srrstatsNA {G1.0} Related academic paper to be published.
#' @srrstatsNA {G1.5, G1.6} No performance claims are made.
#' @srrstatsNA {G3.1} No covariance calculations are done.
#' @srrstatsNA {G3.1a} There are no covariance methods.
#' @srrstatsNA {G4.0} Output to local files is not supported.
#' @srrstatsNA {BS1.3a, BS2.8} This is not supported by the package.
#' @srrstatsNA {BS4.1} There are no sampler comparisons.
#' @srrstatsNA {BS1.4, BS1.5, BS4.3, BS4.4, BS4.5, BS4.6, BS4.7, BS5.4}
#'   Sampling is done in Stan, which has no automatic convergence checkers.
#' @srrstatsNA {BS2.10, BS2.11} Setting of seeds and starting values is handled
#'   by appropriate arguments to `dynamite` which are passed to
#'   `rstan::sampling`.
#' @srrstatsNA {RE2.3} No automatic data transformations are carried out in
#'   this sense, naturally the user can choose to center their data beforehand.
#' @srrstatsNA {RE3.2, RE3.3} Not applicable to Markov chain Monte Carlo
#'   algorithms.
#' @srrstatsNA {RE4.6} Not applicable to a Bayesian model, although posterior
#'   correlations of any parameters can be manually computed based on the
#'   extracted posterior samples.
#' @srrstatsNA {RE4.7, RE4.10, RE4.11, RE4.12} Not applicable.
#' @srrstatsNA {RE4.15, RE6.3, RE4.14, RE7.4} Forecasting is not supported
#' @srrstatsNA {RE7.0, RE7.0a, RE7.1, RE7.1a} Not relevant to the package.
#' @srrstatsNA {RE6.2} Such a plot as a default could be s misleading due to the
#'   lagged dependency structures and overly complicated because there are
#'   multiple channels with various distributional assumptions
#'   (i.e., discrete vs continuous).
#' @srrstatsNA {G5.11, G5.11a} No large data sets are used in the extended
#'   tests.
#' @srrstatsNA {BS7.4a} No assumptions on the scale of the input data is made
#'   and the scales should not matter as long as numerical issues are not
#'   encountered.
#' @srrstatsNA {G5.9, G5.9a, G5.9b} Due to the Monte Carlo variation in the
#'   results, it is not possible to test the explicit effect of negligible
#'   noise in the data.
#' @noRd
NULL

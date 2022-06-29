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
#' @srrstatsNA {G5.11a} Data for tests is not downloaded.
#' @srrstatsNA {BS1.0} The term "Hyperparameter" is not used.
#' @srrstatsNA {BS1.3a, BS2.8} This is not supported by the package.
#' @srrstatsNA {BS2.10, BS2.11} Setting of seeds and starting values
#'   is handled by Stan.
#' @srrstatsNA {BS2.15} Errors can be caught using base R functionalities.
#' @srrstatsNA {BS4.1} There are no sampler comparisons.
#' @srrstatsNA {BS1.4, BS1.5, BS4.3, BS4.4, BS4.5, BS4.6, BS4.7, BS5.4}
#'   Sampling is done in Stan, which has no convergence checkers.
#' @srrstatsNA {RE2.3} No automatic data transformations are carried out in
#'   this sense, naturally the user can choose to center their data beforehand.
#' @srrstatsNA {RE3.2, RE3.3}
#'   Sampling is done in Stan, which has no convergence checkers.
#' @srrstatsNA {RE4.6} Not applicable to a Bayesian model.
#' @srrstatsNA {RE4.7, RE4.10, RE4.11, RE4.12} Not applicable.
#' @srrstatsNA {RE4.15, RE6.3, RE4.14, RE7.4} Forecasting is not supported
#' @srrstatsNA {RE7.0, RE7.0a, RE7.1, RE7.1a} Not relevant to the package.
#' @srrstatsNA {RE6.2} We feel that such a plot as a default could be
#'   misleasing because there are multiple channels and various parameter types
#'   that can potentially be plotted.
#' @srrstatsNA {BS6.2, BS6.3, BS6.5} These can be accomplished with
#'   ggplot and bayesplot, for example.
#' @noRd
NULL

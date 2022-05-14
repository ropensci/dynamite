#' The 'dynamite' package.
#'
#' @description A DESCRIPTION OF THE PACKAGE
#'
#' @docType package
#' @name dynamite-package
#' @import methods
#' @importFrom rstan sampling
#' @importFrom posterior as_draws as_draws_df
#' @importFrom rlang .data :=
#' @importFrom stats drop.terms formula model.matrix model.matrix.lm na.exclude
#' @importFrom stats quantile coef setNames as.formula fitted na.pass terms
#' @importFrom stats plogis rbinom reformulate rnbinom rnorm rpois runif sd

#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan.
#' R package version 2.21.2. https://mc-stan.org
#'
NULL
#' Simulated Data of Gaussian Responses
#'
#' A simulated data containing gaussian response variables with two covariates.
#' The data was generated from a model with time-varying effects of
#' covariate x and the lagged value of the response variable, time-varying
#' intercept, and time-invariant effect of covariate z. The time-varying
#' coefficients vary according to a spline with 20 degrees of freedom.
#'
#' @source The data was generated according to a script in
#' \url{https://github.com/santikka/dynamite/blob/main/data-raw/gaussian_example.R}
#' @format A data frame with 3000 rows and 5 variables:
#' \describe{
#'   \item{y}{The response variable}
#'   \item{x}{Single continuous covariate}
#'   \item{z}{Single binary covariate}
#'   \item{id}{Variable defining individuals (1 to 100)}
#'    \item{time}{Variable defining the time point of the measurement (1 to 30)}
#' }
"gaussian_example"
#' Model Fit for the Simulated Data of Gaussian Responses
#'
#' A `dynamitefit` object obtained by running a `dynamite` on the
#' `gaussian_example` dataset as
#' \preformatted{
#' set.seed(1)
#' gaussian_example_fit <- dynamite(
#'   obs(y ~ -1 + z + varying(~ x + lag(y)), family = gaussian()) +
#'     splines(df = 20),
#'   data = gaussian_example, time = time, group = id,
#'   iter = 2000, chains = 2, cores = 2, refresh = 0)
#' }
#'
#' @source The data was generated according to a script in
#' \url{https://github.com/santikka/dynamite/blob/main/data-raw/gaussian_example_fit.R}
#' @format A `dynamitefit` object.
"gaussian_example_fit"

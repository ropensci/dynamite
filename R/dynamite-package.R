#' The 'dynamite' package.
#'
#' @description A DESCRIPTION OF THE PACKAGE
#'
#' @docType package
#' @name dynamite-package
#' @import methods
#' @importFrom rstan sampling
#' @importFrom stats drop.terms formula model.matrix model.matrix.lm na.exclude
#' na.pass plogis rbinom reformulate rnbinom rnorm rpois runif sd terms
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan.
#' R package version 2.21.2. https://mc-stan.org
#'
NULL
#' Simulated Data of Gaussian Responses
#'
#' A simulated data containing gaussian response variables with one covariate.
#' The data was generated from a model with intercept, effect of
#' covariate x and the lagged value of the response variable, with all effects
#' being time-varying according to a spline with 20 degrees of freedom.
#'
#' @format A data frame with 3000 rows and 4 variables:
#' \describe{
#'   \item{y}{The response variable}
#'   \item{x}{Single continuous covariate}
#'   \item{id}{Variable defining individuals (1 to 100)}
#'    \item{time}{Variable defining the time point of the measurement (1 to 30)}
#' }
"gaussian_example"

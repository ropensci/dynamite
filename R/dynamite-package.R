#' The 'dynamite' package.
#'
#' @description A DESCRIPTION OF THE PACKAGE TODO
#'
#' @docType package
#' @name dynamite-package
#' @import methods
#' @import data.table
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom rstan sampling
#' @importFrom posterior as_draws as_draws_df
#' @importFrom rlang .data caller_env
#' @importFrom stats drop.terms formula model.matrix model.matrix.lm na.exclude
#' @importFrom stats quantile coef setNames as.formula fitted na.pass terms
#' @importFrom stats plogis rbinom reformulate rnbinom rnorm rpois runif sd
#'
#' @srrstats {G5.1} *Data sets created within, and used to test, a package should be exported (or otherwise made generally available) so that users can confirm tests and run examples.*
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
#'   \item{y}{The response variable.}
#'   \item{x}{A continuous covariate.}
#'   \item{z}{A binary covariate.}
#'   \item{id}{Variable defining individuals (1 to 50).}
#'   \item{time}{Variable defining the time point of the measurement (1 to 30).}
#' }
"gaussian_example"
#' Model Fit for the Simulated Data of Gaussian Responses
#'
#' A `dynamitefit` object obtained by running a `dynamite` on the
#' `gaussian_example` dataset as
#' \preformatted{
#' set.seed(1)
#' library(dynamite)
#'
#' #' # note the small number of samples due to size restrictions on CRAN
#' gaussian_example_fit <- dynamite(
#'   obs(y ~ -1 + z + varying(~ x + lag(y)), family = gaussian(),
#'       random_intercept = TRUE) + splines(df = 20),
#'   data = gaussian_example, time = "time", group = "id",
#'   iter = 2000, warmup = 1000, thin = 5,
#'   chains = 2, cores = 2, refresh = 0, save_warmup = FALSE
#' )
#' }
#'
#' @source The data was generated according to a script in
#' \url{https://github.com/santikka/dynamite/blob/main/data-raw/gaussian_example_fit.R}
#' @format A `dynamitefit` object.
"gaussian_example_fit"
#' Simulated Multivariate Panel Data
#'
#' A simulated multichannel data containing multiple individuals with multiple
#' response variables of different distributions.
#'
#' @source The data was generated according to a script in
#' \url{https://github.com/santikka/dynamite/blob/main/data-raw/multichannel_example.R}
#' @format A data frame with 3000 rows and 5 variables:
#' \describe{
#'   \item{id}{Variable defining individuals (1 to 50).}
#'   \item{time}{Variable defining the time point of the measurement (1 to 20).}
#'   \item{g}{Response variable following gaussian distribution.}
#'   \item{p}{Response variable following Poisson distribution.}
#'   \item{z}{Response variable following bernoulli distribution.}

#' }
"multichannel_example"
#' Model Fit for the multichannel_example Data
#'
#'
#' A `dynamitefit` object obtained by running a `dynamite` on the
#' `multichannel_example` dataset as
#' \preformatted{
#' set.seed(1)
#' library(dynamite)
#' f <- obs(g ~ lag(g) + lag(logp), family = gaussian()) +
#'   obs(p ~ lag(g) + lag(logp) + lag(b), family = poisson()) +
#'   obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = bernoulli()) +
#'   aux(numeric(logp) ~ log(p + 1))
#'
#' # note the small number of samples due to size restrictions on CRAN
#' multichannel_example_fit <- dynamite(
#'   f, multichannel_example, "id", "time",
#'   chains = 1, cores = 1, iter = 2000, warmup = 1000, init = 0, refresh = 0,
#'   thin = 5, save_warmup = FALSE)
#' }
#'
#' @source Script in
#' \url{https://github.com/santikka/dynamite/blob/main/data-raw/multichannel_example_fit.R}
#' @format A `dynamitefit` object.
"multichannel_example_fit"
#' Simulated Multivariate Panel Data
#'
#' A simulated data containing multiple individuals with two categorical
#' response variables.
#'
#' @source The data was generated according to a script in
#' \url{https://github.com/santikka/dynamite/blob/main/data-raw/categorical_example.R}
#' @format A data frame with 2000 rows and 5 variables:
#' \describe{
#'   \item{id}{Variable defining individuals (1 to 100).}
#'   \item{time}{Variable defining the time point of the measurement (1 to 20).}
#'   \item{x}{Categorical variable with three levels, A, B, and C.}
#'   \item{x}{Categorical variable with three levels, a, b, and c.}
#'   \item{z}{A continuous covariate.}

#' }
"categorical_example"
#' Model Fit for the categorical_example Data
#'
#' A `dynamitefit` object obtained by running a `dynamite` on the
#' `categorical_example` dataset as
#' \preformatted{
#' set.seed(1)
#' library(dynamite)
#' f <- obs(x ~ z + lag(x) + lag(y), family = categorical()) +
#'   obs(y ~ z + lag(x) + lag(y), family = categorical())
#'
#' # note the small number of samples due to size restrictions on CRAN
#' categorical_example_fit <- dynamite(
#'   f, categorical_example, "id", "time",
#'   chains = 1, refresh = 0, thin = 5, save_warmup = FALSE)
#' }
#'
#' @source Script in
#' \url{https://github.com/santikka/dynamite/blob/main/data-raw/categorical_example_fit.R}
#' @format A `dynamitefit` object.
"categorical_example_fit"

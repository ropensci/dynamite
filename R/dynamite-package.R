#' The 'dynamite' package.
#'
#' @description Easy-to-use and efficient interface for Bayesian inference of
#' complex panel data consisting of multiple individuals with multiple
#' measurements over time. Supports several observational distributions,
#' time-varying effects and realistic counterfactual predictions which take into
#' account the dynamic structure of the model. See the README and the package
#' vignette for details.
#'
#' @docType package
#' @name dynamite-package
#' @importFrom checkmate test_int test_integer test_string test_numeric
#' @importFrom checkmate test_number test_flag test_logical test_character
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom data.table data.table as.data.table is.data.table
#' @importFrom data.table set setDT setDF setkey setkeyv :=
#' @importFrom rstan sampling
#' @importFrom posterior as_draws as_draws_df
#' @importFrom rlang .data caller_env
#' @importFrom stats drop.terms formula model.matrix model.matrix.lm na.exclude
#' @importFrom stats quantile coef setNames as.formula fitted na.pass terms
#' @importFrom stats plogis rbinom reformulate rnbinom rnorm rpois runif sd
#' @importFrom stats nobs
#'
#' @srrstats {G2.0, G2.0a, G2.1, G2.1a, G2.2, G2.3, G2.3a, G2.3b}
#'   Input types are asserted and appropriately restricted and tested
#' @srrstats {G2.7} `data.frame` and its extensions are supported.
#'   Input types are well defined and asserted.
#' @srrstats {G2.8} Package subfunctions receive well defined-inputs.
#' @srrstats {G1.3} Terminology is defined and explained
#' @srrstats {G1.4, G1.4a} roxygen2 is used.
#' @srrstats {G5.1} Package data is exported.
#' @srrstats {G2.10} Extraction of single columns is systematic and robust.
#' @srrstats {RE1.2} Documented across the package.
#' @srrstats {BS1.2, BS1.2a, BS1.2b, BS1.2c} Prior specification is documented.
#' @srrstats {BS1.0} The term "Hyperparameter" is not used.
#' @srrstats {BS2.15} Errors can be caught using base R functionalities.
#' @srrstats {BS2.10, BS2.11} Setting of seeds and starting values is handled
#'   by appropriate arguments to `dynamite` which are passed to
#'   `rstan::sampling`.
#' @srrstats {BS7.1, BS7.2} Parameters used to simulate example datasets are
#'   recovered.
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
#' gaussian_example_fit <- dynamite(
#'   obs(y ~ -1 + z + varying(~ x + lag(y)), family = "gaussian",
#'       random_intercept = TRUE) + splines(df = 20),
#'   data = gaussian_example, time = "time", group = "id",
#'   iter = 2000, warmup = 1000, thin = 5,
#'   chains = 2, cores = 2, refresh = 0, save_warmup = FALSE
#' )
#' }
#' Note the small number of samples due to size restrictions on CRAN.
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
#'   \item{b}{Response variable following Bernoulli distribution.}
#' }
"multichannel_example"
#' Model Fit for the multichannel_example Data
#'
#' A `dynamitefit` object obtained by running a `dynamite` on the
#' `multichannel_example` dataset as
#' \preformatted{
#' set.seed(1)
#' library(dynamite)
#' f <- obs(g ~ lag(g) + lag(logp), family = "gaussian") +
#'   obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
#'   obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
#'   aux(numeric(logp) ~ log(p + 1))
#' multichannel_example_fit <- dynamite(
#'   f, multichannel_example, "id", "time",
#'   chains = 1, cores = 1, iter = 2000, warmup = 1000, init = 0, refresh = 0,
#'   thin = 5, save_warmup = FALSE)
#' }
#' Note the small number of samples due to size restrictions on CRAN.
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
#'   \item{y}{Categorical variable with three levels, a, b, and c.}
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
#' f <- obs(x ~ z + lag(x) + lag(y), family = "categorical") +
#'   obs(y ~ z + lag(x) + lag(y), family = "categorical")
#' categorical_example_fit <- dynamite(
#'   f, categorical_example, "id", "time",
#'   chains = 1, refresh = 0, thin = 5, save_warmup = FALSE)
#' }
#' Note the small number of samples due to size restrictions on CRAN.
#' @source Script in
#' \url{https://github.com/santikka/dynamite/blob/main/data-raw/categorical_example_fit.R}
#' @format A `dynamitefit` object.
"categorical_example_fit"

#' The `dynamite` package.
#'
#' @description Easy-to-use and efficient interface for Bayesian inference of
#' complex panel data consisting of multiple individuals with multiple
#' measurements over time. Supports several observational distributions,
#' time-varying effects and realistic counterfactual predictions which take into
#' account the dynamic structure of the model.
#'
#' # See Also
#'
#' * The package vignettes
#' * [dynamite::dynamiteformula()] for information on defining models.
#' * [dynamite::dynamite()] for information on fitting models.
#' * <https://github.com/ropensci/dynamite/issues/> to submit a bug report
#'   or a feature request.
#'
#' # Authors
#' | Santtu Tikka (author) | <santtuth@@gmail.com>  |
#' | --------------------- | ---------------------- |
#' | Jouni Helske (author) | <jouni.helske@@iki.fi> |
#'
#' @name dynamite-package
#' @importFrom data.table :=
#' @importFrom loo loo
#' @importFrom patchwork wrap_plots
#' @importFrom posterior as_draws as_draws_df ndraws
#' @importFrom rlang .data
#' @importFrom stats formula model.matrix model.matrix.lm na.exclude
#' @importFrom stats quantile coef setNames as.formula fitted na.pass terms
#' @importFrom stats plogis rbinom reformulate rnbinom rnorm rpois runif sd
#' @importFrom stats nobs qlogis update
#'
#' @srrstats {G2.0, G2.0a, G2.1, G2.1a, G2.2, G2.3, G2.3a, G2.3b}
#'   Input types are asserted and appropriately restricted and tested
#' @srrstats {G2.7} `data.frame` and its extensions are supported.
#'   Input types are well defined and asserted.
#' @srrstats {G2.8} Package subfunctions receive well defined-inputs.
#' @srrstats {G1.3} Terminology is defined and explained.
#' @srrstats {G1.4, G1.4a} roxygen2 is used.
#' @srrstats {G5.1} Package data is exported.
#' @srrstats {G2.10} Extraction of single columns is systematic and robust.
#' @srrstats {RE1.2} Documented across the package.
#' @srrstats {BS1.2, BS1.2a, BS1.2b, BS1.2c} Prior specification is documented.
#' @srrstats {BS1.0} The term "Hyperparameter" is not used.
#' @srrstats {BS2.15} Errors can be caught using base R functionalities.
#' @srrstats {BS7.1, BS7.2} Parameters used to simulate example datasets are
#'   recovered.
"_PACKAGE"

#' Simulated Data of Gaussian Responses
#'
#' Simulated data containing gaussian response variables with two covariates.
#' The dataset was generated from a model with time-varying effects of
#' covariate x and the lagged value of the response variable, time-varying
#' intercept, and time-invariant effect of covariate z. The time-varying
#' coefficients vary according to a spline with 20 degrees of freedom.
#'
#' @family examples
#' @source The data was generated via `gaussian_example.R` in
#' <https://github.com/ropensci/dynamite/tree/main/data-raw/>
#' @format A data frame with 3000 rows and 5 variables:
#' \describe{
#'  \item{y}{The response variable.}
#'  \item{x}{A continuous covariate.}
#'  \item{z}{A binary covariate.}
#'  \item{id}{Variable defining individuals (1 to 50).}
#'  \item{time}{Variable defining the time point of the measurement (1 to 30).}
#' }
"gaussian_example"

#' Model Fit for the Simulated Data of Gaussian Responses
#'
#' A `dynamitefit` object obtained by running `dynamite` on the
#' `gaussian_example` dataset as
#' \preformatted{
#' set.seed(1)
#' library(dynamite)
#' gaussian_example_fit <- dynamite(
#'   obs(y ~ -1 + z + varying(~ x + lag(y)) + random(~1), family = "gaussian") +
#'     random_spec() + splines(df = 20),
#'   data = gaussian_example,
#'   time = "time",
#'   group = "id",
#'   iter = 2000,
#'   warmup = 1000,
#'   thin = 10,
#'   chains = 2,
#'   cores = 2,
#'   refresh = 0,
#'   save_warmup = FALSE,
#'   pars = c("omega_alpha_1_y", "omega_raw_alpha_y", "nu_raw", "nu", "L",
#'     "sigma_nu", "a_y"),
#'   include = FALSE
#' )
#' }
#' Note the very small number of samples due to size restrictions on CRAN.
#' @family examples
#' @source The data was generated via `gaussian_example_fit.R` in
#' <https://github.com/ropensci/dynamite/tree/main/data-raw/>
#' @format A `dynamitefit` object.
"gaussian_example_fit"

# #' Model Fit for the time-varying example in the `dynamite_simulation`
# #' Vignette
# #'
# #' A `dynamitefit` object obtained by running `dynamite` with the
# #' `"Fixed_param"` algorithm on the specified `inits` in the example.
# #' \preformatted{
# #' set.seed(1)
# #' library(dynamite)
# #' gaussian_simulation_fit <- dynamite(
# #'   dformula = f,
# #'   data = d,
# #'   time = "time",
# #'   group = "id",
# #'   chains = 1,
# #'   iter = 1,
# #'   algorithm = "Fixed_param",
# #'   init = list(init),
# #' )
# #' }
# #' @family examples
# #' @source The data was generated via to `gaussian_simulation_fit.R` in
# #' <https://github.com/ropensci/dynamite/tree/main/data-raw/>
# #' @format A `dynamitefit` object.
# "gaussian_simulation_fit"

#' Simulated Multivariate Panel Data
#'
#' A simulated multichannel data containing multiple individuals with multiple
#' response variables of different distributions.
#'
#' @family examples
#' @source The data was generated via `multichannel_example.R` in
#' <https://github.com/ropensci/dynamite/tree/main/data-raw/>
#' @format A data frame with 3000 rows and 5 variables:
#' \describe{
#'  \item{id}{Variable defining individuals (1 to 50).}
#'  \item{time}{Variable defining the time point of the measurement (1 to 20).}
#'  \item{g}{Response variable following gaussian distribution.}
#'  \item{p}{Response variable following Poisson distribution.}
#'  \item{b}{Response variable following Bernoulli distribution.}
#' }
"multichannel_example"

#' Model Fit for the Simulated Multivariate Panel Data
#'
#' A `dynamitefit` object obtained by running `dynamite` on the
#' `multichannel_example` dataset as
#' \preformatted{
#' set.seed(1)
#' library(dynamite)
#' f <- obs(g ~ lag(g) + lag(logp), family = "gaussian") +
#'   obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
#'   obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
#'   aux(numeric(logp) ~ log(p + 1))
#' multichannel_example_fit <- dynamite(
#'   f,
#'   data = multichannel_example,
#'   time = "time",
#'   group = "id",
#'   chains = 1,
#'   cores = 1,
#'   iter = 2000,
#'   warmup = 1000,
#'   init = 0,
#'   refresh = 0,
#'   thin = 5,
#'   save_warmup = FALSE
#' )
#' }
#' Note the small number of samples due to size restrictions on CRAN.
#' @family examples
#' @source THe data was generated via `multichannel_example_fit.R` in
#' <https://github.com/ropensci/dynamite/tree/main/data-raw/>
#' @format A `dynamitefit` object.
"multichannel_example_fit"

#' Simulated Categorical Multivariate Panel Data
#'
#' A simulated data containing multiple individuals with two categorical
#' response variables.
#'
#' @family examples
#' @source The data was generated via `categorical_example.R` in
#' <https://github.com/ropensci/dynamite/tree/main/data-raw/>
#' @format A data frame with 2000 rows and 5 variables:
#' \describe{
#'  \item{id}{Variable defining individuals (1 to 100).}
#'  \item{time}{Variable defining the time point of the measurement (1 to 20).}
#'  \item{x}{Categorical variable with three levels, A, B, and C.}
#'  \item{y}{Categorical variable with three levels, a, b, and c.}
#'  \item{z}{A continuous covariate.}
#' }
"categorical_example"

#' Model Fit for the Simulated Categorical Multivariate Panel Data
#'
#' A `dynamitefit` object obtained by running `dynamite` on the
#' `categorical_example` dataset as
#' \preformatted{
#' set.seed(1)
#' library(dynamite)
#' f <- obs(x ~ z + lag(x) + lag(y), family = "categorical") +
#'   obs(y ~ z + lag(x) + lag(y), family = "categorical")
#' categorical_example_fit <- dynamite(
#'   f,
#'   data = categorical_example,
#'   time = "time",
#'   group = "id",
#'   chains = 1,
#'   refresh = 0,
#'   thin = 5,
#'   save_warmup = FALSE
#' )
#' }
#' Note the small number of samples due to size restrictions on CRAN.
#' @family examples
#' @source The data was generated via `categorical_example_fit.R` in
#' <https://github.com/ropensci/dynamite/tree/main/data-raw/>
#' @format A `dynamitefit` object.
"categorical_example_fit"

# #' Simulated Latent Factor Model Panel Data
# #'
# #' A simulated single-channel data containing multiple individuals whose
# #' trajectories are defined by a latent factor and random intercept terms.
# #'
# #' @family examples
# #' @source The data was generated via `latent_factor_example.R` in
# #' <https://github.com/ropensci/dynamite/blob/main/data-raw/>
# #' @format A data frame with 2000 rows and 3 variables:
# #' \describe{
# #'  \item{y}{A continuos variable.}
# #'  \item{id}{Variable defining individuals (1 to 100).}
# #'  \item{time}{Variable defining the time point of the measurement
# #'  (1 to 20).}
# #' }
# "latent_factor_example"

# #' Model Fit for the Simulated Latent Factor Data
# #'
# #' A `dynamitefit` object obtained by running `dynamite` on the
# #' `latent_factor_example` dataset as
# #' \preformatted{
# #' set.seed(1)
# #' library(dynamite)
# #' latent_factor_example_fit <- dynamite(
# #'   obs(y ~ 1, family = "gaussian") + lfactor() + splines(df = 10),
# #'   data = latent_factor_example,
# #'   time = "time",
# #'   group = "id",
# #'   iter = 2000,
# #'   warmup = 1000,
# #'   thin = 10,
# #'   chains = 2,
# #'   cores = 2,
# #'   refresh = 0,
# #'   save_warmup = FALSE,
# #'   pars = c("omega_alpha_1_y", "omega_raw_alpha_y", "omega_raw_psi", "L_lf",
# #'     "lambda_raw_y", "lambda_std_y"),
# #'   include = FALSE
# #' )
# #' }
# #' Note the very small number of samples due to size restrictions on CRAN.
# #' @family examples
# #' @source The data was generated via `latent_factor_example_fit.R` in
# #' <https://github.com/ropensci/dynamite/tree/main/data-raw/>
# #' @format A `dynamitefit` object.
# "latent_factor_example_fit"

#' Extract Samples From a `dynamitefit` Object as a Data Frame
#'
#' Provides a `data.frame` representation of the posterior samples of the model
#' parameters.
#'
#' The arguments `responses` and `types` can be used to extract only a subset
#' of the model parameters (i.e., only certain types of parameters related to a
#' certain response variable).
#'
#' Potential values for the `types` argument are:
#'
#'  * `alpha`\cr Intercept terms (time-invariant or time-varying).
#'  * `beta`\cr Time-invariant regression coefficients.
#'  * `cutpoint`\cr Cutpoints for ordinal regression.
#'  * `delta`\cr Time-varying regression coefficients.
#'  * `nu`\cr Group-level random effects.
#'  * `lambda`\cr Factor loadings.
#'  * `psi`\cr Latent factors.
#'  * `tau`\cr Standard deviations of the spline coefficients of `delta`.
#'  * `tau_alpha`\cr Standard deviations of the spline coefficients of
#'    time-varying `alpha`.
#'  * `sigma_nu`\cr Standard deviations of the random effects `nu`.
#'  * `corr_nu`\cr Pairwise within-group correlations of random effects `nu`.
#'     Samples of the full correlation matrix can be extracted manually as
#'     `rstan::extract(fit$stanfit, pars = "corr_matrix_nu")` if necessary.
#'  * `sigma_lambda`\cr Standard deviations of the latent factor loadings
#'    `lambda`.
#'  * `corr_psi`\cr Pairwise correlations of the latent factors.
#'     Samples of the full correlation matrix can be extracted manually as
#'     `rstan::extract(fit$stanfit, pars = "corr_matrix_psi")` if necessary.
#'  * `sigma`\cr Standard deviations of gaussian responses.
#'  * `corr`\cr Pairwise correlations of multivariate gaussian responses.
#'  * `phi`\cr Describes various distributional parameters, such as:
#'    - Dispersion parameter of the Negative Binomial distribution.
#'    - Shape parameter of the Gamma distribution.
#'    - Precision parameter of the Beta distribution.
#'    - Degrees of freedom of the Student t-distribution.
#'  * `omega`\cr Spline coefficients of the regression coefficients `delta`.
#'  * `omega_alpha`\cr Spline coefficients of time-varying `alpha`.
#'  * `omega_psi`\cr Spline coefficients of the latent factors `psi`.
#'
#' @export
#' @family output
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param types \[`character()`]\cr Type(s) of the parameters for which the
#'   samples should be extracted. See details of possible values. Default is
#'   all values listed in details except spline coefficients `omega`.
#'   This argument is mutually exclusive with `parameters`.
#' @param parameters \[`character()`]\cr Parameter(s) for which the samples
#'   should be extracted. Possible options can be found with function
#'   `get_parameter_names()`. Default is all parameters of specific type for
#'   all responses. This argument is mutually exclusive with `types`.
#' @param responses \[`character()`]\cr Response(s) for which the samples
#'   should be extracted. Possible options are elements of
#'   `unique(x$priors$response)`, and the default is this entire vector.
#'    Ignored if the argument `parameters` is supplied.
#'   `omega_alpha`, and `omega_psi`. See also [dynamite::get_parameter_types()].
#' @param times \[`double()`]\cr Time point(s) to keep. If `NULL`
#'   (the default), all time points are kept.
#' @param groups \[`character()`] Group name(s) to keep. If `NULL`
#'   (the default), all groups are kept.
#' @param summary \[`logical(1)`]\cr If `TRUE`, returns posterior
#'   mean, standard deviation, and posterior quantiles (as defined by the
#'   `probs` argument) for all parameters. If `FALSE` (default), returns the
#'   posterior samples instead.
#' @param probs \[`numeric()`]\cr Quantiles of interest. Default is
#'   `c(0.05, 0.95)`.
#' @param include_fixed \[`logical(1)`]\cr If `TRUE` (default), time-varying
#'   parameters for `1:fixed` time points are included in the output as `NA`
#'   values. If `FALSE`, fixed time points are omitted completely
#'   from the output.
#' @param ... Ignored.
#' @return A `tibble` containing either samples or summary statistics of the
#'   model parameters in a long format. For a wide format, see
#'   [dynamite::as_draws()].
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' as.data.frame(
#'   gaussian_example_fit,
#'   responses = "y",
#'   types = "beta"
#' )
#'
#' # Basic summaries can be obtained automatically with summary = TRUE
#' as.data.frame(
#'   gaussian_example_fit,
#'   responses = "y",
#'   types = "beta",
#'   summary = TRUE
#' )
#'
#' # Time-varying coefficients "delta"
#' as.data.frame(
#'   gaussian_example_fit,
#'   responses = "y",
#'   types = "delta",
#'   summary = TRUE
#' )
#'
#' # Obtain summaries for a specific parameters
#' as.data.frame(
#'   gaussian_example_fit,
#'   parameters = c("tau_y_x", "sigma_y"),
#'   summary = TRUE
#' )
#'
as.data.frame.dynamitefit <- function(x, row.names = NULL, optional = FALSE,
                                      types = NULL, parameters = NULL,
                                      responses = NULL,
                                      times = NULL, groups = NULL,
                                      summary = FALSE, probs = c(0.05, 0.95),
                                      include_fixed = TRUE, ...) {
  out <- as.data.table.dynamitefit(
    x = x,
    keep.rownames = FALSE,
    row.names = row.names,
    optional = optional,
    parameters = parameters,
    responses = responses,
    types = types,
    times = times,
    groups = groups,
    summary = summary,
    probs = probs,
    include_fixed = include_fixed,
    ...
  )
  tibble::tibble(data.table::setDF(out))
}

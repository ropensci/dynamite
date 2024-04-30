#' Extract Regression Coefficients of a Dynamite Model
#'
#' Extracts either time-varying or time-invariant parameters of the model.
#'
#' @export
#' @family output
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param type  \[`character(1)`]\cr Either `beta` (the default) for
#'   time-invariant coefficients, `delta` for time-varying coefficients,
#'   `nu` for random effects, `lambda` for factor loadings, or `psi` for
#'   latent factor. Ignored if the argument `parameters` is supplied.
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), extracts also
#'   time-invariant intercept term alpha if time-invariant parameters beta are
#'   extracted, and time-varying alpha if time-varying delta are extracted.
#'   Ignored if the argument `parameters` is supplied.
#' @param summary \[`logical(1)`]\cr If `TRUE` (default), returns posterior
#'   mean, standard deviation, and posterior quantiles (as defined by the
#'   `probs` argument) for all parameters. If `FALSE`, returns the
#'   posterior samples instead.
#' @inheritParams as.data.frame.dynamitefit
#' @param ... Ignored.
#' @return A `tibble` containing either samples or summary statistics of the
#'   model parameters in a long format.
#' @srrstats {G2.3a, RE4.2} Provides model coefficients.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' betas <- coef(gaussian_example_fit, type = "beta")
#' deltas <- coef(gaussian_example_fit, type = "delta")
#'
coef.dynamitefit <- function(object,
                             parameters = NULL,
                             type = c("beta", "delta", "nu", "lambda", "psi"),
                             responses = NULL, times = NULL, groups = NULL,
                             summary = TRUE, probs = c(0.05, 0.95),
                             include_alpha = TRUE,
                             include_cutpoints = TRUE, ...) {
  stopifnot_(
    !missing(object),
    "Argument {.arg object} is missing."
  )
  stopifnot_(
    is.dynamitefit(object),
    "Argument {.arg object} must be a {.cls dynamitefit} object."
  )
  if (is.null(parameters)) {
    type <- onlyif(is.character(type), tolower(type))
    type <- try(
      match.arg(type, c("beta", "delta", "nu", "lambda", "psi")),
      silent = TRUE
    )
    stopifnot_(
      !inherits(type, "try-error"),
      "Argument {.arg type} must be either {.val beta}, {.val delta},
      {.val nu}, {.val lambda}, or {.val psi}."
    )
    stopifnot_(
      checkmate::test_flag(x = include_alpha),
      "Argument {.arg include_alpha} must be single {.cls logical} value."
    )
    stopifnot_(
      checkmate::test_flag(x = include_cutpoints),
      "Argument {.arg include_cutpoints} must be a single {.cls logical} value."
    )
  }
  if (is.null(parameters) && type %in% c("beta", "delta")) {
    types <- c(
      type,
      onlyif(include_alpha, "alpha"),
      onlyif(include_cutpoints, "cutpoints")
    )
    out <- as.data.frame.dynamitefit(
      object,
      parameters = parameters,
      responses = responses,
      types = types,
      times = times,
      groups = groups,
      summary = summary,
      probs = probs
    )
    # remove extra alphas
    if (identical(type, "delta")) {
      out <- out[!is.na(out$time), ]
    } else {
      out <- out[is.na(out$time), ]
    }
  } else {
    out <- as.data.frame.dynamitefit(
      object,
      parameters = parameters,
      responses = responses,
      types = type,
      times = times,
      groups = groups,
      summary = summary,
      probs = probs
    )
  }
  out
}

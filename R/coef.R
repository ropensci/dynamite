#' Extract Regression Coefficients of a Dynamite Model
#'
#' Extracts either time-varying or time-invariant parameters of the model.
#'
#' @export
#' @family output
#' @param object \[`dynamitefit`]\cr The model fit object.
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
                             types = NULL, parameters = NULL,
                             responses = NULL, times = NULL, groups = NULL,
                             summary = TRUE, probs = c(0.05, 0.95), ...) {
  stopifnot_(
    !missing(object),
    "Argument {.arg object} is missing."
  )
  stopifnot_(
    is.dynamitefit(object),
    "Argument {.arg object} must be a {.cls dynamitefit} object."
  )
  if (is.null(parameters)) {
    types <- onlyif(is.character(types), tolower(types))
    invalid_types <- types[!types %in% all_types]
    stopifnot_(
      length(invalid_types) == 0L,
      c(
        "Argument {.arg types} must contain only valid types.",
        `x` = "Invalid types found: {invalid_types}"
      )
    )
  }
  as.data.frame.dynamitefit(
    x = object,
    types = types,
    parameters = parameters,
    responses = responses,
    times = times,
    groups = groups,
    summary = summary,
    probs = probs
  )
}

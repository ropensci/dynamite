#' Extract Regression Coefficients of Dynamite Model
#'
#' Extracts either time-varying or time-invariant parameters of the model.
#'
#' @export
#' @param object An object of class \code{dynamitefit}.
#' @param type  \[`character(1)`]\cr Either `beta` (the default) for
#'   time-invariant coefficients,  `delta` for time-varying coefficients, or
#'   `nu` for random intercepts.
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), extracts also
#'   time-invariant intercept term alpha if time-invariant parameters beta are
#'   extracted, and time-varying alpha if time-varying delta are extracted.
#' @inheritParams as.data.frame.dynamitefit
#' @param ... Ignored.
#' @examples
#' betas <- coef(gaussian_example_fit, type = "beta")
#' deltas <- coef(gaussian_example_fit, type = "delta")
#'
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @srrstats {RE4.2} *Model coefficients (via `coeff()` / `coefficients()`)*
coef.dynamitefit <- function(object, type = c("beta", "delta", "nu"),
                             summary = TRUE, probs = c(0.05, 0.95),
                             include_alpha = TRUE, ...) {
  type <- match.arg(type)
  if (include_alpha && type != "nu") {
    types <- c("alpha", type)
    out <- as.data.frame(object, types = types, summary = summary,
                         probs = probs)
    # remove extra alphas
    if (type == "delta") {
      out <- out |> dplyr::filter(!is.na(.data$time))
    } else {
      out <- out |> dplyr::filter(is.na(.data$time))
    }
  } else {
    out <- as.data.frame(object, types = type, summary = summary,
                         probs = probs)
  }
  out
}

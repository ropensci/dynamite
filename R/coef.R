#' Extract Regression Coefficients of Dynamite Model
#'
#' Extracts either time-varying or time-invariant parameters of the model.
#'
#' @export
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param type  \[`character(1)`]\cr Either `beta` (the default) for
#'   time-invariant coefficients, `delta` for time-varying coefficients, or
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
#' @srrstats {G2.3a, RE4.2} Provides model coefficients.
coef.dynamitefit <- function(object, type = c("beta", "delta", "nu"),
                             responses = NULL, summary = TRUE,
                             probs = c(0.05, 0.95),
                             include_alpha = TRUE, ...) {
  type <- onlyif(is.character(type), tolower(type))
  type <- try(match.arg(type, c("beta", "delta", "nu")), silent = TRUE)
  stopifnot_(
    !"try-error" %in% class(type),
    "Argument {.arg type} must be either \"beta\", \"delta\", or \"nu\"."
  )
  stopifnot_(
    checkmate::test_flag(x = include_alpha),
    "Argument {.arg include_alpha} must be single {.cls logical} value."
  )
  if (include_alpha && !identical(type, "nu")) {
    types <- c("alpha", type)
    out <- as.data.frame(
      object,
      types = types,
      responses = responses,
      summary = summary,
      probs = probs
    )
    # remove extra alphas
    if (identical(type, "delta")) {
      out <- out |> dplyr::filter(!is.na(.data$time))
    } else {
      out <- out |> dplyr::filter(is.na(.data$time))
    }
  } else {
    out <- as.data.frame(
      object,
      types = type,
      responses = responses,
      summary = summary,
      probs = probs
    )
  }
  out
}

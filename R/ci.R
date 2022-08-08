#' Credible Intervals for Dynamite Model Parameters
#'
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param parm Ignored.
#' @param level \[`numeric(1)`]\cr Credible interval width.
#' @param ... Additional arguments passed to
#'   [dynamite::as.data.frame.dynamitefit()].
#' @return The rows of the resulting matrix will be named using the following
#'   logic: `{parameter}_{time}_{category}_{group}` where `parameter` is the
#'   name of the parameter, `time` is the time index of the parameter,
#'    `category` specifies the level of the response the parameter
#'   is related to if the response is categorical, and `group` determines which
#'   group of observations the parameter is related to in the case of random
#'   intercepts. Non-applicable fields in the this syntax are set to `NA`.
#' @export
#' @examples
#' confint(gaussian_example_fit, level = 0.9)
#'
#' @srrstats {RE4.3} Provides credible intervals.
confint.dynamitefit <- function(object, parm, level = 0.95, ...) {
  stopifnot_(
    !missing(object),
    "Argument {.arg object} is missing."
  )
  stopifnot_(
    is.dynamitefit(object),
    "Argument {.arg object} must be a {.cls dynamitefit} object."
  )
  stopifnot_(
    checkmate::test_number(
      x = level,
      lower = 0.0,
      upper = 1.0,
      na.ok = FALSE,
    ),
    "Argument {.arg level} must be a single
    {.cls numeric} value between 0 and 1."
  )
  a <- (1.0 - level) / 2.0
  d <- as.data.frame.dynamitefit(object, probs = c(a, 1.0 - a))
  row_names <- paste0(
    d$parameter, "_",
    d$time, "_",
    d$category, "_",
    d$group
  )
  out <- d |>
    dplyr::select(
      !c(.data$parameter, .data$time, .data$category,
         .data$group, .data$response, .data$type,
         .data$mean, .data$sd)
    ) |>
    as.matrix()
  colnames(out) <- paste0(100.0 * c(a, 1.0 - a), "%")
  rownames(out) <- row_names
  out
}

#' Credible Intervals for Dynamite Model Parameters
#'
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param parm Ignored.
#' @param level \[`numeric(1)`]\cr Credible interval width.
#' @param ... Additional arguments passed to
#'   [dynamite::as.data.frame.dynamitefit()].
#' @return The rows of the resulting matrix will be named using the following
#'   logic: `{parameter}_{time}_{category}_{group}` where `parameter` is the
#'   name of the parameter, `response` is the name of the response variable
#'   the parameter is related to, `time` is the time index of the parameter,
#'    `category` specifies the level of the response the parameter
#'   is related to if the response is categorical and `group` determined which
#'   group of observations the parameter is related to in the case of random
#'   intercepts. Unrelated fields in the row name syntax are set to `NA`.
#' @export
#' @examples
#' confint(gaussian_example_fit, level = 0.9)
#'
#' @srrstats {RE4.3} *Confidence intervals on those coefficients (via `confint()`)*
confint.dynamitefit <- function(object, parm, level = 0.95, ...) {
  level <- try_type(level, "numeric")[1]
  stopifnot_(
    level > 0.0 && level < 1.0,
    "Argument {.arg level} must be between 0 and 1."
  )
  a <- (1.0 - level)/2.0
  d <- as.data.frame.dynamitefit(object, probs = c(a, 1 - a))
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
  colnames(out) <- paste0(100.0 * c(a, 1 - a), "%")
  rownames(out) <- row_names
  out
}

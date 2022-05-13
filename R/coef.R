#' Extract Regression Coefficients of Dynamite Model
#'
#' @export
#' @param object An object of class \code{dynamitefit}.
#' @param type  \[`character(1)`]\cr Either `beta` (the default) for
#'   time-invariant coefficients or `delta` for time-varying coefficients.
#' @inheritParams as.data.frame.dynamitefit
#' @param ... Ignored.
#' @importFrom stats coef
#' @examples
#' betas <- coef(gaussian_example_fit, type = "beta")
#' deltas <- coef(gaussian_example_fit, type = "delta")
#'
coef.dynamitefit <- function(object, type = c("beta", "delta"),
                             summary = TRUE, probs = c(0.05, 0.95), ...) {
  type <- match.arg(type)
  as.data.frame(object, types = type, summary = summary)
}

#' Visualize Time-varying Regression Coefficients of the Dynamite Model
#'
#' @param model An object of class \code{dynamitefit}.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90\% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity  level for \code{geom_ribbon}.
#'   Default is 0.5.
#' @param scales \[`character(1)`] Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @return A `ggplot` object.
#' @export
#' @examples
#' plot_deltas(gaussian_example_fit)

plot_deltas <- function(model, level = 0.05, alpha = 0.5, scales = "fixed"){

  coef(model, "delta", probs = c(level, 1 - level)) |>
    dplyr::mutate(parameter = gsub("delta_", "", parameter)) |>
    ggplot2::ggplot(aes(.data$time, .data$mean)) +
    ggplot2::geom_ribbon(ggplot2::aes_string(
      ymin = paste0("q", 100 * level),
      ymax = paste0("q", 100 * (1 - level))),
      alpha = alpha) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ parameter, scales = scales) +
    ggplot2::labs(title = "Time-varying coefficients",
         x = "Time", y = "Value")

}

#' Visualize Time-invariant Regression Coefficients of the Dynamite Model
#'
#' @param model An object of class \code{dynamitefit}.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90\% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity  level for \code{geom_ribbon}.
#'   Default is 0.5.
#' @param scales \[`character(1)`] Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @return A `ggplot` object.
#' @export
#' @examples
#' plot_betas(gaussian_example_fit)
plot_betas <- function(model, level = 0.05){
  coef(model, "beta", probs = c(level, 1 - level)) |>
    dplyr::mutate(parameter = gsub("beta_", "", parameter)) |>
    ggplot2::ggplot(aes(.data$mean, .data$parameter)) +
    ggplot2::geom_pointrange(ggplot2::aes_string(
      xmin = paste0("q", 100 * level),
      xmax = paste0("q", 100 * (1 - level)))) +
    ggplot2::labs(title = "Time-invariant coefficients",
         x = "Value", y = "Parameter")
}

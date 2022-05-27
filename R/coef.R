#' Extract Regression Coefficients of Dynamite Model
#'
#' Extracts either time-varying or time-invariant parameters of the model.
#'
#' @export
#' @param object An object of class \code{dynamitefit}.
#' @param type  \[`character(1)`]\cr Either `beta` (the default) for
#'   time-invariant coefficients,  or `delta` for time-varying coefficients.
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), extracts also
#'   time-invariant intercept term alpha if time-invariant parameters beta are
#'   extracted, and time-varying alpha if time-varying delta are extracted.
#' @inheritParams as.data.frame.dynamitefit
#' @param ... Ignored.
#' @examples
#' betas <- coef(gaussian_example_fit, type = "beta")
#' deltas <- coef(gaussian_example_fit, type = "delta")
#'
coef.dynamitefit <- function(object, type = c("beta", "delta"),
                             summary = TRUE, probs = c(0.05, 0.95),
                             include_alpha = TRUE, ...) {
  type <- match.arg(type)
  if (include_alpha) {
    types <- c("alpha", type)
    out <- as.data.frame(object, types = types, summary = summary,
                         probs = probs)
    # remove extra alphas
    if (type == "delta") {
      out <- out |> dplyr::filter(!is.na(time))
    } else {
      out <- out |> dplyr::filter(is.na(time))
    }
  } else {
    out <- as.data.frame(object, types = type, summary = summary,
                         probs = probs)
  }
  out
}

#' Visualize Time-varying Regression Coefficients of the Dynamite Model
#'
#' # TODO Include also alpha it is time-varying?
#' @param model An object of class `dynamitefit`.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity level for `geom_ribbon`.
#'   Default is 0.5.
#' @param scales \[`character(1)`] Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @inheritParams as.data.frame.dynamitefit
#' @return A `ggplot` object.
#' @export
#' @examples
#' plot_deltas(gaussian_example_fit, scales = "free") +
#'   ggplot2::theme_minimal()
plot_deltas <- function(model, level = 0.05, alpha = 0.5, scales = "fixed",
                        include_alpha = TRUE){

  title <- paste0("Posterior mean and ", 100 * (1 - 2 * level),
                  "% intervals of the time-varying coefficients")
  coef(model, "delta", probs = c(level, 1 - level),
       include_alpha = include_alpha) |>
    ggplot2::ggplot(ggplot2::aes(.data$time, .data$mean)) +
    ggplot2::geom_ribbon(ggplot2::aes_string(
      ymin = paste0("q", 100 * level),
      ymax = paste0("q", 100 * (1 - level))),
      alpha = alpha) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ parameter, scales = scales) +
    ggplot2::labs(title = title, x = "Time", y = "Value")

}

#' Visualize Time-invariant Regression Coefficients of the Dynamite Model
#'
#' @inheritParams plot_deltas
#' @return A `ggplot` object.
#' @export
#' @examples
#' plot_betas(gaussian_example_fit)
plot_betas <- function(model, level = 0.05, include_alpha = TRUE){

  title <- paste0("Posterior mean and ", 100 * (1 - 2 * level),
                  "% intervals of the time-invariant coefficients")

  coef(model, "beta", probs = c(level, 1 - level),
       include_alpha = include_alpha) |>
    ggplot2::ggplot(ggplot2::aes(.data$mean, .data$parameter)) +
    ggplot2::geom_pointrange(ggplot2::aes_string(
      xmin = paste0("q", 100 * level),
      xmax = paste0("q", 100 * (1 - level)))) +
    ggplot2::labs(title = title, x = "Value", y = "Parameter")
}

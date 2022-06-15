
#' Visualize Time-varying Regression Coefficients of the Dynamite Model
#'
#' @param model An object of class `dynamitefit`.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity level for `geom_ribbon`.
#'   Default is 0.5.
#' @param scales \[`character(1)`] Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), extracts also
#'   time-invariant intercept term alpha if time-invariant parameters beta are
#'   extracted, and time-varying alpha if time-varying delta are extracted.
#' @return A `ggplot` object.
#' @export
#' @examples
#' plot_deltas(gaussian_example_fit, scales = "free") +
#'   ggplot2::theme_minimal()
#'
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#'
#' @srrstats {BS6.1} *Software should implement a default `plot` method for return objects*
#' @srrstats {BS6.2} *Software should provide and document straightforward abilities to plot sequences of posterior samples, with burn-in periods clearly distinguished*
#' @srrstats {BS6.3} *Software should provide and document straightforward abilities to plot posterior distributional estimates*
#'
#' @srrstats {BS6.5} *Software may provide abilities to plot both sequences of posterior samples and distributional estimates together in single graphic*
#' @srrstats {RE6.0} *Model objects returned by Regression Software (see* **RE4***) should have default `plot` methods, either through explicit implementation, extension of methods for existing model objects, or through ensuring default methods work appropriately.*
#' @srrstats {RE6.1} *Where the default `plot` method is **NOT** a generic `plot` method dispatched on the class of return objects (that is, through an S3-type `plot.<myclass>` function or equivalent), that method dispatch (or equivalent) should nevertheless exist in order to explicitly direct users to the appropriate function.*
#' @srrstats {RE6.2} *The default `plot` method should produce a plot of the `fitted` values of the model, with optional visualisation of confidence intervals or equivalent.*
#' @srrstats {RE6.3} *Where a model object is used to generate a forecast (for example, through a `predict()` method), the default `plot` method should provide clear visual distinction between modelled (interpolated) and forecast (extrapolated) values.*
#' TODO plot_nus
#'
plot_deltas <- function(model, level = 0.05, alpha = 0.5,
                        scales = c("fixed", "free"),
                        include_alpha = TRUE){

  coefs <- coef(model, "delta", probs = c(level, 1 - level),
                include_alpha = include_alpha)
  if (nrow(coefs) == 0) {
    stop_("The model does not contain varying coefficients delta.")
  }

  scales <- match.arg(scales)
  title <- paste0("Posterior mean and ", 100 * (1 - 2 * level),
                  "% intervals of the time-varying coefficients")
  if (any(!is.na(coefs$category))) {
    p <- coefs |>
      ggplot2::ggplot(ggplot2::aes(.data$time, .data$mean,
                                   colour = category, fill = category))
  } else {
    p <- coefs |>
      ggplot2::ggplot(ggplot2::aes(.data$time, .data$mean))
  }
  p + ggplot2::geom_ribbon(ggplot2::aes_string(
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

  coefs <- coef(model, "beta", probs = c(level, 1 - level),
                include_alpha = include_alpha)
  if (nrow(coefs) == 0) {
    stop_("The model does not contain fixed coefficients beta.")
  }
  title <- paste0("Posterior mean and ", 100 * (1 - 2 * level),
                  "% intervals of the time-invariant coefficients")

  if (any(!is.na(coefs$category))) {
    p <- coefs |>
      ggplot2::ggplot(ggplot2::aes(.data$mean, .data$parameter,
                                   colour = category, group = category))
  } else {
    p <- coefs |>
      ggplot2::ggplot(ggplot2::aes(.data$mean, .data$parameter))
  }
  p + ggplot2::geom_pointrange(ggplot2::aes_string(
    xmin = paste0("q", 100 * level),
    xmax = paste0("q", 100 * (1 - level))),
    position = ggplot2::position_dodge(0.5)) +
    ggplot2::labs(title = title, x = "Value", y = "Parameter")
}


#' Visualize Random Intercepts of the Dynamite Model
#'
#' @inheritParams plot_deltas
#' @return A `ggplot` object.
#' @export
#' @examples
#' fit <- dynamite(obs(Reaction ~ lag(Reaction), family = gaussian(),
#'   random_intercept = TRUE), lme4::sleepstudy, "Subject", "Days", chains = 1)
#' nu <-
plot_nus <- function(model, level = 0.05){

  coefs <- coef(model, "nu", probs = c(level, 1 - level))
  if (nrow(coefs) == 0) {
    stop_("The model does not contain random intercepts nu.")
  }
  title <- paste0("Posterior mean and ", 100 * (1 - 2 * level),
                  "% intervals of the random intercepts")

  coefs |>
    ggplot2::ggplot(ggplot2::aes(.data$mean, .data$parameter)) +
    ggplot2::geom_pointrange(ggplot2::aes_string(
      xmin = paste0("q", 100 * level),
      xmax = paste0("q", 100 * (1 - level)))) +
    ggplot2::labs(title = title, x = "Value", y = "Parameter")
}


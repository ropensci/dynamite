#' Plot Regression Coefficients of a Dynamite Model
#'
#' @param x \[`dynamitefit`]\cr The model fit object
#' @param ... Ignored.
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
#' plot(gaussian_example_fit)
#' plot_deltas(gaussian_example_fit, scales = "free") +
#'   ggplot2::theme_minimal()
#' plot_betas(gaussian_example_fit)
#' plot_nus(gaussian_example_fit)
#'
#' @srrstats {G2.3a} Uses match.arg in plot_* functions
#' @srrstats {BS6.1, RE6.0, RE6.1} Implements the `plot` method,
#'   without a default plot
plot.dynamitefit <- function(x, ...) {
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  message_("Please use {.fun plot_deltas}, {.fun plot_betas}, or
           {.fun plot_nus} to produce plots of a {.cls dynamitefit} object.")
}

#' @describeIn plot.dynamitefit
#'   Visualize Time-varying Regression Coefficients of a Dynamite Model
#' @export
plot_deltas <- function(x, level = 0.05, alpha = 0.5,
                        scales = c("fixed", "free"),
                        include_alpha = TRUE) {
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.var x} must be a {.cls dynamitefit} object."
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
  stopifnot_(
    checkmate::test_number(
      x = alpha,
      lower = 0.0,
      upper = 1.0,
      na.ok = FALSE,
    ),
    "Argument {.arg alpha} must be a single
     {.cls numeric} value between 0 and 1."
  )
  coefs <- coef(
    x,
    "delta",
    probs = c(level, 1 - level),
    include_alpha = include_alpha
  )
  stopifnot_(
    nrow(coefs) > 0L,
    "The model does not contain varying coefficients delta."
  )
  scales <- onlyif(is.character(scales), tolower(scales))
  scales <- try(match.arg(scales, c("fixed", "free")), silent = TRUE)
  stopifnot_(
    !"try-error" %in% class(scales),
    "Argument {.arg scales} must be either \"fixed\" or \"free\"."
  )
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the time-varying coefficients"
  )
  if (any(!is.na(coefs$category))) {
    p <- coefs |>
      ggplot2::ggplot(ggplot2::aes(.data$time, .data$mean,
                                   colour = .data$category,
                                   fill = .data$category))
  } else {
    p <- coefs |>
      ggplot2::ggplot(ggplot2::aes(.data$time, .data$mean))
  }
  p + ggplot2::geom_ribbon(ggplot2::aes_string(
    ymin = paste0("q", 100 * level),
    ymax = paste0("q", 100 * (1 - level))),
    alpha = alpha) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ .data$parameter, scales = scales) +
    ggplot2::labs(title = title, x = "Time", y = "Value")
}

#' @describeIn plot.dynamitefit
#'   Time-invariant Regression Coefficients of a Dynamite Model
#' @export
plot_betas <- function(x, level = 0.05, include_alpha = TRUE){
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.var x} must be a {.cls dynamitefit} object."
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
  coefs <- coef(
    x,
    "beta",
    probs = c(level, 1 - level),
    include_alpha = include_alpha
  )
  stopifnot_(
    nrow(coefs) > 0L,
    "The model does not contain fixed coefficients beta."
  )
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the time-invariant coefficients"
  )
  if (any(!is.na(coefs$category))) {
    p <- coefs |>
      ggplot2::ggplot(ggplot2::aes(.data$mean, .data$parameter,
                                   colour = .data$category,
                                   group = .data$category))
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

#' @describeIn plot.dynamitefit
#'   Visualize Random Intercepts of a Dynamite Model
#' @export
plot_nus <- function(x, level = 0.05){
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.var x} must be a {.cls dynamitefit} object."
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
  coefs <- coef(
    x,
    "nu",
    probs = c(level, 1 - level)
  ) |>
    dplyr::mutate(parameter = glue::glue("{parameter}_{group}")) |>
    dplyr::mutate(parameter = factor(.data$parameter,
                                     levels = .data$parameter))
  stopifnot_(
    nrow(coefs) > 0L,
    "The model does not contain random intercepts nu."
  )
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the random intercepts"
  )
  coefs |>
    ggplot2::ggplot(ggplot2::aes(.data$mean, .data$parameter)) +
    ggplot2::geom_pointrange(ggplot2::aes_string(
      xmin = paste0("q", 100 * level),
      xmax = paste0("q", 100 * (1 - level)))) +
    ggplot2::labs(title = title, x = "Value", y = "Parameter")
}


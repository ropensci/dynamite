#' Traceplots and Density Plots of `dynamitefit` Object
#'
#' Plots the traceplots and the density plots of of the model parameters.
#' Possible parameter types are
#'  * `alpha` Intercept terms (time-invariant or time-varying).
#'  * `beta` Time-invariant regression coefficients.
#'  * `delta` Time-varying regression coefficients.
#'  * `nu` Random intercepts.
#'  * `tau` Standard deviations of the spline coefficients of `delta`.
#'  * `tau_alpha` Standard deviations of the spline coefficients of
#'    time-varying `alpha`.
#'  * `sigma_nu` Standard deviation of the random intercepts `nu`.
#'  * `sigma` Standard deviations of gaussian responses.
#'  * `phi` Dispersion parameters of negative binomial responses.
#'  * `omega` Spline coefficients of the regression coefficients `delta`.
#'  * `omega_alpha` Spline coefficients of time-varying `alpha`.
#'
#' Note however typically that drawing these plots for the time-varying
#' parameters `delta` (and `alpha`), spline coefficients, or random
#' intercepts leads to too many plots.
#'
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param responses  \[`character()`]\cr Response(s) for which the plots should
#'   be drawn. Possible options are `unique(x$priors$response)`. Default is
#'   all responses.
#' @param type \[`character(1)`]\cr Type of the parameter for which the plots
#'   should be drawn. See details of possible values.
#' @param ... Further arguments to [bayesplot::mcmc_combo].
#' @return The output object from [bayesplot::mcmc_combo].
#' @export
#' @examples
#' plot(gaussian_example_fit, type = "beta")
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.2, BS6.3, BS6.5} Implements the `plot`
#' method. Further plots can be easily constructed with the help of `as_draws`
#' combined with `ggplot2` and `bayesplot`, for example.
plot.dynamitefit <- function(x, responses = NULL, type, ...) {

  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  stopifnot_(
    !missing(type),
    "Argument {.arg type} is missing while it should be a single
    {.cls character} string."
  )
  stopifnot_(
    checkmate::test_string(x = type, na.ok = FALSE),
    "Argument {.arg type} must be a single {.cls character} string."
  )

  out <- suppressWarnings(as_draws(x, responses = responses, types = type))
  bayesplot::mcmc_combo(out, ...)

}

#' Plot Time-varying Regression Coefficients of a Dynamite Model
#'
#' @param x \[`dynamitefit`]\cr The model fit object
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity level for `geom_ribbon`.
#'   Default is 0.5.
#' @param scales \[`character(1)`] Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), plots also
#'   the time-varying alphas if such parameters exists in the model.
#' @return A `ggplot` object.
#' @examples
#' plot_deltas(gaussian_example_fit, level = 0.025, scales = "free") +
#'   ggplot2::theme_minimal()
#'
#' @srrstats {G2.3a} Uses match.arg.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
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

#' Plot Time-invariant Regression Coefficients of a Dynamite Model
#'
#' @inheritParams plot_deltas
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), plots also
#'   the time-invariant alphas if such parameters exists in the model.
#' @return A `ggplot` object.
#' @examples
#' plot_betas(gaussian_example_fit, level = 0.1)
#'
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
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

#' Plot Random Intercepts of a Dynamite Model
#'
#' @inheritParams plot_deltas
#' @examples
#' plot_nus(gaussian_example_fit)
#'
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
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
  coefs <- try(coef(
    x,
    "nu",
    probs = c(level, 1 - level)
  ), silent = TRUE)
  stopifnot_(
    !"try-error" %in% class(coefs),
    "The model does not contain random intercepts nu."
  )
  coefs <- coefs |>
    dplyr::mutate(parameter = glue::glue("{parameter}_{group}")) |>
    dplyr::mutate(parameter = factor(.data$parameter,
      levels = .data$parameter))

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


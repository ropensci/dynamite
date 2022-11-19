#' Traceplots and Density Plots of a `dynamitefit` Object
#'
#' Produces the traceplots and the density plots of the model parameters.
#' See 'Details' for the available parameter types.
#'
#' Possible parameter types are:
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
#' Note however, that typically drawing these plots for the time-varying
#' parameters `delta` (and `alpha`), spline coefficients, or random
#' intercepts leads to too many plots.
#'
#' @export
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param responses \[`character()`]\cr Response(s) for which the plots should
#'   be drawn. Possible options are `unique(x$priors$response)`. Default is
#'   all responses.
#' @param type \[`character(1)`]\cr Type of the parameter for which the plots
#'   should be drawn. See details of possible values.
#' @param ... Further arguments to [bayesplot::mcmc_combo].
#' @return The output object from [bayesplot::mcmc_combo].
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.2, BS6.3, BS6.5} Implements the `plot`
#' method. Further plots can be easily constructed with the help of `as_draws`
#' combined with `ggplot2` and `bayesplot`, for example.
#' @examples
#' plot(gaussian_example_fit, type = "beta")
#'
plot.dynamitefit <- function(x, responses = NULL, type, ...) {
  stopifnot_(
    !missing(type),
    "Argument {.arg type} is missing while it should be a single
    {.cls character} string."
  )
  stopifnot_(
    checkmate::test_string(x = type, na.ok = FALSE),
    "Argument {.arg type} must be a single {.cls character} string."
  )
  out <- suppressWarnings(
    as_draws_df.dynamitefit(x, responses = responses, types = type)
  )
  bayesplot::mcmc_combo(out, ...)
}

#' Plot Time-varying Regression Coefficients of a Dynamite Model
#'
#' @export
#' @param x \[`dynamitefit`]\cr The model fit object
#' @param responses  \[`character()`]\cr Response(s) for which the coefficients
#'   should be drawn. Possible options are elements of
#'   `unique(x$priors$response)`, and the default is this whole vector.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity level for `geom_ribbon`.
#'   Default is 0.5.
#' @param scales \[`character(1)`] Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), plots also
#'   the time-varying alphas if such parameters exists in the model.
#' @return A `ggplot` object.
#' @srrstats {G2.3a} Uses match.arg.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#' @examples
#' plot_deltas(gaussian_example_fit, level = 0.025, scales = "free") +
#'   ggplot2::theme_minimal()
#'
plot_deltas <- function(x, responses = NULL, level = 0.05, alpha = 0.5,
                        scales = c("fixed", "free"), include_alpha = TRUE) {
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
  coefs <- coef.dynamitefit(
    x,
    "delta",
    responses = responses,
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
    !inherits(scales, "try-error"),
    "Argument {.arg scales} must be either \"fixed\" or \"free\"."
  )
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the time-varying coefficients"
  )
  if (any(!is.na(coefs$category))) {
    p <- ggplot2::ggplot(
      coefs,
      ggplot2::aes_string(
        "time",
        "mean",
        colour = "category",
        fill = "category"
      )
    )
  } else {
    p <- ggplot2::ggplot(coefs, ggplot2::aes_string("time", "mean"))
  }
  p +
    ggplot2::geom_ribbon(
      ggplot2::aes_string(
        ymin = paste0("q", 100 * level),
        ymax = paste0("q", 100 * (1 - level))
      ),
      alpha = alpha
    ) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap("parameter", scales = scales) +
    ggplot2::labs(title = title, x = "Time", y = "Value")
}

#' Plot Time-invariant Regression Coefficients of a Dynamite Model
#'
#' @export
#' @inheritParams plot_deltas
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), plots also
#'   the time-invariant alphas if such parameters exists in the model.
#' @return A `ggplot` object.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#' @examples
#' plot_betas(gaussian_example_fit, level = 0.1)
#'
plot_betas <- function(x, responses = NULL, level = 0.05,
                       include_alpha = TRUE) {
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
  coefs <- coef.dynamitefit(
    x,
    "beta",
    responses = responses,
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
    p <- ggplot2::ggplot(
      coefs,
      ggplot2::aes_string(
        "mean",
        "parameter",
        colour = "category",
        group = "category"
      )
    )
  } else {
    p <- ggplot2::ggplot(coefs, ggplot2::aes_string("mean", "parameter"))
  }
  p +
    ggplot2::geom_pointrange(
      ggplot2::aes_string(
        xmin = paste0("q", 100 * level),
        xmax = paste0("q", 100 * (1 - level))
      ),
      position = ggplot2::position_dodge(0.5)
    ) +
    ggplot2::labs(title = title, x = "Value", y = "Parameter")
}

#' Plot Random Intercepts of a Dynamite Model
#'
#' @export
#' @inheritParams plot_deltas
#' @return A `ggplot` object.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#' @examples
#' plot_nus(gaussian_example_fit)
#'
plot_nus <- function(x, responses = NULL, level = 0.05) {
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
  coefs <- try(
    coef.dynamitefit(
      x,
      "nu",
      responses = responses,
      probs = c(level, 1 - level)
    ),
    silent = TRUE
  )
  stopifnot_(
    !inherits(coefs, "try-error"),
    "The model does not contain random intercepts nu."
  )
  coefs$parameter <- glue::glue("{coefs$parameter}_{coefs$group}")
  coefs$parameter <- factor(coefs$parameter, levels = coefs$parameter)
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the random intercepts"
  )
  ggplot2::ggplot(coefs, ggplot2::aes_string("mean", "parameter")) +
  ggplot2::geom_pointrange(ggplot2::aes_string(
    xmin = paste0("q", 100 * level),
    xmax = paste0("q", 100 * (1 - level))
  )) +
  ggplot2::labs(title = title, x = "Value", y = "Parameter")
}
#' Plot Factor Loadings of a Dynamite Model
#'
#' @export
#' @inheritParams plot_deltas
#' @return A `ggplot` object.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#'
plot_lambdas <- function(x, responses = NULL, level = 0.05) {
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
  coefs <- try(
    coef.dynamitefit(
      x,
      "lambda",
      responses = responses,
      probs = c(level, 1 - level)
    ),
    silent = TRUE
  )
  stopifnot_(
    !inherits(coefs, "try-error"),
    "The model does not contain latent factor psi."
  )
  coefs$parameter <- glue::glue("{coefs$parameter}_{coefs$group}")
  coefs$parameter <- factor(coefs$parameter, levels = coefs$parameter)
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the latent factor loadings"
  )
  ggplot2::ggplot(coefs, ggplot2::aes_string("mean", "parameter")) +
    ggplot2::geom_pointrange(ggplot2::aes_string(
      xmin = paste0("q", 100 * level),
      xmax = paste0("q", 100 * (1 - level))
    )) +
    ggplot2::labs(title = title, x = "Value", y = "Parameter")
}
#' Plot Latent Factors of a Dynamite Model
#'
#' @export
#' @param x \[`dynamitefit`]\cr The model fit object
#' @param responses  \[`character()`]\cr Response(s) for which the coefficients
#'   should be drawn. Possible options are elements of
#'   `unique(x$priors$response)`, and the default is this whole vector.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity level for `geom_ribbon`.
#'   Default is 0.5.
#' @param scales \[`character(1)`] Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @return A `ggplot` object.
#' @srrstats {G2.3a} Uses match.arg.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#'
plot_psis <- function(x, responses = NULL, level = 0.05, alpha = 0.5,
  scales = c("fixed", "free")) {
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
  coefs <- coef.dynamitefit(
    x,
    "psi",
    responses = responses,
    probs = c(level, 1 - level)
  )
  stopifnot_(
    nrow(coefs) > 0L,
    "The model does not contain latent factor psi."
  )
  scales <- onlyif(is.character(scales), tolower(scales))
  scales <- try(match.arg(scales, c("fixed", "free")), silent = TRUE)
  stopifnot_(
    !inherits(scales, "try-error"),
    "Argument {.arg scales} must be either \"fixed\" or \"free\"."
  )
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the latent factors"
  )
  if (any(!is.na(coefs$category))) {
    p <- ggplot2::ggplot(
      coefs,
      ggplot2::aes_string(
        "time",
        "mean",
        colour = "category",
        fill = "category"
      )
    )
  } else {
    p <- ggplot2::ggplot(coefs, ggplot2::aes_string("time", "mean"))
  }
  p +
    ggplot2::geom_ribbon(
      ggplot2::aes_string(
        ymin = paste0("q", 100 * level),
        ymax = paste0("q", 100 * (1 - level))
      ),
      alpha = alpha
    ) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap("parameter", scales = scales) +
    ggplot2::labs(title = title, x = "Time", y = "Value")
}

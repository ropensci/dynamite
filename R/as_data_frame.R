#' Extract Samples From the dynamitefit Object as a Data Frame.
#'
#' You can use the arguments `responses` and `types` to extract only a subset
#' of the model parameters (i.e. only certain types of parameters related to a
#' certain response variable).
#'
#' Potential values for the types argument are
#'  * `alpha` Intercept terms (time-invariant or time-varying).
#'  * `beta` Time-invariant regression coefficients.
#'  * `delta` Time-varying regression coefficients.
#'  * `nu` Random intercepts.
#'  * `tau` Standard deviations of the spline coefficients of `delta`.
#'  * `tau_alpha` Standard deviations of the spline coefficients of
#'    time-varying `alpha`.
#'  * `sigma_nu` Standard deviation of the random intercepts `nu`
#'  * `sigma` Standard deviations of gaussian responses
#'  * `phi` Dispersion parameters of  negative binomial distribution.
#'  * `omega` Spline coefficients of the spline coefficients `delta`.
#'  * `omega_alpha` Spline coefficients of time-varying `alpha`.
#'
#' @param x The estimated \code{dynamite} model.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param responses  \[`character()`]\cr Response(s) for which the  samples
#'   should be extracted. Possible options are `unique(x$priors$response)`
#' @param types \[`character()`]\cr Type(s) of the parameters for which the
#'   samples should be extracted. See details of possible values. Default is
#'   all values listed in details except spline coefficients `omega` and
#'   `omega_alpha`.
#' @param summary \[`logical(1)`]\cr If `TRUE` (default), returns posterior
#'   mean, standard deviation, and posterior quantiles (as defined by the
#'   `probs` argument) for all parameters. If `FALSE`, returns the posterior
#'   samples instead.
#' @param probs \[`numeric()`]\cr Quantiles of interest. Default is
#'   `c(0.05, 0.95)`.
#' @param ... Ignored.
#' @return A `tibble` containing either samples or summary statistics of the
#'   model parameters in a long format. For wide format, see
#'   [dynamite::as_draws()].
#' @export
#' @examples
#' results <- as.data.frame(gaussian_example_fit,
#'   responses = "y", types = "beta", summary = FALSE)
#'
#' results |>
#'   dplyr::group_by(parameter) |>
#'   dplyr::summarise(mean = mean(value), sd = sd(value))
#'
#' # basic summaries can be obtained automatically with summary = TRUE:
#' as.data.frame(gaussian_example_fit,
#'   responses = "y", types = "beta", summary = TRUE)
#'
#' # Compute MCMC diagnostics via posterior package
#' # For this we need to first convert to wide format
#' # and then to draws_df object
#' results |>
#'   dplyr::select(parameter, value, .iteration, .chain) |>
#'   tidyr::pivot_wider(values_from = value, names_from = parameter) |>
#'   posterior::as_draws() |>
#'   posterior::summarise_draws()
#'
#' # Time-varying coefficients delta
#' as.data.frame(gaussian_example_fit,
#'   responses = "y", types = "delta", summary = TRUE)
#'
#' as.data.frame(gaussian_example_fit,
#'   responses = "y", types = "delta", summary = FALSE) |>
#'   dplyr::select(parameter, value, time, .iteration, .chain) |>
#'   tidyr::pivot_wider(
#'     values_from = value,
#'     names_from = c(parameter, time),
#'     names_sep = "_t=") |>
#'   posterior::as_draws() |>
#'   posterior::summarise_draws()
#'
# TODO NA for t <= fixed
as.data.frame.dynamitefit <- function(x, row.names = NULL, optional = FALSE,
                                      responses = NULL, types = NULL,
                                      summary = TRUE, probs = c(0.05, 0.95),
                                      ...) {

  if (is.null(responses)) {
    responses <- unique(x$priors$response)
  } else {
    z <- !(responses %in% unique(x$priors$response))
    if (any(z)) {
      stop_(
        "Model does not contain response variable{?s} {.var {responses[z]}}"
      )
    }
  }
  all_types <- c("alpha", "beta", "delta", "tau", "tau_alpha",
                 "sigma_nu", "sigma", "phi", "nu", "omega", "omega_alpha")

  if (is.null(types)) {
    types <- all_types[1:9]
  } else {
    types <- match.arg(types, all_types, TRUE)
  }

  time_points <- sort(unique(x$data[[x$time_var]]))
  fixed <- x$stan$fixed
  time_points <- time_points[(fixed + 1):length(time_points)]

  values <- function(type, response) {

    draws <- rstan::extract(
      x$stanfit, pars = paste0(type, "_", response), permuted = FALSE)
    n_draws <- prod(dim(draws)[1:2])
    category <- attr(x$stan$responses[[response]], "levels")[-1]
    if (is.null(category)) {
      category <- NA
    }
    if (type == "nu") {
      n_group <- dim(draws)[3]
      d <- data.frame(
        parameter = paste0("nu_", response),
        value = c(draws),
        time = NA,
        category = NA,
        group = rep(1:n_group, each = n_draws),
        .iteration = 1:nrow(draws),
        .chain = rep(1:ncol(draws), each = nrow(draws)))
    }

    if (type == "alpha") {
      if (x$stan$model_vars[[response]]$has_varying_intercept) {
        n_time <- length(time_points)
      } else {
        n_time <- 1
        time_points <- NA
      }
      d <- data.frame(
        parameter = paste0("alpha_", response),
        value = c(draws),
        time = rep(time_points, each = n_draws),
        category = rep(category, each = n_time * n_draws),
        group = NA,
        .iteration = 1:nrow(draws),
        .chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    if (type == "beta") {
      var_names <- paste0("beta_", response, "_",
                          names(x$stan$model_vars[[response]]$J_fixed))
      n_vars <- length(var_names)
      d <- data.frame(
        parameter = rep(var_names, each = n_draws),
        value = c(draws),
        time = NA,
        category = rep(category, each = n_vars * n_draws),
        group = NA,
        .iteration = 1:nrow(draws),
        .chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    if (type == "delta") {
      var_names <- paste0("delta_", response, "_",
                          names(x$stan$model_vars[[response]]$J_varying))
      n_vars <- length(var_names)
      n_time <- length(time_points)
      d <- data.frame(
        parameter = rep(var_names, each = n_time * n_draws),
        value = c(draws),
        time = rep(time_points, each = n_draws),
        category = rep(category, each = n_time * n_vars * n_draws),
        group = NA,
        .iteration = 1:nrow(draws),
        .chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    if (type == "tau") {
      var_names <- paste0("tau_", response, "_",
                          names(x$stan$model_vars[[response]]$J_varying))
      d <- data.frame(
        parameter = rep(var_names, each = n_draws),
        value = c(draws),
        time = NA,
        category = NA,
        group = NA,
        .iteration = 1:nrow(draws),
        .chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    if (type %in% c("tau_alpha", "sigma", "phi", "sigma_nu")) {
      d <- data.frame(
        parameter = paste0(type, "_", response),
        value = c(draws),
        time = NA,
        category = NA,
        group = NA,
        .iteration = 1:nrow(draws),
        .chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    if (type == "omega") {
      D <- x$stan$sampling_vars$D
      var_names <- names(x$stan$model_vars[[response]]$J_varying)
      k <- length(var_names)
      d <- data.frame(
        parameter = rep(
          paste0("omega_", rep(1:D, each = k), "_", var_names), each = n_draws),
        value = c(draws),
        time = NA,
        category = NA,
        group = NA,
        .iteration = 1:nrow(draws),
        .chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    if (type == "omega_alpha") {
      D <- x$stan$sampling_vars$D
      d <- data.frame(
        parameter = rep(paste0("omega_alpha_", 1:D), each = n_draws),
        value = c(draws),
        time = NA,
        category = NA,
        group = NA,
        .iteration = 1:nrow(draws),
        .chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    d
  }
  out <- tidyr::expand_grid(type = types, response = responses) |>
    dplyr::mutate(parameter = glue::glue("{type}_{response}")) |>
    dplyr::rowwise() |>
    dplyr::filter(any(grepl(paste0("^", .data$parameter),
                            x$stanfit@sim$pars_oi))) |>
    dplyr::select(.data$response, .data$type) |>
    dplyr::mutate(value = list(values(.data$type, .data$response))) |>
    tidyr::unnest(cols = .data$value)
  if (summary) {
    out <- out |>
      dplyr::group_by(.data$parameter, .data$time, .data$category,
                      .data$response, .data$type) |>
      dplyr::summarise(
        mean = mean(.data$value),
        sd = sd(.data$value),
        # use quantile2 from posterior for simpler (more R-friendly) names
        dplyr::as_tibble(
          as.list(posterior::quantile2(.data$value, probs = probs)))
        ) |>
      dplyr::ungroup()
  }
  out
}

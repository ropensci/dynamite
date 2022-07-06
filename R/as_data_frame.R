#' Extract Samples From the `dynamitefit` Object as a Data Frame.
#'
#' You can use the arguments `responses` and `types` to extract only a subset
#' of the model parameters (i.e., only certain types of parameters related to a
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
#'  * `sigma_nu` Standard deviation of the random intercepts `nu`.
#'  * `sigma` Standard deviations of gaussian responses.
#'  * `phi` Dispersion parameters of negative binomial responses.
#'  * `omega` Spline coefficients of the regression coefficients `delta`.
#'  * `omega_alpha` Spline coefficients of time-varying `alpha`.
#'
#' @param x  \[`dynamitefit`]\cr The model fit object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param responses  \[`character()`]\cr Response(s) for which the samples
#'   should be extracted. Possible options are elements of
#'   `unique(x$priors$response)`, and the default is this whole vector.
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
#' @param include_fixed \[`logical(1)`]\cr If `TRUE` (default), time-varying
#'   parameters for `1:fixed` time points are included in the output as `NA`
#'   values. If `FALSE`, fixed time points are omitted completely
#'   from the output.
#' @param ... Ignored.
#' @return A `tibble` containing either samples or summary statistics of the
#'   model parameters in a long format. For a wide format, see
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
as.data.frame.dynamitefit <- function(x, row.names = NULL, optional = FALSE,
                                      responses = NULL, types = NULL,
                                      summary = TRUE, probs = c(0.05, 0.95),
                                      include_fixed = TRUE, ...) {
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  stopifnot_(
    checkmate::test_character(
      x = responses,
      any.missing = FALSE,
      min.len = 1L,
      null.ok = TRUE
    ),
    "Argument {.arg responses} must be a {.cls character} vector."
  )
  stopifnot_(
    checkmate::test_character(
      x = types,
      any.missing = FALSE,
      min.len = 1L,
      null.ok = TRUE
    ),
    "Argument {.arg types} must be a {.cls character} vector."
  )
  stopifnot_(
    checkmate::test_flag(x = summary),
    "Argument {.arg summary} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_numeric(
      x = probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values between
     0 and 1."
  )
  stopifnot_(
    checkmate::test_flag(x = include_fixed),
    "Argument {.arg include_fixed} must be a single {.cls logical} value."
  )
  if (is.null(responses)) {
    responses <- unique(x$priors$response)
  } else {
    z <- responses %in% unique(x$priors$response)
    stopifnot_(
      all(z),
      "Model does not contain response variable{?s} {.var {responses[!z]}}."
    )
  }
  all_types <- c(
    "alpha", "beta", "delta", "tau", "tau_alpha",
    "sigma_nu", "sigma", "phi", "nu", "omega", "omega_alpha"
  )
  if (is.null(types)) {
    types <- all_types[1L:9L]
  } else {
    types <- onlyif(is.character(types), tolower(types))
    types <- try(match.arg(types, all_types, TRUE), silent = TRUE)
    stopifnot_(
      !"try-error" %in% class(types),
      "Argument {.arg type} contains unknown types."
    )
  }
  all_time_points <- sort(unique(x$data[[x$time_var]]))
  fixed <- x$stan$fixed
  values <- function(type, response) {
    draws <- rstan::extract(
      x$stanfit,
      pars = paste0(type, "_", response),
      permuted = FALSE
    )
    n_draws <- prod(dim(draws)[1L:2L])
    category <- attr(x$stan$responses[[response]], "levels")[-1L]
    if (is.null(category)) {
      category <- NA
    }
    d <- switch(type,
      `nu` = {
        n_group <- dim(draws)[3L]
        data.frame(
          parameter = paste0("nu_", response),
          value = c(draws),
          time = NA,
          category = NA,
          group = rep(seq_len(n_group), each = n_draws)
        )
      },
      `alpha` = {
        if (x$stan$model_vars[[response]]$has_varying_intercept) {
          time_points <- ifelse_(
            include_fixed,
            all_time_points,
            all_time_points[seq.int(fixed + 1L, length(all_time_points))]
          )
          n_time <- length(time_points)
          n_na <- include_fixed * fixed * n_draws
        } else {
          n_time <- 1L
          time_points <- NA
          n_na <- 0L
        }
        data.frame(
          parameter = paste0("alpha_", response),
          value = c(rep(NA, n_na), c(draws)),
          time = rep(time_points, each = n_draws),
          category = rep(category, each = n_time * n_draws),
          group = NA
        )
      },
      `beta` = {
        var_names <- paste0(
          "beta_", response, "_",
          names(x$stan$model_vars[[response]]$J_fixed)
        )
        n_vars <- length(var_names)
        data.frame(
          parameter = rep(var_names, each = n_draws),
          value = c(draws),
          time = NA,
          category = rep(category, each = n_vars * n_draws),
          group = NA
        )
      },
      `delta` = {
        var_names <- paste0(
          "delta_", response, "_",
          names(x$stan$model_vars[[response]]$J_varying)
        )
        n_vars <- length(var_names)
        time_points <- ifelse_(
          include_fixed,
          all_time_points,
          all_time_points[seq.int(fixed + 1L, length(all_time_points))]
        )
        n_na <- include_fixed * fixed * n_draws
        n_time <- length(time_points)
        n_time2 <- n_time - include_fixed * fixed
        d <- do.call(dplyr::bind_rows, lapply(seq_len(n_vars), function(i) {
          data.frame(
            parameter = var_names[i],
            value = c(rep(NA, n_na),
                      c(draws[, , (i - 1L) * n_time2 + seq_len(n_time2)])),
            time = rep(time_points, each = n_draws),
            category = rep(category, each = n_time * n_draws),
            group = NA
          )
        }))
      },
      `tau` = {
        var_names <- paste0(
          "tau_", response, "_",
          names(x$stan$model_vars[[response]]$J_varying)
        )
        d <- data.frame(
          parameter = rep(var_names, each = n_draws),
          value = c(draws),
          time = NA,
          category = NA,
          group = NA
        )
      },
      `omega` = {
        D <- x$stan$sampling_vars$D
        var_names <- names(x$stan$model_vars[[response]]$J_varying)
        k <- length(var_names)
        data.frame(
          parameter = rep(
            paste0("omega_",
                   rep(seq_len(D), each = k), "_", var_names), each = n_draws),
          value = c(draws),
          time = NA,
          category = NA,
          group = NA
        )
      },
      `omega_alpha` = {
        D <- x$stan$sampling_vars$D
        d <- data.frame(
          parameter = rep(paste0("omega_alpha_", seq_len(D)), each = n_draws),
          value = c(draws),
          time = NA,
          category = NA,
          group = NA
        )
      },
      { # default case for tau_alpha, sigma, phi and sigma_nu
        data.frame(
          parameter = paste0(type, "_", response),
          value = c(draws),
          time = NA,
          category = NA,
          group = NA
        )
      }
    )
    d$.iteration <- seq_len(nrow(draws))
    d$.chain <- rep(seq_len(ncol(draws)), each = nrow(draws))
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
      dplyr::group_by(.data$parameter, .data$time, .data$category, .data$group,
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

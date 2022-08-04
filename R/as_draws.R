#' Convert `dynamite` Output to `draws_df` Format
#'
#' Converts the output from [dynamite::dynamite()] call to a
#' `draws_df` format of the `posterior` package, enabling the use
#' of diagnostics and plotting methods of `posterior` and `bayesplot`
#' packages. Note that this function returns all variables in a wide format,
#' whereas [dynamite::as.data.frame()] uses the long format.
#'
#' You can use the arguments `responses` and `types` to extract only a subset
#' of the model parameters (i.e., only certain types of parameters related to a
#' certain response variable).
#'
#' Potential values for the types argument are:
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
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @inheritParams as.data.frame.dynamitefit
#' @return A `draws_df` object.
#' @aliases as_draws as_draws_df
#' @export
#' @export as_draws_df
#' @rdname as_draws-dynamitefit
#' @method as_draws_df dynamitefit
#' @examples
#' as_draws(gaussian_example_fit, types = c("sigma", "beta"))
#'
as_draws_df.dynamitefit <- function(x, responses = NULL, types = NULL, ...) {
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  d <- as.data.frame.dynamitefit(
    x,
    responses = responses,
    types = types,
    summary = FALSE,
    include_fixed = FALSE
  ) |>
    dplyr::select(.data$parameter, .data$value, .data$time, .data$category,
                  .data$group, .data$.iteration, .data$.chain) |>
    dplyr::arrange(.data$parameter, .data$time, .data$category, .data$group,
      .data$.chain, .data$.iteration) |>
    tidyr::pivot_wider(
      values_from = .data$value,
      names_from = c(.data$parameter, .data$time, .data$category,
                     .data$group),
      names_glue = "{parameter}[{time}]_{category}_id{group}"
    )
  # remove NAs from time-invariant parameter names
  colnames(d) <- gsub("\\[NA\\]", "", colnames(d))
  # remove NAs from parameters which are not category specific
  colnames(d) <- gsub("_NA", "", colnames(d))
  # remove NAs from parameters which are not group specific
  colnames(d) <- gsub("_idNA", "", colnames(d))
  d |> posterior::as_draws()
}

#' @export
#' @export as_draws
#' @rdname as_draws-dynamitefit
#' @method as_draws dynamitefit
#' @inheritParams as_draws_df.dynamitefit
as_draws.dynamitefit <- function(x, responses = NULL, types = NULL, ...) {
  as_draws_df.dynamitefit(x, responses, types, ...)
}

#' Extract Samples From the dynamitefit Object as a Data Frame.
#'
#' You can use the arguments `responses` and `types` to extract only a subset
#' of the model parameters (i.e. only certain types of parameters related to a
#' certain response variable). The prior data frame (from `get_priors`) shows
#' potential values for these variables.
#'
#' @note The spline coefficients `alpha` are never returned, but they can be
#' obtained with `as_draws` function among all the other parameters.
#'
#' @param x The estimated \code{dynamite} model.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param responses  \[`character()`]\cr Response(s) for which the  samples
#' should be extracted.
#' @param types \[`character()`]\cr Type(s) of the parameters for which the
#' samples should be extracted.
#' @param ... Ignored.
#' @importFrom tidyr unnest
#' @export
as.data.frame.dynamitefit <- function(x, row.names = NULL, optional = FALSE,
                                      responses = NULL, types = NULL, ...) {

  if (is.null(responses)) {
    responses <- unique(x$priors$response)
  }
  if (is.null(types)) {
    types <- unique(x$priors$type)
  }


  values <- function(type, response) {

    draws <- rstan::extract(
      x$stanfit, pars = paste0(type, "_", response), permuted = FALSE)
    n_draws <- prod(dim(draws)[1:2])
    category <- attr(x$responses[[response]], "levels")[-1]
    if(is.null(category)) category <- NA

    if (type == "beta") {
      var_names <- paste0("beta_", response, "_",
                          names(x$model_vars[[response]]$L_fixed))
      n_vars <- length(var_names)
      d <- data.frame(
        parameter = rep(var_names, each = n_draws),
        value = c(draws),
        time = NA,
        category = rep(category, each = n_vars * n_draws),
        iter = 1:nrow(draws),
        chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    if (type == "delta") {
      var_names <- paste0("delta_", response, "_",
                          names(x$model_vars[[response]]$L_varying))
      n_vars <- length(var_names)
      n_time <- length(x$time)
      d <- data.frame(
        parameter = rep(var_names, each = n_time * n_draws),
        value = c(draws),
        time = rep(x$time, each = n_draws),
        category = rep(category, each = n_time * n_vars * n_draws),
        iter = 1:nrow(draws),
        chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    if (type == "tau") {
      var_names <- paste0("tau_", response, "_",
                          names(x$model_vars[[response]]$L_varying))
      d <- data.frame(
        parameter = rep(var_names, each = n_draws),
        value = c(draws),
        time = NA,
        category = NA,
        iter = 1:nrow(draws),
        chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    if (type %in% c("sigma", "phi")) {
      d <- data.frame(
        parameter = paste0(type, "_", response),
        value = c(draws),
        time = NA,
        category = NA,
        iter = 1:nrow(draws),
        chain = rep(1:ncol(draws), each = nrow(draws)))
    }
    d
  }
  x$priors |>
    dplyr::select(.data$response, .data$type) |>
    dplyr::filter(.data$response %in% responses & .data$type %in% types) |>
    dplyr::distinct() |>
    dplyr::rowwise() |>
    dplyr::mutate(value = list(values(.data$type, .data$response))) |>
    tidyr::unnest(cols = .data$value)

}

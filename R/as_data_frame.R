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
#' @importFrom stats quantile
#' @param x The estimated \code{dynamite} model.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param responses  \[`character()`]\cr Response(s) for which the  samples
#' should be extracted. Possible options are `unique(x$priors$response)`
#' @param types \[`character()`]\cr Type(s) of the parameters for which the
#' samples should be extracted. Possible options are `unique(x$priors$type)`.
#' @param summary \[`logical(1)`]\cr If `TRUE` (default), returns posterior
#'   mean and lower and upper limits of the 95% posterior intervals for all
#'   parameters. If `FALSE`, returns all the posterior samples instead.
#' @param ... Ignored.
#' @importFrom tidyr unnest
#' @export
#' @examples
#' results <- as.data.frame(gaussian_example_fit,
#'   responses = "y", types = "beta")
#' results |>
#'   dplyr::group_by(parameter) |>
#'   dplyr::summarise(mean = mean(value), sd = sd(value))
#'
as.data.frame.dynamitefit <- function(x, row.names = NULL, optional = FALSE,
                                      responses = NULL, types = NULL,
                                      summary = TRUE, ...) {

  if (is.null(responses)) {
    responses <- unique(x$priors$response)
  } else {
    z <- !(responses %in% unique(x$priors$response))
    if (sum(!z) == 0) {
      stop("Model does not contain any of the input responses. ")
    }
    if (any(z)) {
      warning("Model does not contain response(s) ",
              paste0(responses[z], collapse = ", "))
    }

  }
  if (is.null(types)) {
    types <- unique(x$priors$type)
  } else {
    z <- !(types %in% unique(x$priors$type))
    if (sum(!z) == 0) {
      stop("Model does not contain any variables of chosen type(s). ")
    }
    if (any(z)) {
      warning("Model does not contain variable type(s) ",
              paste0(types[z], collapse = ", "))
    }

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
  out <- x$priors |>
    dplyr::select(.data$response, .data$type) |>
    dplyr::filter(.data$response %in% responses & .data$type %in% types) |>
    dplyr::distinct() |>
    dplyr::rowwise() |>
    dplyr::mutate(value = list(values(.data$type, .data$response))) |>
    tidyr::unnest(cols = .data$value)
  if (summary) {
    out <- out |>
      dplyr::group_by(.data$parameter, .data$time,
                      .data$response, .data$type) |>
      dplyr::summarise(mean = mean(.data$value),
                `2.5%` = quantile(.data$value, 0.025),
                `97.5%` = quantile(.data$value, 0.975)) |>
      dplyr::ungroup()
  }
  out
}

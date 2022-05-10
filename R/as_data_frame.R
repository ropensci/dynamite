#' Extract Samples From the dynamitefit Object as a Data Frame.
#'
#' @param x The estimated \code{dynamite} model.
#' @param row.names \code{NULL} (default) or a character vector giving the row
#' names for the data frame.
#' @param optional Ignored.
#' @param responses TODO
#' @param types TODO
#' @param ... Ignored.
#' @export
as.data.frame.dynamitefit <- function(x, row.names = NULL,
                                      optional = FALSE, responses = NULL,
  types = NULL, ...) {

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
  x$priors %>%
    dplyr::select(response, type) %>%
    dplyr::filter(response %in% responses & type %in% types) %>%
    dplyr::distinct() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(value = list(values(type, response))) %>%
    tidyr::unnest(cols = value)

}

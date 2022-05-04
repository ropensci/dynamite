#' Extract samples from the dynamitefit object as data frame.
#'
#' @param x The estimated \code{dynamite} model.
#' @param row.names \code{NULL} (default) or a character vector giving the row
#' names for the data frame.
#' @param optional Ignored.
#' @param parameter_types What type of parameters should be returned? Possible
#' choices are \code{"beta"}, \code{"tau"}, \code{"a"}, \code{"beta"}, or
#' \code{"all"} (default). #TODO sigmas etc? others? Or just return all of these...
#' @param ... Ignored.
#' @export
as.data.frame.dynamitefit <- function(x, row.names = NULL,
                                      optional = FALSE, parameter_types, ...) {
  # should there be an option to be more specific in what parameters to extract?
  # TODO distributional parameters (sigma, shape)
  parameter_types <- match.arg(
    parameter_types,
    c("beta", "tau", "a", "lambda", "all"), TRUE
  ) # TODO all?, naming? a = spline_coefficients? alpha (greek as others)?
  if (parameter_types == "all") {
    parameter_types <- c("beta", "tau", "a", "lambda")
  }
  all_pars <- x$stanfit@sim$pars_oi
  d <- NULL
  coef_names <- unlist(x$coef_names)
  if ("tau" %in% parameter_types) {
    idx <- grep("^tau", all_pars)
    samples <- rstan::extract(x$stanfit, pars = all_pars[idx], permuted = FALSE)
    var_names <- paste0("tau_", coef_names)
    n <- nrow(samples) # samples per chain
    k <- ncol(samples) # number of chains
    d <- data.frame(
      iter = 1:n,
      chain = rep(1:k, each = n),
      time = NA,
      value = c(samples),
      variable = rep(var_names, each = n * k),
      row.names = row.names
    )
  }
  if ("lambda" %in% parameter_types) {
    samples <- rstan::extract(x$stanfit, pars = "lambda", permuted = FALSE)
    n <- nrow(samples) # samples per chain
    k <- ncol(samples) # number of chains
    d <- rbind(d, data.frame(
      iter = 1:n,
      chain = rep(1:k, each = n),
      time = NA,
      value = c(samples),
      variable = "lambda",
      row.names = row.names
    ))
  }
  if ("a" %in% parameter_types) { # TODO name
    idx <- grep("^a_(?!r)", all_pars, perl = TRUE)
    samples <- rstan::extract(x$stanfit, pars = all_pars[idx], permuted = FALSE)
    # TODO categorical case, how to name S-1 cases?
    # a is K x D or S x K x D depending on channel
    a <- rep(paste0("a[", 1:x$spline$D, "]"), each = length(coef_names))
    var_names <- paste0(a, "_", coef_names)
    n <- nrow(samples) # samples per chain
    k <- ncol(samples) # number of chains
    d <- rbind(d, data.frame(
      iter = 1:n,
      chain = rep(1:k, each = n),
      time = NA,
      value = c(samples),
      variable = rep(var_names, each = n * k),
      row.names = row.names
    ))
  }
  if ("beta" %in% parameter_types) { # TODO name
    idx <- grep("^beta", all_pars, perl = TRUE)
    samples <- rstan::extract(x$stanfit, pars = all_pars[idx], permuted = FALSE)
    # beta is T x K or T x K x S
    fixed <- x$prediction_basis$fixed
    time <- if (fixed > 0) x$time[-(1:fixed)] else x$time
    var_names <- rep(paste0("beta_", coef_names),
      each = length(time)
    )
    n <- nrow(samples) # samples per chain
    k <- ncol(samples) # number of chains
    d <- rbind(d, data.frame(
      iter = 1:n,
      chain = rep(1:k, each = n),
      time = rep(time, each = n * k), # times number of coefficients
      value = c(samples),
      variable = rep(var_names, each = n * k),
      row.names = row.names
    ))
  }
  d
}

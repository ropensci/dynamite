#' Extract samples from the btvcmfit object as data frame.
#'
#' @export
as.data.frame.btvcmfit <- function(x, row.names = NULL, optional = FALSE, parameter_types, ...) {
    # should there be an option to be more specific in what parameters to extract?
    # TODO distributional parameters (sigma, shape)
    parameter_types <- match.arg(parameter_types,
        c("beta", "tau", "a", "lambda"), TRUE) #TODO a = spline_coefficients? alpha (greek as others)?
    all_pars <- x$stanfit@sim$pars_oi
    d <- NULL
    coef_names <- unlist(x$coef_names)
    if ("tau" %in% parameter_types) {
        idx <-  grep("^tau", all_pars)
        samples <- rstan::extract(x$stanfit, pars = all_pars[idx], permuted = FALSE)
        # TODO variable names should be something like tau_response_coefficient?
        var_names <- paste0("tau_", coef_names)
        n <- nrow(samples) # samples per chain
        k <- ncol(samples) # number of chains
        d <- data.frame(
            iter = 1:n,
            chain = rep(1:k, each = n),
            time = NA,
            value = c(samples),
            variable = rep(var_names, each = n * k),
            row.names = row.names)
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
            row.names = row.names))
    }
    if ("a" %in% parameter_types) { #TODO name
        idx <- grep("^a_(?!r)", all_pars, perl = TRUE)
        samples <- rstan::extract(x$stanfit, pars = all_pars[idx], permuted = FALSE)
        # TODO categorical case, how to name S-1 cases?
        # a is K x D or S x K x D depending on channel
        a <- rep(paste0("a[", 1:x$model_data$D, "]"), each = length(coef_names))
        var_names <- paste0(a, "_", coef_names)
        n <- nrow(samples) # samples per chain
        k <- ncol(samples) # number of chains
        d <- rbind(d, data.frame(
            iter = 1:n,
            chain = rep(1:k, each = n),
            time = NA,
            value = c(samples),
            variable = rep(var_names, each = n * k),
            row.names = row.names))
    }
    if ("beta" %in% parameter_types) { #TODO name
        idx <- grep("^beta", all_pars, perl = TRUE)
        samples <- rstan::extract(x$stanfit, pars = all_pars[idx], permuted = FALSE)
        # TODO categorical case, how to name S-1 cases (or S if we don't drop the zeros)?
        # beta is T x K or T x K x S
        var_names <- rep(paste0("beta_", coef_names), each = x$model_data$T)
        n <- nrow(samples) # samples per chain
        k <- ncol(samples) # number of chains
        d <- rbind(d, data.frame(
            iter = 1:n,
            chain = rep(1:k, each = n),
            time = rep(x$time, each = n * k), # times number of coefficients
            value = c(samples),
            variable = rep(var_names, each = n * k),
            row.names = row.names))
    }
    d
}

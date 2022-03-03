#' Extract fitted values of btvcmfit
#'
#'
#' @export
#' @param object An object of class \code{btvcmfit}.
#' @param ... Ignored.
#' @importFrom stats fitted
fitted.btvcmfit <- function(object, newdata = NULL,...) {

    # TODO By default use original data, stored where and in which format?
    # see bssm::fitted.mcmc_output, brms::fitted.brmsfit etc
    # TODO need to check that newdata is compatible with the model fit

    samples <- rstan::extract(object$stanfit)
    # TODO compute expected values of the posterior predictive distribution

}

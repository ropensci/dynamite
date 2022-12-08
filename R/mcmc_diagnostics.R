#' Diagnostic Values of a Dynamite Model
#'
#' Prints HMC diagnostics, and lists parameters with smallest effective sample
#' sizes and largest Rhat values. See [rstan::check_hmc_diagnostics()] and
#' [posterior::default_convergence_measures()] for details.
#'
#' @export
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param n \[`integer(1)`]\cr How many rows to print in
#'   parameter-specific convergence measures. The default is 1. Should be a
#'   positive (unrestricted) integer.
#' @return Returns `x` (invisibly).
#' @examples
#' mcmc_diagnostics(gaussian_example_fit)
#'
mcmc_diagnostics <- function(x, n) {
  UseMethod("mcmc_diagnostics", x)
}

#' @export
#' @rdname mcmc_diagnostics
mcmc_diagnostics.dynamitefit <- function(x, n = 1L) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  stopifnot_(
    checkmate::test_int(
      x = n,
      lower = 1L
    ),
    "Argument {.arg n} must be a single {.cls integer}."
  )
  if (!is.null(x$stanfit)) {
    if (identical(x$stanfit@stan_args[[1L]]$algorithm, "NUTS")) {
      cat("NUTS sampler diagnostics:\n")
      invisible(utils::capture.output(msg <-
        utils::capture.output(rstan::check_hmc_diagnostics(x$stanfit),
          type = "message"
        )))
      cat(msg, sep = "\n")
    }
    init <- seq_len(n)
    sumr <- posterior::summarise_draws(
      suppressWarnings(as_draws(x)),
      posterior::default_convergence_measures()
    )
    cat("\nSmallest bulk-ESS values: \n")
    bulk <- sumr[order(sumr$ess_bulk), c("variable", "ess_bulk")][init, ]
    var <- bulk$variable
    out <- bulk$ess_bulk
    names(out) <- var
    print(tibble::as_tibble(t(out)))
    cat("\nSmallest tail-ESS values: \n")
    tail <- sumr[order(sumr$ess_tail), c("variable", "ess_tail")][init, ]
    var <- tail$variable
    out <- tail$ess_tail
    names(out) <- var
    print(tibble::as_tibble(t(out)))
    cat("\nLargest Rhat values: \n")
    rhat <- sumr[
      order(sumr$rhat, decreasing = TRUE),
      c("variable", "rhat")
    ][init, ]
    var <- rhat$variable
    out <- rhat$rhat
    names(out) <- var
    print(tibble::as_tibble(t(out)))
  } else {
    cat("No Stan model fit is available.")
  }
  invisible(x)
}

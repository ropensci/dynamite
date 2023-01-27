#' Diagnostic Values of a Dynamite Model
#'
#' Prints HMC diagnostics, and lists parameters with smallest effective sample
#' sizes and largest Rhat values. See [rstan::check_hmc_diagnostics()] and
#' [posterior::default_convergence_measures()] for details.
#'
#' @export
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param n \[`integer(1)`]\cr How many rows to print in
#'   parameter-specific convergence measures. The default is 3. Should be a
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
mcmc_diagnostics.dynamitefit <- function(x, n = 3L) {
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
      n_draws <- ndraws(x)
      n_divs <- rstan::get_num_divergent(x$stanfit)
      n_trees <- rstan::get_num_max_treedepth(x$stanfit)
      bfmis <- rstan::get_bfmi(x$stanfit)
      all_ok <- n_divs == 0 && n_trees == 0 && all(bfmis > 0.2)
      cat("\nNUTS sampler diagnostics:\n")
      if (all_ok) {
        cat("\nNo divergences, saturated max treedepths or low E-BFMIs.\n")
      }
      if (n_divs > 0) {
        cat("\n", n_divs, " out of ", n_draws, " iterations ended with a ",
          "divergence. See Stan documentation for details.\n", sep = "")
      }
      if (n_trees > 0) {
        mt <- x$stanfit@stan_args[[1]]$control$max_treedepth
        mt <- ifelse_(is.null(mt), 10, mt)
        cat("\n", n_trees, " out of ", n_draws, " saturated the maximum ",
          "tree depth of ", mt, ". See Stan documentation for details.\n",
          sep = "")
      }
      if (any(bfmis < 0.2)) {
        cat("\nChain(s)", cs(which(bfmis < 0.2)), "had E-BFMI below 0.2,",
          "indicating possible issues. See Stan documentation for details.\n")
      }
    }
    init <- seq_len(n)
    sumr <- posterior::summarise_draws(
      suppressWarnings(as_draws(x)),
      posterior::default_convergence_measures()
    )
    cat("\nSmallest bulk-ESS values: \n")
    bulk <- sumr[order(sumr$ess_bulk), c("variable", "ess_bulk")][init, ]
    out <- matrix(bulk$ess_bulk, dimnames = list(bulk$variable, ""))
    print(out, digits = 1)
    cat("\nSmallest tail-ESS values: \n")
    tail <- sumr[order(sumr$ess_tail), c("variable", "ess_tail")][init, ]
    out <- matrix(tail$ess_tail, dimnames = list(tail$variable, ""))
    print(out, digits = 1)
    cat("\nLargest Rhat values: \n")
    rhat <- sumr[
      order(sumr$rhat, decreasing = TRUE),
      c("variable", "rhat")
    ][init, ]
    out <- matrix(rhat$rhat, dimnames = list(rhat$variable, ""))
    print(out, digits = 3)
  } else {
    cat("No Stan model fit is available.")
  }
  invisible(x)
}

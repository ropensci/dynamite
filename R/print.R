#' Print Summary of Dynamite Object
#'
#' Prints the summary information of the estimated dynamite model.
#'
#' @param x An output from \code{\link{dynamite}}.
#' @param hmc_diagnostics \[`logical(1)`]\cr if `TRUE` (default), prints the
#' summary of Stan's HMC sampler diagnostics. See [rstan::check_hmc_diagnostics()]
#' for details.
#' @method print dynamitefit
#' @export
print.dynamitefit <- function(x, hmc_diagnostics = TRUE,...) {
  if (!is.null(x$stanfit)) {
    if (hmc_diagnostics) {
      rstan::check_hmc_diagnostics(x$stanfit)
    }
    draws <- suppressWarnings(as_draws(x))
    sumr <- posterior::summarise_draws(draws,
              posterior::default_convergence_measures())
    min_ess <- which.min(sumr$ess_bulk)
    cat("\nSmallest bulk-ESS: ", round(sumr$ess_bulk[min_ess]), " (",
        sumr$variable[min_ess], ")", sep = "")
    min_ess <- which.min(sumr$ess_tail)
    cat("\nSmallest tail-ESS: ", round(sumr$ess_tail[min_ess]), " (",
        sumr$variable[min_ess], ")", sep = "")
    max_rhat <- which.max(sumr$rhat)
    cat("\nLargest Rhat: ", round(sumr$rhat[max_rhat], 3), " (",
        sumr$variable[max_rhat], ")", sep = "")
    cat("\n\nElapsed time (seconds):\n")
    print(rstan::get_elapsed_time(x$stanfit))

    cat("\nSummary statistics of the time-invariant parameters:\n")
    print(draws |> dplyr::select(matches("([^\\]])$", perl = TRUE)) |>
      posterior::summarise_draws())

  } else {
    message("No Stan model fit is available.")
  }
  invisible(x)
}


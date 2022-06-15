#' Diagnostic Values of Dynamite Model
#'
#' Prints HMC diagnostics, and lists parameters with smallest effective sample
#' sizes and largest Rhat values. See [rstan::check_hmc_diagnostics()] and
#' [posterior::default_convergence_measures()] for details.
#'
#' @param x A `dynamitefit`object.
#' @param n \[`integer(1)`]\cr How many rows to print in
#'   parameter-specific convergence measures. Default is 1.
#' @export
mcmc_diagnostics <- function(x, n = 1) {
  if (!is.null(x$stanfit)) {
    if (x$stanfit@stan_args[[1]]$algorithm == "NUTS") {
      rstan::check_hmc_diagnostics(x$stanfit)
    }

    sumr <- posterior::summarise_draws(suppressWarnings(as_draws(x)),
      posterior::default_convergence_measures())

    cat("\nSmallest bulk-ESS values: \n")
    sumr |>
      dplyr::select(variable, ess_bulk) |>
      dplyr::arrange(ess_bulk) |>
      head(n_vars) |>
      tidyr::pivot_wider(names_from = variable, values_from = ess_bulk) |>
      print()
    cat("\nSmallest tail-ESS values: \n")
    sumr |>
      dplyr::select(variable, ess_tail) |>
      dplyr::arrange(ess_tail) |>
      head(n_vars) |>
      tidyr::pivot_wider(names_from = variable, values_from = ess_tail) |>
      print()
    cat("\nLargest Rhat values: \n")
    sumr |>
      dplyr::select(variable, rhat) |>
      dplyr::arrange(dplyr:::desc(rhat)) |>
      head(n_vars) |>
      tidyr::pivot_wider(names_from = variable, values_from = rhat) |>
      print()
  } else {
    message("No Stan model fit is available.")
  }

}

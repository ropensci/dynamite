#' Diagnostic Values of a Dynamite Model
#'
#' Prints HMC diagnostics, and lists parameters with smallest effective sample
#' sizes and largest Rhat values. See [rstan::check_hmc_diagnostics()] and
#' [posterior::default_convergence_measures()] for details.
#'
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param n \[`integer(1)`: \sQuote{1L}]\cr How many rows to print in
#'   parameter-specific convergence measures. Default is 1.
#' @export
#' @examples
#' mcmc_diagnostics(gaussian_example_fit)
mcmc_diagnostics <- function(x, n = 1L) {
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  n <- try_type(n, "integer")[1L]
  if (!is.null(x$stanfit)) {
    if (x$stanfit@stan_args[[1L]]$algorithm == "NUTS") {
      rstan::check_hmc_diagnostics(x$stanfit)
    }
    #TODO ok to suppress warnings?
    sumr <- posterior::summarise_draws(suppressWarnings(as_draws(x)),
      posterior::default_convergence_measures())

    cat("\nSmallest bulk-ESS values: \n")
    sumr |>
      dplyr::select(.data$variable, .data$ess_bulk) |>
      dplyr::arrange(.data$ess_bulk) |>
      utils::head(n) |>
      tidyr::pivot_wider(names_from = .data$variable,
                         values_from = .data$ess_bulk) |>
      print()
    cat("\nSmallest tail-ESS values: \n")
    sumr |>
      dplyr::select(.data$variable, .data$ess_tail) |>
      dplyr::arrange(.data$ess_tail) |>
      utils::head(n) |>
      tidyr::pivot_wider(names_from = .data$variable,
                         values_from = .data$ess_tail) |>
      print()
    cat("\nLargest Rhat values: \n")
    sumr |>
      dplyr::select(.data$variable, .data$rhat) |>
      dplyr::arrange(dplyr::desc(.data$rhat)) |>
      utils::head(n) |>
      tidyr::pivot_wider(names_from = .data$variable,
                         values_from = .data$rhat) |>
      print()
  } else {
    message_("No Stan model fit is available.")
  }
}

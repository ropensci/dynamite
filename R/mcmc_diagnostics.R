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
    sumr <- posterior::summarise_draws(
      suppressWarnings(as_draws(x)),
      posterior::default_convergence_measures()
    )

    cat("\nSmallest bulk-ESS values: \n")
    sumr |>
      dplyr::select(.data$variable, .data$ess_bulk) |>
      dplyr::arrange(.data$ess_bulk) |>
      utils::head(n) |>
      tidyr::pivot_wider(
        names_from = .data$variable,
        values_from = .data$ess_bulk
      ) |>
      print()
    cat("\nSmallest tail-ESS values: \n")
    sumr |>
      dplyr::select(.data$variable, .data$ess_tail) |>
      dplyr::arrange(.data$ess_tail) |>
      utils::head(n) |>
      tidyr::pivot_wider(
        names_from = .data$variable,
        values_from = .data$ess_tail
      ) |>
      print()
    cat("\nLargest Rhat values: \n")
    sumr |>
      dplyr::select(.data$variable, .data$rhat) |>
      dplyr::arrange(dplyr::desc(.data$rhat)) |>
      utils::head(n) |>
      tidyr::pivot_wider(
        names_from = .data$variable,
        values_from = .data$rhat
      ) |>
      print()
  } else {
    cat("No Stan model fit is available.")
  }
  invisible(x)
}

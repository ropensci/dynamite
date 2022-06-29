#' Print a Summary of a Dynamite Object
#'
#' Prints the summary information of the estimated dynamite model: The smallest
#' effective sample sizes, largest Rhat and summary statistics of the
#' time-invariant model parameters.
#'
#' @param x \[`dynamitefit`] The model fit object.
#' @param ... Further arguments to the print method for tibbles.
#'   See [tibble::formatting].
#' @method print dynamitefit
#' @export
#' @srrstats {BS6.0} *Software should implement a default `print` method for return objects*
#' @srrstats {BS5.3} *Bayesian Software should return convergence statistics or equivalent*
#' @srrstats {BS5.5} *Appropriate diagnostic statistics to indicate absence of convergence should either be returned or immediately able to be accessed.*
#' @srrstats {RE4.17} *Model objects returned by Regression Software should implement or appropriately extend a default `print` method which provides an on-screen summary of model (input) parameters and (output) coefficients.*
#' @examples
#' print(gaussian_example_fit)
print.dynamitefit <- function(x, ...) {
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  if (!is.null(x$stanfit)) {
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
    print(draws |> dplyr::select(dplyr::matches("([^\\]])$", perl = TRUE)) |>
      posterior::summarise_draws(), ...)

  } else {
    message_("No Stan model fit is available.")
  }
  invisible(x)
}

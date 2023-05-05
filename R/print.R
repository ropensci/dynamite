#' Print a Summary of a Dynamite Model Fit Object
#'
#' Information on the estimated dynamite model can be obtained via
#'  `print` including the following: The model formula, the data, the smallest
#' effective sample sizes, largest Rhat and summary statistics of the
#' time-invariant model parameters.
#'
#' @export
#' @rdname dynamite
#' @return `print` returns `x` invisibly.
#' @srrstats {BS6.0, RE4.17} Implements the `print` method for the model fit.
#' @srrstats {BS5.3, BS5.5} Contains convergence statistics in the output.
#' @examples
#' print(gaussian_example_fit)
#'
print.dynamitefit <- function(x, ...) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  cat("Model:\n")
  attr(x$dformulas$all, "random") <- attr(x$dformulas$stoch, "random")
  print.dynamiteformula(x$dformulas$all)
  cat(
    "\nData: ", x$data_name, " (Number of observations: ", nobs(x), ")",
    sep = ""
  )
  if (!is.null(x$group_var)) {
    cat(
      "\nGrouping variable: ", x$group_var, " (Number of groups: ",
      n_unique(x$data[[x$group_var]]), ")",
      sep = ""
    )
  }
  cat(
    "\nTime index variable: ", x$time_var, " (Number of time points: ",
    n_unique(x$data[[x$time_var]]), ")\n",
    sep = ""
  )
  if (!is.null(x$stanfit)) {
    draws <- suppressWarnings(as_draws(x))
    sumr <- posterior::summarise_draws(
      draws,
      posterior::default_convergence_measures()
    )
    min_ess <- which.min(sumr$ess_bulk)
    cat("\nSmallest bulk-ESS: ", round(sumr$ess_bulk[min_ess]), " (",
      sumr$variable[min_ess], ")",
      sep = ""
    )
    min_ess <- which.min(sumr$ess_tail)
    cat("\nSmallest tail-ESS: ", round(sumr$ess_tail[min_ess]), " (",
      sumr$variable[min_ess], ")",
      sep = ""
    )
    max_rhat <- which.max(sumr$rhat)
    cat("\nLargest Rhat: ", round(sumr$rhat[max_rhat], 3), " (",
      sumr$variable[max_rhat], ")",
      sep = ""
    )
    cat("\n\nElapsed time (seconds):\n")
    print(rstan::get_elapsed_time(x$stanfit))

    cat(
      "\nSummary statistics of the time- and group-invariant parameters:\n"
    )
    match_names <- grepl(
      pattern = "^(?!.*^nu|^omega|^lambda|.*\\[.*]).*",
      x = names(draws),
      perl = TRUE
    )
    print(posterior::summarise_draws(draws[, match_names]), ...)
  } else {
    cat("No Stan model fit is available.\n")
  }
  invisible(x)
}

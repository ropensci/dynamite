#' Print a Summary of a \pkg{dynamite} Model Fit Object
#'
#' Information on the estimated `dynamite` model can be obtained via
#' `print()` including the following: The model formula, the data,
#' the smallest effective sample sizes, largest Rhat and summary statistics of
#' the time-invariant and group-invariant model parameters.
#'
#' @export
#' @rdname dynamite
#' @param full_diagnostics By default, the effective sample size (ESS) and Rhat
#' are computed only for the time- and group-invariant parameters
#' (`full_diagnostics = FALSE`). Setting this to `TRUE` computes ESS and Rhat
#' values for all model parameters, which can take some time for complex models.
#' @return `print` returns `x` invisibly.
#' @srrstats {BS6.0, RE4.17} Implements the `print` method for the model fit.
#' @srrstats {BS5.3, BS5.5} Contains convergence statistics in the output.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' print(gaussian_example_fit)
#'
print.dynamitefit <- function(x, full_diagnostics = FALSE, ...) {
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
    cat("\n")
    mcmc_algorithm <- get_algorithm(x$stanfit) %in% c("NUTS", "hmc")
    if (mcmc_algorithm) {
      hmc_diagnostics(x)
    }
    # draws <- suppressWarnings(as_draws(x))
    # match_names <- grepl(
    #   pattern = "^(?!.*^nu|^omega|^lambda|.*\\[.*]).*",
    #   x = names(draws),
    #   perl = TRUE
    # )
    param_types <-  setdiff(all_types, c("nu", "omega", "lambda"))
    if (full_diagnostics && mcmc_algorithm) {
      # compute only the convergence measures for all variables
      draws <- suppressWarnings(as_draws(x))
      sumr <- posterior::summarise_draws(
        draws,
        posterior::default_convergence_measures()
      )
    } else {
      draws <- suppressWarnings(as_draws(x, types = param_types))
      sumr <- posterior::summarise_draws(draws)
    }
    if (mcmc_algorithm) {
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
      runtimes <- get_elapsed_time(x$stanfit)
      if (nrow(runtimes) > 2L) {
        rs <- rowSums(runtimes)
        cat("\n\nElapsed time (seconds) for fastest and slowest chains:\n")
        print(runtimes[c(which.min(rs), which.max(rs)), ])
      } else {
        cat("\n\nElapsed time (seconds):\n")
        print(runtimes)
      }
    }
    cat(
      "\nSummary statistics of the time- and group-invariant parameters:\n"
    )
    if (full_diagnostics) {
      draws <- suppressWarnings(as_draws(x, types = param_types))
      sumr <- posterior::summarise_draws(draws)
    }
    print(sumr, ...)
  } else {
    cat("No Stan model fit is available.\n")
  }
  invisible(x)
}

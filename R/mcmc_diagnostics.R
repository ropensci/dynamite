#' Diagnostic Values of a \pkg{dynamite} Model
#'
#' Prints HMC diagnostics and lists parameters with smallest effective sample
#' sizes and largest Rhat values. See [hmc_diagnostics()] and
#' [posterior::default_convergence_measures()] for details.
#'
#' @export
#' @family diagnostics
#' @rdname mcmc_diagnostics
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param n \[`integer(1)`]\cr How many rows to print in
#'   parameter-specific convergence measures. The default is 3. Should be a
#'   positive (unrestricted) integer.
#' @param ... Ignored.
#' @return Returns `x` (invisibly).
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' mcmc_diagnostics(gaussian_example_fit)
#'
mcmc_diagnostics <- function(x, ...) {
  UseMethod("mcmc_diagnostics", x)
}

#' @export
#' @rdname mcmc_diagnostics
mcmc_diagnostics.dynamitefit <- function(x, n = 3L, ...) {
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
  if (is.null(x$stanfit)) {
    cat("No Stan model fit is available.")
  } else {
    algorithm <- get_algorithm(x$stanfit)
    stopifnot_(
      algorithm %in% c("NUTS", "hmc"),
      "MCMC diagnostics are only meaningful for samples from MCMC.
       The model was estimated using the ", algorithm, "algorithm."
    )
    hmc_diagnostics(x)
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
  }
  invisible(x)
}

#' HMC Diagnostics for a \pkg{dynamite} Model
#'
#' Prints the divergences, saturated treedepths, and low E-BFMI warnings.
#'
#' @export
#' @family diagnostics
#' @rdname hmc_diagnostics
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param ... Ignored.
#' @return Returns `x` (invisibly).
#' data.table::setDTthreads(1) # For CRAN
#' hmc_diagnostics(gaussian_example_fit)
#'
hmc_diagnostics <- function(x, ...) {
  UseMethod("hmc_diagnostics", x)
}

#' @export
#' @rdname hmc_diagnostics
hmc_diagnostics.dynamitefit <- function(x, ...) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  if (is.null(x$stanfit)) {
    cat("No Stan model fit is available.")
  } else {
    algorithm <- get_algorithm(x$stanfit)
    stopifnot_(
      algorithm %in% c("NUTS", "hmc"),
      "MCMC diagnostics are only meaningful for samples from MCMC.
      Model was estimated using the ", algorithm, "algorithm."
    )
    n_draws <- ndraws(x)
    diags <- get_diagnostics(x$stanfit)
    n_divs <- diags$num_divergent
    n_trees <- diags$num_max_treedepth
    bfmis <- diags$ebfmi
    all_ok <- all(n_divs == 0L) && all(n_trees == 0L) && all(bfmis > 0.2)
    cat("NUTS sampler diagnostics:\n")
    all_ok_str <- ifelse_(
      all_ok,
      "\nNo divergences, saturated max treedepths or low E-BFMIs.\n",
      ""
    )
    cat(all_ok_str)
    div_str <- ifelse_(
      any(n_divs > 0L),
      paste0(
        "\n", sum(n_divs), " out of ", n_draws, " iterations ended with a ",
        "divergence. See Stan documentation for details.\n"
      ),
      ""
    )
    cat(div_str)
    mt <- get_max_treedepth(x$stanfit)
    mt <- ifelse_(is.null(mt), 10, mt)
    trees_str <- ifelse_(
      any(n_trees > 0L),
      paste0(
        "\n", sum(n_trees), " out of ", n_draws, " saturated the maximum ",
        "tree depth of ", mt, ". See Stan documentation for details.\n"
      ),
      ""
    )
    cat(trees_str)
    bfmis_str <- ifelse_(
      any(bfmis < 0.2),
      paste0(
        "\nChain(s) ", cs(which(bfmis < 0.2)), " had E-BFMI below 0.2, ",
        "indicating possible issues. See Stan documentation for details.\n"
      ),
      ""
    )
    cat(bfmis_str)
  }
  invisible(x)
}

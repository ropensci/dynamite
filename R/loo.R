#' Approximate Leave-One-Out (LOO) Cross-validation
#'
#' Estimates the leave-one-out (LOO) information criterion for `dynamite`
#' models using Pareto smoothed importance sampling with the `loo` package.
#'
#' @export
#' @export loo
#' @family diagnostics
#' @aliases loo
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param separate_channels \[`logical(1)`]\cr If `TRUE`, computes LOO
#'   separately for each channel. This can be useful in diagnosing where the
#'   model fails. Default is `FALSE`, in which case the likelihoods of
#'   different channels are combined, i.e., all channels of are left out.
#' @param thin \[`integer(1)`]\cr Use only every `thin` posterior sample when
#'   computing LOO. This can be beneficial with when the model object contains
#'   large number of samples. Default is `1` meaning that all samples are used.
#' @param ... Ignored.
#' @return An output from [loo::loo()] or a list of such outputs (if
#'   `separate_channels` was `TRUE`).
#' @references Aki Vehtari, Andrew, Gelman, and Johah Gabry (2017).
#' Practical Bayesian model evaluation using leave-one-out cross-validation and
#' WAIC. Statistics and Computing. 27(5), 1413--1432.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' \donttest{
#' # Please update your rstan and StanHeaders installation before running
#' # on Windows
#' if (!identical(.Platform$OS.type, "windows")) {
#'   # this gives warnings due to the small number of iterations
#'   suppressWarnings(loo(gaussian_example_fit))
#'   suppressWarnings(loo(gaussian_example_fit, separate_channels = TRUE))
#' }
#' }
#'
loo.dynamitefit <- function(x, separate_channels = FALSE, thin = 1L, ...) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  stopifnot_(
    is.null(x$imputed),
    "Leave-one-out cross-validation is not supported for models estimated using
     multiple imputation."
  )
  stopifnot_(
    !is.null(x$stanfit),
    "No Stan model fit is available."
  )
  stopifnot_(
    checkmate::test_flag(x = separate_channels),
    "Argument {.arg separate_channels} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_int(x = thin, lower = 1L, upper = ndraws(x)),
    "Argument {.arg thin} must be a single positive {.cls integer}."
  )
  n_chains <- x$stanfit@sim$chains
  n_draws <- ndraws(x)
  idx_draws <- seq.int(1L, n_draws, by = thin)
  # need equal number of samples per chain
  idx_draws <- idx_draws[seq_len(n_draws %/% thin - n_draws %% n_chains)]
  n_draws <- length(idx_draws) %/% n_chains
  out <- initialize_predict(
    x,
    newdata = NULL,
    type = "mean",
    eval_type = "loglik",
    funs = list(),
    impute = "none",
    new_levels = "none",
    global_fixed = FALSE,
    idx_draws,
    expand = FALSE,
    df = FALSE
  )$simulated
  # avoid NSE notes from R CMD check
  patterns <- NULL
  loo_ <- function(ll, n_draws, n_chains) {
    ll <- t(matrix(ll, ncol = n_draws * n_chains))
    reff <- loo::relative_eff(
      exp(ll),
      chain_id = rep(seq_len(n_chains), each = nrow(ll) / n_chains)
    )
    loo::loo(ll, r_eff = reff)
  }
  if (separate_channels) {
    ll <- split(
      x = data.table::melt(
        out,
        measure.vars = patterns("_loglik$"),
        na.rm = TRUE
      ),
      by = "variable"
    )
    lapply(ll, function(x) loo_(x$value, n_draws, n_chains))
  } else {
    temp <- out[, .SD, .SDcols = patterns("_loglik$")]
    ll <- temp[is.finite(rowSums(temp))][,
      rowSums(.SD),
      .SDcols = names(temp)
    ]
    loo_(ll, n_draws, n_chains)
  }
}

#' Approximate Leave-One-Out (LOO) Cross-validation
#'
#' Estimates the leave-one-out (LOO) information criterion for `dynamite`
#' models using Pareto smoothed importance sampling with the `loo` package.
#'
#' @export
#' @export loo
#' @aliases loo
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param separate_channels \[`logical(1)`]\cr If `TRUE`, computes LOO
#'   separately for each channel. This can be useful in diagnosing where the
#'   model fails. Default is `FALSE`, in which case the likelihoods of
#'   different channels are combined, i.e. all channels of are left out.
#' @param ... Ignored.
#' @return An output from [loo::loo()] or a list of such outputs (if
#'   `separate_channels` was `TRUE`).
#' @examples
#' \dontrun{
#' # this gives warnings due to the small number of iterations
#' loo(gaussian_example_fit)
#' loo(gaussian_example_fit, separate_channels = TRUE)
#' }
loo.dynamitefit <- function(x, separate_channels = FALSE, ...) {
  stopifnot_(
    !is.null(x$stanfit),
    "No Stan model fit is available."
  )
  stopifnot_(
    checkmate::test_flag(x = separate_channels),
    "Argument {.arg separate_channels} must be a single {.cls logical} value."
  )
  out <- initialize_predict(
    x,
    newdata = NULL,
    type = "mean",
    eval_type = "loglik",
    funs = list(),
    impute = "none",
    new_levels = "none",
    global_fixed = FALSE,
    n_draws = NULL,
    expand = FALSE
  )$simulated

  n_chains <- x$stanfit@sim$chains
  n_draws <- ndraws(x) / n_chains

  loo_ <- function(ll, n_draws, n_chains) {
    ll <- t(matrix(ll, ncol = n_draws * n_chains))
    reff <- loo::relative_eff(exp(ll), chain_id = rep(1:n_chains, each = n_draws))
    loo::loo(ll, r_eff = reff)
  }

  if (separate_channels) {
    ll <- out |>
      tidyr::pivot_longer(
        dplyr::ends_with("_loglik"), names_pattern  = "(.*)_loglik",
        values_drop_na = TRUE
      ) |>
      split(f = ~name)
    lapply(ll, function(x) loo_(x$value, n_draws, n_chains))
  } else {
    ll <- out |> tidyr::drop_na() |>
      dplyr::mutate(loglik = rowSums(dplyr::across(
        dplyr::ends_with("_loglik")))) |>
      dplyr::pull(.data$loglik)
    loo_(ll, n_draws, n_chains)
  }

}

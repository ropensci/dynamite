#' Approximate Leave-Future-Out (LFO) Cross-validation
#'
#' Estimates the leave-future-out (LFO) information criterion for `dynamite`
#' models using Pareto smoothed importance sampling.
#'
#' For multichannel models, the log-likelihoods of all channels are combined.
#' For models with groups, expected log predictive densities (ELPDs) are
#' computed independently for each group, but the re-estimation of the model
#' is triggered if pareto k values of any group exceeds the threshold.
#'
#' @export
#' @export lfo
#' @aliases lfo
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param L  \[`integer(1)`]\cr Positive integer defining how many time points
#'   should be used for the initial fit.
#' @param verbose \[`logical(1)`]\cr If `TRUE` (default), print the progress of
#'   the LFO computations to the console.
#' @param k_threshold \[`numeric(1)`]\cr Threshold for the pareto k estimate
#'   triggering refit. Default is 0.7.
#' @param ... Additional parameters to `dynamite`.
#' @return An `lfo` object which is a `list` with the following components:
#'   * `ELPD` Expected log predictive density estimate.
#'   * `ELPD_SE` Standard error of ELPD. This is a crude approximation which
#'      does not take into account potential serial correlations.
#'   *  `pareto_k` Pareto k values.
#'   *  `refits` Time points where model was re-estimated.
#'   *  `L` L value used in the LFO estimation.
#'   *  `k_threshold` Threshold used in the LFO estimation.
#' @examples
#' \dontrun{
#' # this gives warnings due to the small number of iterations
#' out <- lfo(gaussian_example_fit, L = 20)
#' out$ELPD
#' out$ELPD_SE
#' }
lfo <- function(x, L, verbose = TRUE, k_threshold = 0.7, ...) {

  stopifnot_(
    !is.null(x$stanfit),
    "No Stan model fit is available."
  )

  stopifnot_(
    checkmate::test_flag(x = verbose),
    "Argument {.arg verbose} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_number(x = k_threshold),
    "Argument {.arg k_threshold} must be a single {.cls numeric} value."
  )

  log_sum_exp <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  }
  log_mean_exp <- function(x) {
    log_sum_exp(x) - log(length(x))
  }

  T_ <- x$stan$sampling_vars$T

  stopifnot_(
    checkmate::test_int(x = L,
      lower = 0,
      upper = T_),
    "Argument {.arg L} must be a single {.cls integer} value between 0 and {T_}."
  )

  responses <- get_responses(x$dformulas$stoch)
  time <- x$time_var
  id <- x$group_var
  timepoints <- sort(unique(x$data[[time]]))
  d <- data.table::copy(x$data)
  set_na <- d[[x$time_var]] > timepoints[L]
  d[set_na, (responses) := NA]

  if (verbose) message_(paste0("Estimating model with ", L, " time points."))
  fit <- update(x, data = d, refresh = 0, ...)

  # would be faster to use only data
  # x$data[eval(time) >= timepoints[L] - x$stan$fixed]
  # but in a case of missing data this is not necessarily enough
  out <- initialize_predict(
    fit,
    newdata = x$data,
    type = "mean",
    eval_type = "loglik",
    funs = list(),
    impute = "none",
    new_levels = "none",
    global_fixed = FALSE,
    n_draws = NULL,
    expand = FALSE
  )$simulated

  n_draws <- ndraws(x)
  # sum the log-likelihood over the channels and non-missing time points
  # for each group, time, and draw
  # drop those id&time pairs which contain NA
  lls <- out |>
    dplyr::filter(.data[[time]] > timepoints[L]) |>
    dplyr::mutate(loglik = rowSums(dplyr::across(
      dplyr::ends_with("_loglik")))) |>
    dplyr::select(-dplyr::ends_with("_loglik")) |>
    tidyr::drop_na()

  elpds <- vector("list", T_ - L)
  elpds[[1]] <- lls |>
    dplyr::filter(.data[[time]] == timepoints[L + 1]) |>
    dplyr::group_by(.data[[time]], onlyif(!is.null(id), .data[[id]])) |>
    dplyr::summarise(elpd = log_mean_exp(.data$loglik), .groups = "keep") |>
    dplyr::pull(.data$elpd)

  i_refit <- L
  refits <- timepoints[L]
  ks <- vector("list", T_ - L - 1)

  for (i in (L + 1):(T_ - 1)) {

    if (nrow(lls |> dplyr::filter(.data[[time]]== i + 1)) > 0) {
      logratio <- lls |>
        dplyr::filter(.data[[time]] > timepoints[i_refit] & .data[[time]] <= timepoints[i]) |>
        dplyr::group_by(onlyif(!is.null(id), .data[[id]]), .data$.draw) |>
        dplyr::summarise(logratio = sum(.data$loglik), .groups = "keep")

      psis_obj <- suppressWarnings(loo::psis(
        matrix(logratio$logratio, nrow = n_draws), r_eff = NA
      ))
      k <- loo::pareto_k_values(psis_obj)
      ks[[i - L]] <- k
      if (any(k > k_threshold)) {
        if (verbose) message_(paste0("Estimating model with ", i, " time points"))
        # refit the model based on the first i time points
        i_refit <- i
        refits <- c(refits, timepoints[i])
        d <- data.table::copy(x$data)
        set_na <- d[[x$time_var]] > timepoints[i]
        d[set_na, (responses) := NA]
        fit <- update(fit, data = d, refresh = 0, ...)

        out <- initialize_predict(
          fit,
          newdata = x$data,
          type = "mean",
          eval_type = "loglik",
          funs = list(),
          impute = "none",
          new_levels = "none",
          global_fixed = FALSE,
          n_draws = NULL,
          expand = FALSE
        )$simulated

        lls <- out |>
          dplyr::filter(.data[[time]] > timepoints[L]) |>
          tidyr::drop_na() |>
          dplyr::mutate(loglik = rowSums(dplyr::across(
            dplyr::ends_with("_loglik")))) |>
          dplyr::select(-dplyr::ends_with("_loglik"))

        elpds[[i - L + 1]] <-  lls |>
          dplyr::filter(.data[[time]] == timepoints[i + 1]) |>
          dplyr::group_by(.data[[time]], onlyif(!is.null(id), .data[[id]])) |>
          dplyr::summarise(elpd = log_mean_exp(.data$loglik), .groups = "keep") |>
          dplyr::pull(.data$elpd)

      } else {
        lw <- loo::weights.importance_sampling(psis_obj, normalize = TRUE)
        ll <- lls |>
          dplyr::filter(.data[[time]] == timepoints[i + 1]) |>
          dplyr::pull(.data$loglik)
        elpds[[i - L + 1]] <- log_sum_exp_rows(t(lw) + matrix(ll, ncol = n_draws))
      }
    } else {
      # no observations
      ks[[i - L]] <- NA
      elpds[[i - L + 1]] <- NA
    }
  }
  elpds <- unlist(elpds)
  out <- list(ELPD = sum(elpds), ELPD_SE = sd(elpds) * sqrt(length(elpds)),
    pareto_k = ks, ELPDs = elpds, refit_times = timepoints[refits], L = L,
    k_threshold = k_threshold)

  class(out) <- "lfo"
  out
}
#' Print the results from the LFO
#'
#' Prints the summary of the leave-future-out cross-validation.
#' @param x x \[`lfo`]\cr Output from `lfo` function.
#' @param ... Ignored.
#' @return Returns `x` invisibly.
#' @export
#' @examples
#' \dontrun{
#' lfo(gaussian_example_fit, L = 20)
#' }
print.lfo <- function(x, ...) {
  cat("\nApproximate LFO starting from time point", x$L)
  cat("\nModel was re-estimated at time points ",
    paste(x$refit_times, collapse = ", "),
    " (Based on Pareto k threshold of ", x$k_threshold, ")\n", sep = "")
  cat("\nEstimated expected log predictive density (ELPD):", x$ELPD)
  cat("\nStandard error estimate of the ELPD:", x$ELPD_SE)
  invisible(x)
}

#' Diagnostic Plot for Pareto k Values from LFO
#'
#' Plots Pareto k values per each time point (with one point per group),
#' together with the horizontal line representing the used threshold.
#' @param x x \[`lfo`]\cr Output from `lfo` function.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' plot(lfo(gaussian_example_fit, L = 20))
#' }
plot.lfo <- function(x, ...) {
  d <- data.frame(k = unlist(x$pareto_k),
    time = rep(x$L + 1:length(x$pareto_k),
      times = lengths(x$pareto_k)))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$time, y = .data$k)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$k > x$k_threshold),
      shape = 3, show.legend = FALSE, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = x$k_threshold,
      linetype = 2, color = "red2") +
    ggplot2::scale_color_manual(values = c("cornflowerblue", "darkblue")) +
    ggplot2::labs(x = "Time", y = "Pareto k")
}


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
#' @family diagnostics
#' @aliases lfo
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param L  \[`integer(1)`]\cr Positive integer defining how many time points
#'   should be used for the initial fit.
#' @param verbose \[`logical(1)`]\cr If `TRUE` (default), print the progress of
#'   the LFO computations to the console.
#' @param k_threshold \[`numeric(1)`]\cr Threshold for the Pareto k estimate
#'   triggering refit. Default is 0.7.
#' @param ... Additional arguments passed to [rstan::sampling()] or
#'   [cmdstanr::sample()], such as `chains` and `cores` (`parallel_chains` in
#'   `cmdstanr`).
#' @return An `lfo` object which is a `list` with the following components:
#'
#'   * `ELPD`\cr Expected log predictive density estimate.
#'   * `ELPD_SE`\cr Standard error of ELPD. This is a crude approximation which
#'      does not take into account potential serial correlations.
#'   * `pareto_k`\cr Pareto k values.
#'   * `refits`\cr Time points where model was re-estimated.
#'   * `L`\cr L value used in the LFO estimation.
#'   * `k_threshold`\cr Threshold used in the LFO estimation.
#'
#' @references Paul-Christian Bürkner, Jonah Gabry, and Aki Vehtari (2020).
#' Approximate leave-future-out cross-validation for Bayesian time series
#' models, Journal of Statistical Computation and Simulation, 90:14, 2499-2523.
#' @examples
#' \donttest{
#' # Please update your rstan and StanHeaders installation before running
#' # on Windows
#' if (!identical(.Platform$OS.type, "windows")) {
#'   # this gives warnings due to the small number of iterations
#'   out <- suppressWarnings(
#'     lfo(gaussian_example_fit, L = 20, chains = 1, cores = 1)
#'   )
#'   out$ELPD
#'   out$ELPD_SE
#' }
#' }
#'
lfo <- function(x, L, verbose = TRUE, k_threshold = 0.7, ...) {
  stopifnot_(
    !is.null(x$stanfit),
    "No Stan model fit is available."
  )
  stopifnot_(
    !missing(L) && checkmate::test_int(x = L, lower = 1L),
    "Argument {.arg L} must be a single positive {.cls integer}."
  )
  stopifnot_(
    checkmate::test_flag(x = verbose),
    "Argument {.arg verbose} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_number(x = k_threshold),
    "Argument {.arg k_threshold} must be a single {.cls numeric} value."
  )
  time_var <- x$time_var
  group_var <- x$group_var
  timepoints <- sort(unique(x$data[[time_var]]))
  T_ <- length(timepoints)
  stopifnot_(
    checkmate::test_int(x = L, lower = 0, upper = T_),
    "Argument {.arg L} must be a single {.cls integer} between 0 and {T_}."
  )
  responses <- get_responses(x$dformulas$stoch)
  d <- data.table::copy(x$data)
  set_na_ <- d[[time_var]] > timepoints[L]
  # d[set_na, (responses) := NA, env = list(set_na = set_na)]
  d[set_na_, (responses) := NA]
  if (verbose) {
    message_("Estimating model with {L} time points.")
  }
  fit <- update_(x, data = d, refresh = 0, ...)

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
    expand = FALSE,
    df = FALSE
  )$simulated
  # avoid NSE notes from R CMD check
  loglik <- patterns <- .draw <- NULL
  n_draws <- ndraws(x)
  # sum the log-likelihood over the channels and non-missing time points
  # for each group, time, and draw
  # drop those id&time pairs which contain NA
  subset_indices_ <- out[[time_var]] > timepoints[L]
  lls <- #stats::na.omit(
    out[
    subset_indices_
    # time > timepoints[L], #,
    # env = list(time = time, timepoints = timepoints, L = L)
  ][,
    loglik := rowSums(.SD),
    .SDcols = patterns("_loglik$")
  ][,
    .SD,
    .SDcols = !patterns("_loglik$")
  ]
  #)
  elpds <- vector("list", T_ - L)
  subset_index_ <- lls[[time_var]] == timepoints[L + 1L]
  elpds[[1L]] <- stats::na.omit(lls[
    subset_index_
    # time == timepoints[L + 1L], ,
    # env = list(time = time, timepoints = timepoints, L = L)
  ][,
    list(elpd = log_mean_exp(loglik)),
    by = c(time_var, group_var)
    # by = list(time, id)#,
    # env = list(log_mean_exp = "log_mean_exp", time = time, id = id)
  ][["elpd"]])
  i_refit <- L
  refits <- timepoints[L]
  ks <- vector("list", T_ - L - 1L)
  for (i in seq.int(L + 1L, T_ - 1L)) {
    subset_index_ <- lls[[time_var]] == timepoints[i + 1L]
    if (lls[subset_index_, .N] > 0L) {
      # .N,
      # env = list(time = time, i = i)] > 0L) {
      logratio_subset_index_ <-
        lls[[time_var]] > timepoints[i_refit] &
          lls[[time_var]] <= timepoints[i]
      logratio <- lls[
        logratio_subset_index_, # ,
        # env = list(
        #  time = time,
        #  timepoints = timepoints,
        #  i = i,
        #  i_refit = i_refit
        # )
      ][,
        list(logratio = sum(loglik, na.rm = TRUE)),
        # by = list(id, .draw),
        by = c(group_var, ".draw")
      ]
      lr <- matrix(logratio$logratio, nrow = n_draws, byrow = TRUE)
      lw_subset_index_ <- lls[[time_var]] == timepoints[i + 1]
      ll <- matrix(
        lls[
          lw_subset_index_,
          loglik # ,
          # env = list(time = time, timepoints = timepoints, i = i)
        ], nrow = n_draws, byrow = TRUE
      )
      non_na_idx <- intersect(which(!is.na(colSums(lr))), which(!is.na(colSums(ll))))
      lr <- lr[, non_na_idx]
      ll <- ll[, non_na_idx]
      psis_obj <- suppressWarnings(
        loo::psis(
          lr,
          r_eff = NA
        )
      )
      k <- loo::pareto_k_values(psis_obj)
      ks[[i - L]] <- k
      if (any(k > k_threshold)) {
        if (verbose) {
          message_("Estimating model with {i} time points.")
        }
        # refit the model based on the first i time points
        i_refit <- i
        refits <- c(refits, timepoints[i])
        d <- data.table::copy(x$data)
        set_na_ <- d[[x$time_var]] > timepoints[i]
        # d[set_na, (responses) := NA, env = list(set_na = set_na)]
        d[set_na_, (responses) := NA]
        fit <- update_(fit, data = d, refresh = 0, ...)
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
          expand = FALSE,
          df = FALSE
        )$simulated
        threshold_subset_index_ <- out[[time_var]] > timepoints[L]
        lls <- out[
          threshold_subset_index_, # ,
          # env = list(time = time, timepoint = timepoints, L = L)
        ][,
          loglik := rowSums(.SD),
          .SDcols = patterns("_loglik$")
        ][,
          .SD,
          .SDcols = !patterns("_loglik$")
        ]
        elpds_subset_index_ <- lls[[time_var]] == timepoints[i + 1L]
        elpds[[i - L + 1L]] <- stats::na.omit(lls[
          elpds_subset_index_ # ,
          # env = list(time = time, timepoints = timepoints, i = i)
        ][,
          list(elpd = log_mean_exp(loglik)),
          # by = list(time, id)#,
          by = c(time_var, group_var)
          # env = list(log_mean_exp = "log_mean_exp", time = time, id = id)
        ][["elpd"]])
      } else {
        lw <- loo::weights.importance_sampling(psis_obj, normalize = TRUE)
        elpds[[i - L + 1L]] <-
          log_sum_exp_rows(
            t(lw + ll),
            ncol(lw),
            n_draws
          )
      }
    } else {
      # no observations
      ks[[i - L]] <- NA
      elpds[[i - L + 1L]] <- NA
    }
  }
  elpds <- unlist(elpds)
  structure(
    list(
      ELPD = sum(elpds),
      ELPD_SE = sd(elpds) * sqrt(length(elpds)),
      pareto_k = ks,
      ELPDs = elpds,
      refit_times = refits,
      L = L,
      k_threshold = k_threshold
    ),
    class = "lfo"
  )
}
#' Print the results from the LFO
#'
#' Prints the summary of the leave-future-out cross-validation.
#' @param x x \[`lfo`]\cr Output of the `lfo` method.
#' @param ... Ignored.
#' @return Returns `x` invisibly.
#' @export
#' @examples
#' \donttest{
#' # Please update your rstan and StanHeaders installation before running
#' # on Windows
#' if (!identical(.Platform$OS.type, "windows")) {
#'   # This gives warnings due to the small number of iterations
#'   suppressWarnings(lfo(gaussian_example_fit, L = 20))
#' }
#' }
#'
print.lfo <- function(x, ...) {
  cat("\nApproximate LFO starting from time point", x$L)
  cat(
    "\nModel was re-estimated at time points ",
    paste(x$refit_times, collapse = ", "),
    " (Based on Pareto k threshold of ", x$k_threshold, ")\n",
    sep = ""
  )
  cat("\nEstimated expected log predictive density (ELPD):", x$ELPD)
  cat("\nStandard error estimate of the ELPD:", x$ELPD_SE)
  invisible(x)
}

#' Diagnostic Plot for Pareto k Values from LFO
#'
#' Plots Pareto k values per each time point (with one point per group),
#' together with a horizontal line representing the used threshold.
#'
#' @param x \[`lfo`]\cr Output from the `lfo` function.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
#' @examples
#' \donttest{
#' # Please update your rstan and StanHeaders installation before running
#' # on Windows
#' if (!identical(.Platform$OS.type, "windows")) {
#'   # This gives warnings due to the small number of iterations
#'   plot(suppressWarnings(
#'     lfo(gaussian_example_fit, L = 20, chains = 1, cores = 1)
#'   ))
#' }
#' }
#'
plot.lfo <- function(x, ...) {
  d <- data.frame(
    k = unlist(x$pareto_k),
    time = rep(
      x$L + seq_len(length(x$pareto_k)),
      times = lengths(x$pareto_k)
    )
  )
  d$threshold <- d$k > x$k_threshold
  # avoid NSE notes from R CMD check
  time <- k <- threshold <- NULL
  ggplot2::ggplot(d, ggplot2::aes(x = time, y = k)) +
    ggplot2::geom_point(
      ggplot2::aes(color = threshold),
      shape = 3,
      show.legend = FALSE,
      alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = x$k_threshold,
      linetype = 2,
      color = "red2"
    ) +
    ggplot2::scale_color_manual(values = c("cornflowerblue", "darkblue")) +
    ggplot2::labs(x = "Time", y = "Pareto k")
}

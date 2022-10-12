#' Approximate Leave-Future-Out (LFO) Cross-validation
#'
#' Estimates the leave-future-out (LFO) information criterion for `dynamite`
#' models using Pareto smoothed importance sampling.
#'
#' @export
#' @export lfo
#' @aliases lfo
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param L Positive integer defining how many observations should be used for
#'   the initial fit.
#' @param verbose If \code{TRUE} (default), print the progress of the LFO
#'   computations to the console.
#' @param k_threshold Threshold for the pareto k estimate triggering refit. Default
#'   is 0.7.
#' @param ... Additional parameters to `dynamite`.
#' @return List with components \code{ELPD} (Expected log predictive density),
#'   \code{ELPDs} (observation-specific ELPDs for each time point),
#'   \code{ks} (Pareto k values), and \code{refits} (time points where model
#'   was re-estimated).
#' @examples
#' \dontrun{
#' # this gives warnings due to the small number of iterations
#' lfo(gaussian_example_fit, L = 20)
#' }
lfo <- function(x, L, verbose = TRUE, k_threshold = 0.7, ...) {

  log_sum_exp <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  }
  log_mean_exp <- function(x) {
    log_sum_exp(x) - log(length(x))
  }

  stopifnot_(
    !is.null(x$stanfit),
    "No Stan model fit is available."
  )
  T_ <- x$stan$sampling_vars$T
  responses <- get_responses(x$dformulas$stoch)
  time <- x$time_var
  id <- x$group_var


  d <- data.table::copy(x$data)
  d[eval(time) > L, (responses) := NA]


  if (verbose) message_(paste0("Estimating model with ", L, " time points."))
  fit <- update(x, data = d, refresh = 0, ...)

  out <- initialize_predict(
    fit,
    newdata = x$data[eval(time) >= L],
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
  # sum the log-likelihood over the channels for each group and draw
  lls <- out |> tidyr::drop_na() |>
    dplyr::mutate(loglik = rowSums(dplyr::across(
      dplyr::ends_with("_loglik")))) |>
      dplyr::select(-dplyr::ends_with("_loglik"))

  elpds <- vector("list", T_ - L)
  elpds[[1]] <- lls |>
    dplyr::filter(.data[[time]] == L + 1) |>
    dplyr::group_by(.data[[time]], onlyif(!is.null(id), .data[[id]])) |>
    dplyr::summarise(elpd = log_mean_exp(.data$loglik), .groups = "keep") |>
    dplyr::pull(elpd)

  i_refit <- L
  refits <- L
  ks <- vector("list", T_ - L - 1)

  for (i in (L + 1):(T_ - 1)) {
    logratio <- lls |>
      dplyr::filter(.data[[time]] > i_refit & .data[[time]] <= i) |>
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
      refits <- c(refits, i)
      d <- data.table::copy(x$data)
      d[eval(time) > i, (responses) := NA]
      fit <- update(fit, data = d, refresh = 0, ...)

      out <- initialize_predict(
        fit,
        newdata = x$data[eval(time) >= L],
        type = "mean",
        eval_type = "loglik",
        funs = list(),
        impute = "none",
        new_levels = "none",
        global_fixed = FALSE,
        n_draws = NULL,
        expand = FALSE
      )$simulated

      lls <- out |> tidyr::drop_na() |>
        dplyr::mutate(loglik = rowSums(dplyr::across(
          dplyr::ends_with("_loglik")))) |>
        dplyr::select(-dplyr::ends_with("_loglik"))

      elpds[[i - L + 1]] <-  lls |>
        dplyr::filter(.data[[time]] == i + 1) |>
        dplyr::group_by(.data[[time]], onlyif(!is.null(id), .data[[id]])) |>
        dplyr::summarise(elpd = log_mean_exp(.data$loglik), .groups = "keep") |>
        dplyr::pull(elpd)

    } else {
      lw <- loo::weights.importance_sampling(psis_obj, normalize = TRUE)
      ll <- lls |>
        dplyr::filter(.data[[time]] == i + 1) |>
        dplyr::pull(.data$loglik)
      elpds[[i - L + 1]] <- log_sum_exp_rows(t(lw) + matrix(ll, ncol = n_draws))
    }
  }

  out <- list(
    ELPD = sum(unlist(elpds)),
    ELPDs = elpds,
    pareto_k = ks,
    refit_times = refits,
    L = L, k_threshold = k_threshold)

  class(out) <- "lfo"
  out
}


plot.lfo <- function(x, ...) {
  d <- data.frame(k = unlist(x$pareto_k),
    time = rep(x$L + 1:length(x$pareto_k), times = lengths(x$pareto_k)))
  ggplot2::ggplot(d, ggplot2::aes(x = time, y = k)) +
    ggplot2::geom_point(ggplot2::aes(color = k > x$threshold),
      shape = 3, show.legend = FALSE, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = x$threshold,
      linetype = 2, color = "red2") +
    ggplot2::scale_color_manual(values = c("cornflowerblue", "darkblue")) +
    ggplot2::labs(x = "Time", y = "Pareto k")
}


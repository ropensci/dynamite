#' Update a Dynamite Model
#'
#' Note that using a different backend for the original model fit and when
#' updating can lead to an error due to different naming in `cmdstanr` and
#' `rstan` sampling arguments.
#'
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param dformula \[`dynamiteformula`]\cr Updated model formula. By default
#'   the original formula is used.
#' @param data
#'   \[`data.frame`, `tibble::tibble`, or `data.table::data.table`]\cr Data for
#'   the updated model. By default original data is used.
#' @param priors \[`data.frame`]\cr Updated priors. By default the priors of
#'   the original model are used.
#' @param recompile \[`logical(1)`]\cr Should the model be recompiled? If
#' `NULL` (default), tries to avoid recompilation. Recompilation is forced when
#'  the model formula or the priors are changed, or if the new data contains
#'  missing values in a channel which did not contain missing values in the
#'  original data. Recompilation is also forced in case the backend previous or
#'  new backend is `cmdstanr`.
#' @param ... Additional parameters to `dynamite`.
#' @return An updated `dynamitefit` object.
#' @export
#' @family fitting
#' @examples
#' \dontrun{
#' # re-estimate the example fit without thinning:
#' # As the model is compiled on Windows, this will fail on other platforms
#' if (identical(.Platform$OS.type, "windows")) {
#'   fit <- update(gaussian_example_fit, thin = 1)
#' }
#' }
#'
update.dynamitefit <- function(object, dformula = NULL, data = NULL,
                               priors = NULL, recompile = NULL, ...) {
  call <- object$call
  if (!is.null(dformula)) {
    call$dformula <- dformula
    recompile <- TRUE
  }
  if (!is.null(priors)) {
    call$priors <- priors
    recompile <- TRUE
  } else {
    if (is.null(dformula)) {
      # use original priors by default in case the formula did not change
      call$priors <- get_priors(object)
    }
  }
  if (!is.null(data)) {
    call$data <- substitute(data)
  }
  extras <- match.call(expand.dots = FALSE)$...
  # always ask for recompilation if using cmdstanr
  # however, this does not necessarily lead to recompilation if original model
  # was compiled in the same session
  if (object$backend == "cmdstanr" || identical(extras$backend, "cmdstanr")) {
    recompile <- TRUE
  }
  if (length(extras) > 0L) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) {
      call[[a]] <- extras[[a]]
    }
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  call$group <- ifelse_(
    is.null(call$group),
    object$group_var,
    call$group
  )
  # check that the model code does not change
  if (!isTRUE(recompile)) {
    call0 <- call
    call0$debug <- list(no_compile = TRUE, model_code = TRUE)
    recompile <- !identical(
      eval(call0)$model_code, as.character(get_code(object))
    )
  }
  if (!recompile) {
    call$debug <- list(stanfit = object$stanfit@stanmodel)
  }
  eval(call)
}

#' Internal `update` Method For LFO
#'
#' @inheritParams update.dynamitefit
#' @noRd
update_ <- function(object, data, refresh, ...) {
  call <- object$call
  call$dformula <- formula(object)
  call$data <- data
  call$time <- object$time_var
  call$group <- object$group_var
  call$priors <- get_priors(object)
  call$refresh <- refresh
  call$debug <- NULL
  recompile <- NULL
  extras <- match.call(expand.dots = FALSE)$...
  if (object$backend == "cmdstanr" || identical(extras$backend, "cmdstanr")) {
    recompile <- TRUE
  }
  # remove dynamite arguments
  extras[c("dformula", "data", "time", "group", "priors", "debug")] <- NULL
  if (length(extras) > 0L) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) {
      call[[a]] <- extras[[a]]
    }
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  e <- new.env()
  if (!isTRUE(recompile)) {
    call0 <- call
    call0$debug <- list(no_compile = TRUE, model_code = TRUE)
    recompile <- !identical(
      eval(call0, env = e)$model_code, as.character(get_code(object))
    )
  }
  if (!recompile) {
    call$debug <- list(stanfit = object$stanfit@stanmodel)
  }
  eval(call, env = e)
}

#' Update Dynamite Model
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
#'  the model formula or priors are changed, or if the new data contains
#'  missing values in a channel which did not contain missing values in the
#'  original data. Recompilation is also forced in case the backend previous or
#'  new backend is `cmdstanr`.
#' @param ... Additional parameters to `dynamite`.
#' @return Updated `dynamitefit` object.
#' @export
#' @examples
#' # re-estimate the example fit without thinning:
#' fit <- update(gaussian_example_fit, thin = 1)
update.dynamitefit <- function(object, dformula = NULL, data = NULL, priors = NULL,
  recompile = NULL, ...) {

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
    attr(data, "data_name") <- deparse1(substitute(data))
    call$data <- data
  }
  if (object$backend == "cmdstanr") {
    recompile <- TRUE
  }

  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
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

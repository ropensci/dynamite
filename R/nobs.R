#' Extract the Number of Observations Used to Fit a Dynamite Model
#'
#' @export
#' @family output
#' @param object \[`dynamitefit`]\cr The model fit object.
#' @param ... Not used.
#' @return Total number of non-missing observations as an `integer`.
#' @srrstats {RE4.5} Provides the number of observations.
#' @examples
#' nobs(gaussian_example_fit)
#'
nobs.dynamitefit <- function(object, ...) {
  stopifnot_(
    !missing(object),
    "Argument {.arg object} is missing."
  )
  stopifnot_(
    is.dynamitefit(object),
    "Argument {.var object} must be a {.cls dynamitefit} object."
  )
  sampling_vars <- get_data(object)
  cg <- attr(object$dformulas$all, "channel_groups")
  n_cg <- n_unique(cg)
  n_obs <- 0L
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    n_cg <- ifelse_(
      length(cg_idx) > 1L,
      paste0("n_obs_", object$stan$channel_group_vars[[i]]$y_cg),
      paste0("n_obs_", object$dformulas$all[[cg_idx]]$name)
    )
    n_obs <- n_obs + sum(sampling_vars[[n_cg]])
  }
  n_obs
}

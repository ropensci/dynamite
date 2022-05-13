#' Convert \code{dynamite} Output to \code{draws_df} Format
#'
#' Converts the output from \code{dynamite} call to a
#' \code{draws_df} format of the \code{posterior} package, enabling the use
#' of diagnostics and plotting methods of \code{posterior} and \code{bayesplot}
#' packages. Note that this function returns all variables in a wide format,
#' whereas [dynamite::as.data.frame()] uses the long format.
#'
#' @param x An object of class \code{dynamitefit}.
#' @inheritParams as.data.frame.dynamitefit
#' @return A \code{draws_df} object.
#' @importFrom posterior as_draws as_draws_df
#' @aliases as_draws as_draws_df
#' @export
#' @export as_draws_df
#' @rdname as_draws-dynamitefit
#' @method as_draws_df dynamitefit
#' @examples
#' as_draws(gaussian_example_fit, types = c("sigma", "beta"))
#'
as_draws_df.dynamitefit <- function(x, responses = NULL, types = NULL, ...) {
  d <- as.data.frame(x, responses, types, summary = FALSE) |>
    dplyr::select(.data$parameter, .data$value, .data$time, .data$category,
                  .data$.iteration, .data$.chain) |>
     tidyr::pivot_wider(
       values_from = .data$value,
       names_from = c(.data$parameter, .data$time, .data$category),
       names_glue = "{parameter}[{time}]_{category}")
  # remove NAs from time-invariant parameter names
  colnames(d) <- gsub("\\[NA\\]", "", colnames(d))
  # remove NAs from parameters which are not category specific
  colnames(d) <- gsub("_NA", "", colnames(d))
  d |> posterior::as_draws()
}
#' @export
#' @export as_draws
#' @rdname as_draws-dynamitefit
#' @method as_draws dynamitefit
as_draws.dynamitefit <- function(x, responses = NULL, types = NULL, ...) {
  as_draws_df.dynamitefit(x, responses, types, ...)
}

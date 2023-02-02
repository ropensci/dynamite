#' Convert `dynamite` Output to `draws_df` Format
#'
#' Converts the output from [dynamite::dynamite()] call to a
#' `draws_df` format of the \pkg{posterior} package, enabling the use
#' of diagnostics and plotting methods of \pkg{posterior} and \pkg{bayesplot}
#' packages. Note that this function returns variables in a wide format,
#' whereas [dynamite::as.data.frame()] uses the long format.
#'
#' You can use the arguments `parameters`, `responses` and `types` to extract
#' only a subset of the model parameters (i.e., only certain types of parameters
#'  related to a certain response variable).
#'
#' See potential values for the types argument in
#' [dynamite::as.data.frame.dynamitefit()] and
#' [dynamite::get_parameter_names()] for potential values for `parameters`
#' argument.
#' @export
#' @aliases as_draws_df
#' @export as_draws_df
#' @rdname as_draws-dynamitefit
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @inheritParams as.data.frame.dynamitefit
#' @return A `draws_df` object.
#' @examples
#' as_draws(gaussian_example_fit, types = c("sigma", "beta"))
#'
as_draws_df.dynamitefit <- function(x, parameters = NULL, responses = NULL,
  types = NULL, ...) {
  # avoid NSE notes from R CMD check
  .chain <- .draw <- .iteration <- NULL
  category <- group <- parameter <- response <- time <- type <- value <- NULL
  d <- as.data.table.dynamitefit(
    x,
    parameters = parameters,
    responses = responses,
    types = types,
    summary = FALSE,
    include_fixed = FALSE
  )[,
    .SD,
    .SDcols = c(
      "parameter",
      "value",
      "time",
      "category",
      "group",
      ".iteration",
      ".chain"
    )
  ][
    order(parameter, time, category, group, .chain, .iteration),
  ]
  dn <- unique(glue::glue("{d$parameter}[{d$time}]_{d$category}_id{d$group}"))
  d <- data.table::dcast(
    data = d,
    formula = .chain + .iteration ~ parameter + time + category + group,
    value.var = "value"
  )
  # remove NAs from time-invariant parameter names
  dn <- gsub("\\[NA\\]", "", dn)
  # remove NAs from parameters which are not category specific
  dn <- gsub("_NA", "", dn)
  # remove NAs from parameters which are not group specific
  dn <- gsub("_idNA", "", dn)
  colnames(d) <- c(".chain", ".iteration", dn)
  posterior::as_draws(
    data.table::setDF(d)
  )
}

#' @export
#' @export as_draws
#' @aliases as_draws
#' @rdname as_draws-dynamitefit
#' @return A `draws_df` object.
#' @inheritParams as_draws_df.dynamitefit
as_draws.dynamitefit <- function(x, responses = NULL, types = NULL, ...) {
  as_draws_df.dynamitefit(x, responses, types, ...)
}

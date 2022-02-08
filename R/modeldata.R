#' Define the data used by the model.
#' @param data A \code{data.frame} in long format,
#' @param override A logical value. If \code{FALSE} (the default), an existing
#'    definition for the model data will not be overridden by another call to \code{use_data}.
#'    If \code{TRUE}, any existing definitions will be replaced.
#' @export
# use_data <- function(data, replace = FALSE) {
#     # TODO some verification that data is in correct format
#     structure(data, replace = replace, class = "modeldata")
# }
#
# # Checks if argument represent model data
# is.modeldata <- function(x) {
#     inherits(x, "modeldata")
# }
#
# TODO Verify that this is unnecessary

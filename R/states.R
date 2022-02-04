#' Define the hidden states of the model.
#' @param n An integer value giving the number of hidden states.
#' @param fixed An integer vector with two elements giving
#'    the number of hidden states to fix at the head and tail.
#' @param individual A logical value. If \code{FALSE} (the default), then
#'    the subjects are governed by the same hidden state process. If \code{TRUE},
#'    each subject will have an individual process.
#' @param override A logical value. If \code{FALSE} (the default), an existing
#'    definition for the hidden states will not be overridden by another call to \code{states}.
#'    If \code{TRUE}, any existing definitions will be replaced.
#' @export
states <- function(n, fixed = c(1L, 1L), individual = FALSE, replace = FALSE) {
    if (missing(n)) {
        stop_("Please enter the number of hidden states 'n'")
    }
    n <- try_(n, type = "integer")[1]
    individual <- try_(individual, type = "logical")[1]
    fixed <- try_(fixed, type = "integer")[1:2]
    structure(
        list(n = n, fixed = fixed, individual = individual),
        replace = replace,
        class = "hiddenstates"
    )
}

# Checks if the argument represents a hidden states definition
is.hiddenstates <- function(x) {
    inherits(x, "hiddenstates")
}

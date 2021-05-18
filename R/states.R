#' @export
states <- function(n, fixed = c(1L, 1L), individual = FALSE) {
    if (missing(n)) {
        stop_("Please enter the number of hidden states 'n'")
    }
    if (!is.logical(individual)) {
        stop_("Argument 'individual' must be logical")
    }
    if (!is.integer(fixed)) {
        stop("Argument 'fixed' must be of type 'integer'")
    }
    structure(
        list(n = n[1], fixed = fixed[1:2], individual = individual[1]),
        class = "hiddenstates"
    )
}

is.hiddenstates <- function(x) {
    inherits(x, "hiddenstates")
}

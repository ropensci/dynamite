# TODO implement families

#' Family Functions for \pkg{because} Models
#'
#' @export
becausefamily <- function(name, ...) {
    # A wrapper to allow both:
    #    gaussian(...) and becausefamily("gaussian", ...)
    becausefamily_(name, ...)
}

# Internal use
becausefamily_ <- function(name, ...) {
    name <- tolower(as.character(name)[1])
    if (!is_supported(name)) {
        stop_(name, " is not a supported family")
    }
    # do something
    out <- list(name = name)
    class(out) <- c("becausefamily")
    out
}

# Checks if argument is a becausefamily object
is.becausefamily <- function(x) {
    inherits(x, "becausefamily")
}

#' @rdname becausefamily
#' @export
categorical <- function(...) {
    # do something different
    becausefamily_("categorical", ...)
}

#' @rdname becausefamily
#' @export
gaussian <- function(...) {
    # do something else
    becausefamily_("gaussian", ...)
}

# Hardcoded families for now

# Check if response of a family is continuous
is_continuous <- function(x) {
    x$name %in% continuous_distributions
}

# Check if response of a family is discrete
is_discrete <- function(x) {
    x$name %in% discrete_distributions
}

# Check if family is supported
is_supported <- function(name) {
    name %in% supported_families
}

continuous_distributions <- c(
    "gaussian",
    "gamma"
)

discrete_distributions <- c(
    "categorical",
    "binomial",
    "poisson"
)

supported_families <- c(
    "categorical",
    "gaussian"
)

# Generate is_x convenience functions for all supported families x
for (family in supported_families) {
    assign(paste0("is_", family), (function(y) {
        force(y)
        function(x) {
            identical(x$name, y)
        }
    })(family))
}

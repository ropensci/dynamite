# TODO implement families

#' Family Functions for \pkg{btvcm} Models
#'
#' @export
btvcmfamily <- function(name, ...) {
    # A wrapper to allow both:
    #    gaussian(...) and btvcmfamily("gaussian", ...)
    btvcmfamily_(name, ...)
}

# Internal use
btvcmfamily_ <- function(name, ...) {
    name <- tolower(as.character(name)[1])
    if (!is_supported(name)) {
        stop_(name, " is not a supported family")
    }
    # do something
    out <- list(name = name)
    class(out) <- c("btvcmfamily")
    out
}

# Checks if argument is a btvcmfamily object
is.btvcmfamily <- function(x) {
    inherits(x, "btvcmfamily")
}
# TODO could have an option to define the reference category? Or maybe simpler (for us)
# to ask user to relevel the factor beforehand
#' @rdname btvcmfamily
#' @export
categorical <- function(...) {
    # do something different
    btvcmfamily_("categorical", ...)
}

# TODO some families might conflict with stats such as stats::gaussian, not sure if matters
# Perhaps check how brms combines the use of own custom families and stats::gamma etc?
# What features we need for the family functions?
#' @rdname btvcmfamily
#' @export
gaussian <- function(...) {
    # do something else
    btvcmfamily_("gaussian", ...)
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
    "binomial",
    "categorical",
    "gaussian",
    "gamma",
    "poisson"
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
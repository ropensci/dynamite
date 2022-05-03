# TODO implement families

#' Family Functions for \pkg{dynamite} Models
#'
#' @export
dynamitefamily <- function(name, ...) {
    # A wrapper to allow both:
    #    gaussian(...) and dynamitefamily("gaussian", ...)
    dynamitefamily_(name, ...)
}

# Internal use
dynamitefamily_ <- function(name, ...) {
    name <- tolower(as.character(name)[1])
    if (!is_supported(name)) {
        stop_(name, " is not a supported family")
    }
    # do something
    out <- list(name = name)
    class(out) <- c("dynamitefamily")
    out
}

# Checks if argument is a dynamitefamily object
is.dynamitefamily <- function(x) {
    inherits(x, "dynamitefamily")
}

# TODO could have an option to define the reference category? Or maybe simpler (for us)
# to ask user to relevel the factor beforehand
categorical_ <- function(...) {
    # do something different
    dynamitefamily_("categorical", ...)
}

gaussian_ <- function(...) {
    # do something else
    dynamitefamily_("gaussian", ...)
}

binomial_ <- function(...) {
    # do something else
    dynamitefamily_("binomial", ...)
}

bernoulli_ <- function(...) {
    # do something else
    dynamitefamily_("bernoulli", ...)
}

poisson_ <- function(...) {
    # do something else
    dynamitefamily_("poisson", ...)
}

negbin_ <- function(...) {
    # do something else
    dynamitefamily_("negbin", ...)
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

# Check is a function call is of the form family(...),
# where 'family' is supported
is_valid_family_call <- function(x) {
    family <- as.character(x[[1]])[1]
    x[[1]] <- as.symbol(paste0(family, "_"))
    return(list(supported = is_supported(family), call = x))
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
    "bernoulli", # separate as stan has separate (more efficient) pmf for it
    "categorical",
    "negbin",
    "gaussian",
    "poisson" #TODO add gamma
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

# TODO implement families

#' Family Functions for \pkg{dynamite} Models
#'
#' @param \[`character(1)`]\cr Name of the family.
#' @param ... TODO is this needed?
#'
#' @export
dynamitefamily <- function(name) {
  # A wrapper to allow both:
  #    gaussian(...) and dynamitefamily("gaussian", ...)
  dynamitefamily_(name)
}

#' Construct a dynamitefamily object
#'
#' Checks that a given family is supported by the package
#'
#' @param name \[`character(1)`]\cr Name of the family.
#' @param ... TODO is this needed?
#'
#' @noRd
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

#' Is an object dynamitefamily object?
#'
#' Checks if argument is a dynamitefamily object
#'
#' @param x An R object.
#'
#' @noRd
is.dynamitefamily <- function(x) {
  inherits(x, "dynamitefamily")
}

# TODO could have an option to define the reference category?
# Or maybe simpler (for us) to ask user to relevel the factor beforehand
categorical_ <- function(...) {
  # do something different
  dynamitefamily_("categorical", ...)
}

#* Gaussian family
#*
#' @noRd
gaussian_ <- function(...) {
  # do something else
  dynamitefamily_("gaussian", ...)
}

#* Binomial family
#*
#' @noRd
binomial_ <- function(...) {
  # do something else
  dynamitefamily_("binomial", ...)
}

#* Bernoulli family
#*
#' @noRd
bernoulli_ <- function(...) {
  # do something else
  dynamitefamily_("bernoulli", ...)
}

#* Poisson family
#*
#' @noRd
poisson_ <- function(...) {
  # do something else
  dynamitefamily_("poisson", ...)
}

#* Negative binomial family
#*
#' @noRd
negbin_ <- function(...) {
  # do something else
  dynamitefamily_("negbin", ...)
}

#' Check if response of a family is continuous
#'
#' @param x A `dynamitefamily` object.
#'
#' @noRd
is_continuous <- function(x) {
  x$name %in% continuous_distributions
}

#' Check if response of a family is discrete
#'
#' @param x A `dynamitefamily` object.
#'
#' @noRd
is_discrete <- function(x) {
  x$name %in% discrete_distributions
}

#' Validate calls to family functions
#'
#' Checks if a function call is of the form family(...),
#' where 'family' is supported by the package.
#'
#' @param x A language object.
#'
#' @noRd
is_valid_family_call <- function(x) {
  family <- as.character(x[[1]])[1]
  x[[1]] <- as.symbol(paste0(family, "_"))
  return(list(supported = is_supported(family), call = x))
}

#' Check if a family is supported
#'
#' @param \[`character(1)`]\cr Name of the family.
#'
#' @noRd
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
  "poisson" # TODO add gamma
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

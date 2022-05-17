# TODO implement families

#' Family Functions for \pkg{dynamite} Models
#'
#' @param name \[`character(1)`]\cr Name of the family.
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
#' @noRd
dynamitefamily_ <- function(name) {
  name <- tolower(as.character(name)[1])
  if (!is_supported(name)) {
    stop_(name, " is not a supported family")
  }
  # do something
  structure(
    list(name = name),
    class = "dynamitefamily"
  )
}

#' Is an object a dynamitefamily object?
#'
#' Checks if argument is a dynamitefamily object
#'
#' @param x An R object.
#'
#' @noRd
is.dynamitefamily <- function(x) {
  inherits(x, "dynamitefamily")
}

#' Validate calls to family functions
#'
#' Checks if a function call is of the form family(...),
#' where 'family' is supported by the package.
#'
#' @param x A language object
#'
#' @noRd
is_valid_family_call <- function(x) {
  family <- as.character(x[[1]])[1]
  x[[1]] <- as.symbol(paste0(family, "_"))
  list(
    supported = is_supported(family),
    call = x,
    call_str = family
  )
}

#' Check if a family is supported
#'
#' @param \[`character(1)`]\cr Name of the family
#'
#' @noRd
is_supported <- function(name) {
  name %in% supported_families
}

# TODO add gamma
supported_families <- c(
  "binomial",
  "bernoulli", # separate as Stan has more efficient pmf for it
  "categorical",
  "negbin",
  "gaussian",
  "poisson",
  "deterministic"
)

# Generate 'family_' and 'is_family' convenience functions
# for all supported families
for (family in supported_families) {
  force(family)
  assign(paste0("is_", family), (function(y) {
    force(y)
    function(x) {
      identical(x$name, y)
    }
  })(family))
  assign(paste0(family, "_"), (function(y) {
    force(y)
    function() {
      dynamitefamily_(y)
    }
  })(family))
}

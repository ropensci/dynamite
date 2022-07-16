#' Construct a `dynamitefamily` Object
#'
#' Family Functions for \pkg{dynamite} Models
#'
#' @param name \[`character(1)`]\cr Name of the family.
#' @noRd
dynamitefamily <- function(name) {
  name <- tolower(as.character(name)[1])
  if (!is_supported(name)) {
    stop_("{.val {name}} is not a supported family.")
  }
  structure(
    list(name = name),
    class = "dynamitefamily"
  )
}

#' Is an Object a `dynamitefamily` Object?
#'
#' @param x An R object.
#' @noRd
is.dynamitefamily <- function(x) {
  inherits(x, "dynamitefamily")
}

#' Check if a Family is Supported
#'
#' @param \[`character(1)`]\cr Name of the family.
#' @noRd
is_supported <- function(name) {
  name %in% supported_families
}

supported_families <- c(
  "binomial",
  "bernoulli", # separate as Stan has more efficient pmf for it
  "categorical",
  "negbin",
  "gaussian",
  "poisson",
  "deterministic",
  "gamma",
  "exponential",
  "beta"
)

# Generate `family_` and `is_family` convenience functions
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
      dynamitefamily(y)
    }
  })(family))
}

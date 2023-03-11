#' Construct a `dynamitefamily` Object
#'
#' Family Functions for \pkg{dynamite} Models
#'
#' @param name \[`character(1)`]\cr Name of the family.
#' @noRd
dynamitefamily <- function(name) {
  name <- tolower(as.character(name)[1])
  stopifnot_(
    is_supported(name),
    "{.val {name}} is not a supported family."
  )
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

#' Check If a Family Is Supported
#'
#' @param name \[`character(1)`]\cr Name of the family.
#' @noRd
is_supported <- function(name) {
  name %in% supported_families
}

#' Check If a Family Is Multivariate
#'
#' @param x \[`dynamitefamily`]\cr A family object.
#' @noRd
is_multivariate <- function(x) {
  x$name %in% c("mvgaussian", "multinomial", "categorical")
}

supported_families <- c(
  "binomial",
  "bernoulli", # separate as Stan has more efficient pmf for it
  "categorical",
  "multinomial",
  "negbin",
  "gaussian",
  "mvgaussian",
  "poisson",
  "deterministic",
  "gamma",
  "exponential",
  "beta",
  "student"
)

#' Test If Multivariate Family Uses Univariate Components
#'
#' @param x \[`dynamitefamily`]\cr A family object.
#' @noRd
has_univariate <- function(x) {
  x$name %in% setdiff(supported_families, "multinomial")
}

#' Get Univariate Version of a Multivariate Family
#'
#' @param x \[`dynamitefamily`]\cr A family object.
#' @noRd
get_univariate <- function(x) {
  out <- stats::setNames(supported_families, supported_families)
  out["mvgaussian"] <- "gaussian"
  out["categorical"] <- "category"
  unname(out[x$name])
}

#' Is the GLM Likelihood Variant Supported By Stan for a Family
stan_supports_glm_likelihood <- function(x) {
  x$name %in% c("bernoulli", "gaussian", "poisson", "negbin")
}

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

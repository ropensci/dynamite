#' Construct a `dynamitefamily` Object
#'
#' Family Functions for \pkg{dynamite} Models
#'
#' @param name \[`character(1)`]\cr Name of the family.
#' @param link \[`character(1)`]\cr Name of the link function.
#' @noRd
dynamitefamily <- function(name, link) {
  name <- tolower(as.character(name)[1L])
  stopifnot_(
    is_supported_family(name),
    "{.val {name}} is not a supported family."
  )
  link <- ifelse_(
    missing(link) || is.null(link),
    default_link(name),
    tolower(as.character(link)[1L])
  )
  stopifnot_(
    is_supported_link(name, link),
    "{.val {link}} is not a supported link function for a {.val {name}} channel."
  )
  structure(
    list(name = name, link = link),
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
is_supported_family <- function(name) {
  name %in% supported_families
}

#' Check If a Link Function Is Supported for a Family
#'
#' @param name \[`character(1)`]\cr Name of the family.
#' @param link \[`character(1)`]\cr Name of the link function.
#' @noRd
is_supported_link <- function(name, link) {
  link %in% supported_links[[name]]
}

#' Check If a Family Is Multivariate
#'
#' @param x \[`dynamitefamily`]\cr A family object.
#' @noRd
is_multivariate <- function(x) {
  x$name %in% c("mvgaussian", "multinomial")
}

#' Get the Default Link of a Family
#'
#' @param name \[`character(1)`]\cr Name of the family.
#' @noRd
default_link <- function(name) {
  supported_links[[name]][1L]
}

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
  unname(out[x$name])
}

supported_families <- c(
  "bernoulli", # separate as Stan has more efficient pmf for it
  "beta",
  "binomial",
  "categorical",
  "cumulative", # ordered probit and logit
  "deterministic",
  "exponential",
  "gamma",
  "gaussian",
  "multinomial",
  "mvgaussian",
  "negbin",
  "poisson",
  "student"
)

supported_links <- list(
  bernoulli = c("logit"),
  beta = c("logit"),
  binomial = c("logit"),
  categorical = c("softmax"),
  cumulative = c("logit", "probit"),
  deterministic = c("identity"),
  exponential = c("log"),
  gamma = c("log"),
  gaussian = c("identity"),
  multinomial = c("softmax"),
  mvgaussian = c("identity"),
  negbin = c("log"),
  poisson = c("log"),
  student = c("identity")
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
    function(link) {
      dynamitefamily(y, link)
    }
  })(family))
}

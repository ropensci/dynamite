#' Create A Backward Compatible Stan Array Declaration
#'
#' @param backend \[`character(1)`]\cr Either `"rstan"` or `"cmdstanr"`
#' @param type \[`character(1)`]\cr Variable type, e.g., `"int"`, `"real"` etc.
#' @param name \[`character(1)`]\cr Name of the variable.
#' @param arr_dis \[`character(1)`]\cr Dimensions of the array
#'   (without brackets).
#' @param bounds \[`character(1)`]\cr Bounds of the variable (without < or >)
#' @param dims \[`character(1)`]\cr Dimensions of the array elements
#'   (without brackets).
#' @param comment \[`character(1)`]\cr Comment string to append to the end of
#'   the line (wihtout // prefix).
#' @noRd
stan_array <- function(backend, type, name, arr_dims,
                       bounds = "", dims = "", comment = "") {
  dims <- ifelse_(
    nzchar(dims),
    paste0("[", dims, "]"),
    ""
  )
  bounds <- ifelse_(
    nzchar(bounds),
    paste0("<", bounds, ">"),
    ""
  )
  comment <- ifelse_(
    nzchar(comment),
    paste0(" // ", comment),
    ""
  )
  ifelse_(
    stan_supports_array_keyword(backend),
    paste0(
      "array[", arr_dims, "] ", type, bounds, dims, " ", name, ";", comment
    ),
    paste0(
      type, dims, bounds, " ", name, "[", arr_dims, "];", comment
    )
  )
}

#' Is Array Keyword Syntax Supported By Current Stan Version
#'
#' @param backend Either `"rstan"` or `"cmdstanr"`.
#' @noRd
stan_supports_array_keyword <- function(backend) {
  utils::compareVersion(stan_version(backend), "2.26") >= 0
}

#' Is Categorical Logit GLM Supported By Current Stan Version
#'
#' @noRd
stan_supports_categorical_logit_glm <- function(backend) {
  utils::compareVersion(stan_version(backend), "2.23") >= 0
}

#' Get Stan Version
#'
#' @param backend Either `"rstan"` or `"cmdstanr"`.
#' @noRd
stan_version <- function(backend) {
  ifelse_(
    backend == "rstan",
    as.character(rstan::stan_version()),
    as.character(cmdstanr::cmdstan_version())
  )
}

#' Is the GLM Likelihood Variant Supported By Stan for a Family
#'
#' @param x \[`dynamitefamily`]\cr A family object.
#' @noRd
stan_supports_glm_likelihood <- function(x) {
  x$name %in% c("bernoulli", "gaussian", "poisson", "negbin")
}

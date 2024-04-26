#' Create A Backward Compatible Stan Array Declaration
#'
#' @param backend \[`character(1)`]\cr Either `"rstan"` or `"cmdstanr"`
#' @param type \[`character(1)`]\cr Variable type, e.g., `"int"`, `"real"` etc.
#' @param name \[`character(1)`]\cr Name of the variable.
#' @param arr_dims \[`character(1)`]\cr Dimensions of the array
#'   (without brackets).
#' @param bounds \[`character(1)`]\cr Bounds of the variable (without < or >)
#' @param dims \[`character(1)`]\cr Dimensions of the array elements
#'   (without brackets).
#' @param comment \[`character(1)`]\cr Comment string to append to the end of
#'   the line (without // prefix).
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

#' Create A Backward Compatible Stan Array for Function Arguments
#' @noRd
stan_array_arg <- function(backend, type, name, n_dims = 0, data = FALSE) {
  commas <- paste(rep(",", n_dims), collapse = " ")
  data <- ifelse(
    data,
    "data ",
    ""
  )
  ifelse_(
    stan_supports_array_keyword(backend),
    paste0(data, "array[", commas, "] ", type, " ", name),
    paste0(data, type, "[",commas, "] ", name)
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
stan_supports_categorical_logit_glm <- function(backend,
                                                common_intercept = TRUE) {
  utils::compareVersion(stan_version(backend), "2.23") >= 0 && common_intercept
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

#' Check That Stan Installation Is Functional
#'
#' @noRd
stan_version_is_functional <- function() {
  !is_windows() ||
    R_version() < "4.2.0" ||
    utils::compareVersion(stan_version("rstan"), "2.26") >= 0
}

#' Is the GLM Likelihood Variant Supported By Stan for a Family
#'
#' @param family \[`dynamitefamily`]\cr A family object.
#' @param backend Either `"rstan"` or `"cmdstanr"`.
#' @param common_intercept \[`logical(1)`]\cr Does the intercept vary between
#'   groups?
#' @noRd
stan_supports_glm_likelihood <- function(family, backend, common_intercept) {
  ifelse_(
    is_categorical(family),
    stan_supports_categorical_logit_glm(backend, common_intercept),
    (family$name %in% c("bernoulli", "gaussian", "poisson", "negbin")) ||
      (identical(family$name, "cumulative") && identical(family$link, "logit"))
  )
}

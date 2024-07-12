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
    paste0(data, type, "[", commas, "] ", name)
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


# Wrapper methods for backends --------------------------------------------

#' Get `pars_oi` of a Stan model fit
#'
#' @param x A `stanfit` (from `rstan`) or a `CmdStanMCMC`
#'   (from `cmdstanr`) object.
#' @noRd
get_pars_oi <- function(x) {
  UseMethod("get_pars_oi")
}

#' Get the model code of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @noRd
get_model_code <- function(x) {
  UseMethod("get_model_code")
}

#' Get the number of chains of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @noRd
get_nchains <- function(x) {
  UseMethod("get_nchains")
}

#' Get the algorithm used in a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @noRd
get_algorithm <- function(x) {
  UseMethod("get_algorithm")
}

#' Get the diagnostics of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @noRd
get_diagnostics <- function(x) {
  UseMethod("get_diagnostics")
}

#' Get the maximum treedepth of chains of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @noRd
get_max_treedepth <- function(x) {
  UseMethod("get_max_treedepth")
}

#' Get the number of draws of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @noRd
get_ndraws <- function(x) {
  UseMethod("get_ndraws")
}

#' Get the draws of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @noRd
get_draws <- function(x, ...) {
  UseMethod("get_draws")
}

#' Get the elapsed time of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @noRd
get_elapsed_time <- function(x) {
  UseMethod("get_elapsed_time")
}

get_pars_oi.stanfit <- function(x) {
  x@sim$pars_oi
}

get_pars_oi.CmdStanMCMC <- function(x) {
  x$metadata()$stan_variables
}

get_pars_oi.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_pars_oi")
}

get_model_code.stanfit <- function(x) {
  x@stanmodel@model_code[1L]
}

get_model_code.CmdStanMCMC <- function(x) {
  x$code()
}

get_model_code.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_model_code")
}

get_nchains.stanfit <- function(x) {
  x@sim$chains
}

get_nchains.CmdStanMCMC <- function(x) {
  x$num_chains()
}

get_nchains.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_nchains")
}

get_algorithm.stanfit <- function(x) {
  x@stan_args[[1L]]$algorithm
}

get_algorithm.CmdStanMCMC <- function(x) {
  x$metadata()$algorithm
}

get_algorithm.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_algorithm")
}

get_diagnostics.stanfit <- function(x) {
  list(
    num_divergent = rstan::get_num_divergent(x),
    num_max_treedepth = rstan::get_num_max_treedepth(x),
    ebfmi = rstan::get_bfmi(x)
  )
}

get_diagnostics.CmdStanMCMC <- function(x) {
  x$diagnostic_summary()
}

get_diagnostics.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_diagnostics")
}

get_max_treedepth.stanfit <- function(x) {
  x@stan_args[[1L]]$control$max_treedepth
}

get_max_treedepth.CmdStanMCMC <- function(x) {
  x$metadata()$max_treedepth
}

get_max_treedepth.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_max_treedepth")
}

get_ndraws.stanfit <- function(x) {
  (x@sim$n_save[1L] - x@sim$warmup2[1L]) * x@sim$chains
}

get_ndraws.CmdStanMCMC <- function(x) {
  m <- x$metadata()
  x$metadata()$iter_sampling * get_nchains(x)
}

get_ndraws.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_ndraws")
}

get_draws.stanfit <- function(x, pars) {
  posterior::as_draws(
    rstan::extract(
      x,
      pars = pars,
      permuted = FALSE
    )
  )
}

get_draws.CmdStanMCMC <- function(x, pars) {
  x$draws(variables = pars)
}

get_draws.CMdStanMCMC_CSV <- function(x) {
  NextMethod("get_draws")
}

get_elapsed_time.stanfit <- function(x) {
  rstan::get_elapsed_time(x)
}

get_elapsed_time.CmdStanMCMC <- function(x) {
  x$time()$chains
}

get_elapsed_time.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_elapsed_time")
}

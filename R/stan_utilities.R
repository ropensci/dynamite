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

stan_reserved_keywords <- c(
  "int", "real", "vector", "row_vector", "matrix", "ordered",
  "positive_ordered", "simplex", "unit_vector", "cholesky_factor_corr",
  "cholesky_factor_cov", "corr_matrix", "cov_matrix", "functions", "model",
  "parameters", "transformed", "generated", "quantities", "data", "var",
  "return", "if", "else", "while", "for", "in", "break", "continue", "void",
  "reject", "print", "target", "T"
)

#' Ensure that a character string is a valid Stan variable name
#'
#' This function prepares a name such that it is valid for Stan. From
#' Stan Reference Manual: "A variable by itself is a well-formed expression of
#' the same type as the variable. Variables in Stan consist of ASCII strings
#' containing only the basic lower-case and upper-case Roman letters, digits,
#' and the underscore (_) character. Variables must start with a letter
#' (a--z and A--Z) and may not end with two underscores (__)". Adds a prefix
#' when the first character is not a letter and a suffix for reserved keywords.
#'
#' @param x A `character` vector.
#' @noRd
stan_name <- function(x, check_first = TRUE) {
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^a-zA-Z0-9_]", "", x)
  x <- gsub("_{2,}$", "", x)
  for (i in seq_along(x)) {
    if (check_first && !grepl("^[a-zA-Z]", x[i])) {
      x[i] <- paste0("v_", x[i])
    }
    if (tolower(x[i]) %in% stan_reserved_keywords) {
      x[i] <- paste0(x[i], "_var")
    }
  }
  x
}

# Wrapper methods for backends --------------------------------------------

#' Get `pars_oi` of a Stan model fit
#'
#' @param x A `stanfit` (from `rstan`) or a `CmdStanMCMC`
#'   (from `cmdstanr`) object.
#' @keywords internal
#' @export
get_pars_oi <- function(x) {
  UseMethod("get_pars_oi")
}

#' Get the model code of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @keywords internal
#' @export
get_model_code <- function(x) {
  UseMethod("get_model_code")
}

#' Get the number of chains of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @keywords internal
#' @export
get_nchains <- function(x) {
  UseMethod("get_nchains")
}

#' Get the algorithm used in a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @keywords internal
#' @export
get_algorithm <- function(x) {
  UseMethod("get_algorithm")
}

#' Get the diagnostics of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @keywords internal
#' @export
get_diagnostics <- function(x) {
  UseMethod("get_diagnostics")
}

#' Get the maximum treedepth of chains of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @keywords internal
#' @export
get_max_treedepth <- function(x) {
  UseMethod("get_max_treedepth")
}

#' Get the number of draws of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @keywords internal
#' @export
get_ndraws <- function(x) {
  UseMethod("get_ndraws")
}

#' Get the draws of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @param pars A `character` vector of parameter names.
#' @keywords internal
#' @export
get_draws <- function(x, pars) {
  UseMethod("get_draws")
}

#' Get the elapsed time of a Stan model fit
#'
#' @inheritParams get_pars_oi
#' @keywords internal
#' @export
get_elapsed_time <- function(x) {
  UseMethod("get_elapsed_time")
}

#' @rdname get_pars_oi
#' @keywords internal
#' @export
get_pars_oi.stanfit <- function(x) {
  x@sim$pars_oi
}

#' @rdname get_pars_oi
#' @keywords internal
#' @export
get_pars_oi.CmdStanMCMC <- function(x) {
  x$metadata()$stan_variables
}

#' @rdname get_pars_oi
#' @keywords internal
#' @export
get_pars_oi.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_pars_oi")
}

#' @rdname get_model_code
#' @keywords internal
#' @export
get_model_code.stanfit <- function(x) {
  x@stanmodel@model_code[1L]
}

#' @rdname get_model_code
#' @keywords internal
#' @export
get_model_code.CmdStanMCMC <- function(x) {
  paste0(x$code(), collapse = "\n")
}

#' @rdname get_model_code
#' @keywords internal
#' @export
get_model_code.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_model_code")
}

#' @rdname get_nchains
#' @keywords internal
#' @export
get_nchains.stanfit <- function(x) {
  x@sim$chains
}

#' @rdname get_nchains
#' @keywords internal
#' @export
get_nchains.CmdStanMCMC <- function(x) {
  x$num_chains()
}

#' @rdname get_nchains
#' @keywords internal
#' @export
get_nchains.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_nchains")
}

#' @rdname get_algorithm
#' @keywords internal
#' @export
get_algorithm.stanfit <- function(x) {
  x@stan_args[[1L]]$algorithm
}

#' @rdname get_algorithm
#' @keywords internal
#' @export
get_algorithm.CmdStanMCMC <- function(x) {
  x$metadata()$algorithm
}

#' @rdname get_algorithm
#' @keywords internal
#' @export
get_algorithm.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_algorithm")
}

#' @rdname get_diagnostics
#' @keywords internal
#' @export
get_diagnostics.stanfit <- function(x) {
  list(
    num_divergent = rstan::get_num_divergent(x),
    num_max_treedepth = rstan::get_num_max_treedepth(x),
    ebfmi = rstan::get_bfmi(x)
  )
}

#' @rdname get_diagnostics
#' @keywords internal
#' @export
get_diagnostics.CmdStanMCMC <- function(x) {
  x$diagnostic_summary()
}

#' @rdname get_diagnostics
#' @keywords internal
#' @export
get_diagnostics.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_diagnostics")
}

#' @rdname get_max_treedepth
#' @keywords internal
#' @export
get_max_treedepth.stanfit <- function(x) {
  x@stan_args[[1L]]$control$max_treedepth
}

#' @rdname get_max_treedepth
#' @keywords internal
#' @export
get_max_treedepth.CmdStanMCMC <- function(x) {
  x$metadata()$max_treedepth
}

#' @rdname get_max_treedepth
#' @keywords internal
#' @export
get_max_treedepth.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_max_treedepth")
}

#' @rdname get_ndraws
#' @keywords internal
#' @export
get_ndraws.stanfit <- function(x) {
  (x@sim$n_save[1L] - x@sim$warmup2[1L]) * x@sim$chains
}

#' @rdname get_ndraws
#' @keywords internal
#' @export
get_ndraws.CmdStanMCMC <- function(x) {
  m <- x$metadata()
  x$metadata()$iter_sampling * get_nchains(x)
}

#' @rdname get_ndraws
#' @keywords internal
#' @export
get_ndraws.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_ndraws")
}

#' @rdname get_draws
#' @keywords internal
#' @export
get_draws.stanfit <- function(x, pars) {
  posterior::as_draws(
    rstan::extract(
      x,
      pars = pars,
      permuted = FALSE
    )
  )
}

#' @rdname get_draws
#' @keywords internal
#' @export
get_draws.CmdStanMCMC <- function(x, pars) {
  x$draws(variables = pars)
}

#' @rdname get_draws
#' @keywords internal
#' @export
get_draws.CMdStanMCMC_CSV <- function(x, pars) {
  NextMethod("get_draws")
}

#' @rdname get_elapsed_time
#' @keywords internal
#' @export
get_elapsed_time.stanfit <- function(x) {
  rstan::get_elapsed_time(x)
}

#' @rdname get_elapsed_time
#' @keywords internal
#' @export
get_elapsed_time.CmdStanMCMC <- function(x) {
  x$time()$chains
}

#' @rdname get_elapsed_time
#' @keywords internal
#' @export
get_elapsed_time.CmdStanMCMC_CSV <- function(x) {
  NextMethod("get_elapsed_time")
}

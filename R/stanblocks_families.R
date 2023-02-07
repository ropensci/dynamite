#' Wrapper to Parse Stan Model Blocks Based on the Family of the Response
#'
#' @param prefix \[`character(1)`]\cr Stan model block name, e.g., "model".
#' @param family \[`character(1)`]\cr Supported family name
#' @param args Channel specific component of `model_vars`
#'   (see [create_blocks()])
#' @noRd
lines_wrap <- function(prefix, family, args) {
  do.call(what = paste0(prefix, "_lines_", family), args = args)
}

#' Is a Prior Definition Vectorizable
#'
#' @param x A `character` vector of length one.
#' @noRd
vectorizable_prior <- function(x) {
  length(x) == 1L && !grepl("\\(", x)
}

#' Parameters for Stan Code Generation
#'
#' Each of the following functions of the form `(stan_block)_lines_(family)`
#' uses a subset from a  common set of arguments that are documented here.
#'
#' @param backend \[`character(1)`]\cr `"rstan"` or `"cmdstanr"`.
#' @param y \[`character(1)`]\cr The name of the response of the channel.
#' @param idt \[`function`]\cr An indentation function, see [indenter_()].
#' @param obs \[`integer()`]\cr A vector of indices indicating non-missing
#'   values of the channel.
#' @param has_missing \[`logical(1)`]\cr Does the channel contain missing
#'   observations?
#' @param has_offset \[`logical(1)`]\cr Does the channel contain an offset term
#' @param has_fixed \[`logical(1)`]\cr Does the channel have time-invariant
#'   predictors?
#' @param has_varying \[`logical(1)`]\cr Does the channel have time-varying
#'   predictors?
#' @param has_random \[`logical(1)`]\cr  Does the channel have random
#' (group-specific) effects?
#' @param has_fixed_intercept \[`logical(1)`]\cr Does the channel have a
#'   time-invariant intercept?
#' @param has_varying_intercept \[`logical(1)`]\cr Does the channel have a
#'   time-varying intercept?
#' @param has_random_intercept \[`logical(1)`]\cr Does the channel have a
#'   random intercept?
#' @param noncentered \[`logical(1)`]\cr Should the noncentered parametrization
#'   be used for splines?
#' @param shrinkage \[`logical(1)`]\cr Should the common global shrinkage
#'   parameter be used?
#' @param lb \[`double(1)`]\cr Lower bound for the `tau` parameter.
#' @param J \[`integer()`]\cr Model matrix column indices of the predictors
#'   of the channel
#' @param J_fixed \[`integer()`]\cr Model matrix column indices of the
#'   time-invariant predictors of the channel
#' @param J_varying \[`integer()`]\cr Model matrix column indices of the
#'   time-varying predictors of the channel
#' @param K \[`integer(1)`]\cr Total number of predictors of the channel
#' @param K_fixed \[`integer(1)`]\cr Number of time-invariant predictors of
#'   the channel.
#' @param K_varying \[`integer(1)`]\cr Number of time-varying predictors of
#'   the channel.
#' @param K_random \[`integer(1)`]\cr Number of random effect predictors of
#'   the channel.
#' @param L_fixed \[`integer(1)`]\cr Indices of the time-invariant predictors
#'   of the channel in the joint parameter vector `gamma`.
#' @param L_varying \[`integer(1)`]\cr Indices of the time-varying predictors
#'   of the channel in the joint parameter vector `gamma`.
#' @param S \[`integer(1)`]\cr Number of categories for a categorical channel.
#' @param write_alpha \[`logical(1)`]\cr If `TRUE`, use vectorized prior for
#'   `alpha` parameters of the model by hard-coding the parameters of the prior
#'   distribution to the model code?
#' @param write_beta \[`logical(1)`]\cr If `TRUE`, use vectorized prior for
#'   time-invariant predictors `beta`.
#' @param write_delta \[`logical(1)`]\cr If `TRUE`, use vectorized prior for
#'   time-varying predictors `delta`.
#' @param write_tau \[`logical(1)`]\cr If `TRUE`, use vectorized prior for
#'   `tau` parameters.
#' @param write_sigma_nu \[`logical(1)`]\cr If `TRUE`, use vectorized prior for
#'   `sigma_nu` parameters.
#' @param sigma_prior_distr \[`character(1)`]\cr `sigma` parameter prior
#'   specification.
#' @param sigma_nu_prior_distr \[`character(1)`]\cr `sigma_nu` parameter prior
#'   specification.
#' @param alpha_prior_distr \[`character(1)`]\cr `alpha` parameter prior
#'   specification.
#' @param alpha_prior_npars \[`integer(1)`]\cr Number of parameters for the
#'   prior of `alpha`.
#' @param tau_alpha_prior_distr \[`character(1)`]\cr `tau_alpha` parameter
#'   prior specification.
#' @param beta_prior_distr \[`character(1)`]\cr `beta` parameter prior
#'   specification.
#' @param beta_prior_npars \[`integer(1)`]\cr Number of parameters for the
#'   prior of `beta`.
#' @param delta_prior_distr \[`character(1)`]\cr `delta` parameter prior
#'   specification.
#' @param delta_prior_npars \[`integer(1)`]\cr Number of parameters for the
#'   prior of `delta`.
#' @param tau_prior_distr \[`character(1)`]\cr `tau` parameter prior
#'   specification.
#' @param tau_prior_npars \[`integer(1)`]\cr Number of parameters for the
#'   prior of `tau`.
#' @param sigma_nu_prior_distr \[`character(1)`]\cr `sigma_nu` parameter prior
#'   specification.
#' @param sigma_nu_prior_npars \[`integer(1)`]\cr Number of parameters for the
#'   prior of `tau`.
#' @param phi_prior_distr \[`character(1)`]\cr `phi` parameter prior
#'   specification.
#' @noRd
NULL

# Data block --------------------------------------------------------------

data_lines_default <- function(y, idt, has_missing,
                               write_beta, write_delta, write_tau,
                               write_sigma_nu,
                               beta_prior_npars = 1L, delta_prior_npars = 1L,
                               tau_prior_npars = 1L, sigma_nu_prior_npars = 1L,
                               has_random_intercept, ...) {
  icpt <- ifelse(
    has_random_intercept,
    " - 1",
    ""
  )
  paste_rows(
    onlyif(has_missing, "// Missing data indicators"),
    onlyif(has_missing, "int<lower=0> obs_{y}[N, T];"),
    onlyif(has_missing, "int<lower=0> n_obs_{y}[T];"),
    "// number of fixed, varying and random coefficients, and related indices",
    "int<lower=0> K_fixed_{y};",
    "int<lower=0> K_varying_{y};",
    "int<lower=0> K_random_{y}; // includes the potential random intercept",
    "int<lower=0> K_{y}; // K_fixed + K_varying",
    "int J_fixed_{y}[K_fixed_{y}];",
    "int J_varying_{y}[K_varying_{y}];",
    "int J_random_{y}[K_random_{y}{icpt}]; // no intercept",
    "int J_{y}[K_fixed_{y} + K_varying_{y}]; // fixed and varying",
    "int L_fixed_{y}[K_fixed_{y}];",
    "int L_varying_{y}[K_varying_{y}];",
    "int L_random_{y}[K_random_{y}{icpt}];",
    onlyif(
      write_beta || write_delta || write_tau || write_sigma_nu,
      "// Parameters of vectorized priors"
    ),
    onlyif(
      write_beta && beta_prior_npars > 0L,
        "matrix[K_fixed_{y}, {beta_prior_npars}] beta_prior_pars_{y};"
    ),
    onlyif(
      write_delta && delta_prior_npars > 0L,
      "matrix[K_varying_{y}, {delta_prior_npars}] delta_prior_pars_{y};"
    ),
    onlyif(
      write_tau && tau_prior_npars > 0L,
      "matrix[K_varying_{y}, {tau_prior_npars}] tau_prior_pars_{y};"
    ),
    onlyif(
      write_sigma_nu && sigma_nu_prior_npars > 0L,
      "matrix[K_random_{y}, {sigma_nu_prior_npars}] sigma_nu_prior_pars_{y};"
    ),
    "// Responses",
    .indent = idt(1)
  )
}

data_lines_categorical <- function(y, idt,
                                   has_missing, write_alpha, write_beta,
                                   write_delta, write_tau,
                                   alpha_prior_npars = 1L,
                                   beta_prior_npars = 1L,
                                   delta_prior_npars = 1L, tau_prior_npars = 1L,
                                   has_random_intercept, ...) {
  icpt <- ifelse(
    has_random_intercept,
    " - 1",
    ""
  )
  dtext_def <- paste_rows(
    "int<lower=0> S_{y}; // number of categories",
    onlyif(has_missing, "// Missing data indicators"),
    onlyif(has_missing, "int<lower=0> obs_{y}[N, T];"),
    onlyif(has_missing, "int<lower=0> n_obs_{y}[T];"),
    "// number of fixed, varying and random coefficients, and related indices",
    "int<lower=0> K_fixed_{y};",
    "int<lower=0> K_varying_{y};",
    "int<lower=0> K_random_{y}; // includes the potential random intercept",
    "int<lower=0> K_{y}; // K_fixed + K_varying",
    "int J_fixed_{y}[K_fixed_{y}];",
    "int J_varying_{y}[K_varying_{y}];",
    "int J_random_{y}[K_random_{y}{icpt}]; // no intercept",
    "int J_{y}[K_fixed_{y} + K_varying_{y}]; // fixed and varying",
    "int L_fixed_{y}[K_fixed_{y}];",
    "int L_varying_{y}[K_varying_{y}];",
    "int L_random_{y}[K_random_{y}{icpt}];",
    onlyif(
      write_alpha || write_beta || write_delta || write_tau,
      "// Parameters of vectorized priors"
    ),
    onlyif(
      write_alpha && alpha_prior_npars > 0L,
      "matrix[S_{y} - 1, {alpha_prior_npars}] alpha_prior_pars_{y};"
    ),
    onlyif(
      write_beta && beta_prior_npars > 0L,
      "matrix[K_fixed_{y} * (S_{y} - 1), {beta_prior_npars}] beta_prior_pars_{y};"
    ),
    onlyif(
      write_delta && delta_prior_npars > 0L,
      "matrix[K_varying_{y} * (S_{y} - 1), {delta_prior_npars}] delta_prior_pars_{y};"
    ),
    onlyif(
      write_tau && tau_prior_npars > 0L,
      "matrix[K_varying_{y}, {tau_prior_npars}] tau_prior_pars_{y};"
    ),
    "// Responses",
    "int<lower=0> y_{y}[T, N];",
    .indent = idt(1)
  )
  paste_rows(dtext_def, .parse = FALSE)
}

data_lines_gaussian <- function(y, idt, has_missing, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(data_lines_default))]
  paste_rows(
    do.call(what = data_lines_default, args = args),
    "matrix[N, T] y_{y};",
    .indent = idt(c(0, 1))
  )
}

data_lines_binomial <- function(y, idt, has_missing, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(data_lines_default))]
  paste_rows(
    do.call(what = data_lines_default, args = args),
    "int<lower=0> y_{y}[T, N];",
    "// Trials for binomial response {y}",
    "int<lower=1> trials_{y}[T, N];",
    .indent = idt(c(0, 1, 1, 1))
  )
}

data_lines_bernoulli <- function(y, idt, has_missing, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(data_lines_default))]
  paste_rows(
    do.call(what = data_lines_default, args = args),
    "int<lower=0,upper=1> y_{y}[T, N];",
    .indent = idt(c(0, 1))
  )
}

data_lines_poisson <- function(y, idt, has_missing, has_offset, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(data_lines_default))]
  paste_rows(
    do.call(what = data_lines_default, args = args),
    "int<lower=0> y_{y}[T, N];",
    "// Offset term",
    onlyif(has_offset, "real offset_{y}[T, N];"),
    .indent = idt(c(0, 1, 1, 1))
  )
}

data_lines_negbin <- function(y, idt, has_missing, has_offset, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(data_lines_default))]
  paste_rows(
    do.call(what = data_lines_default, args = args),
    "int<lower=0> y_{y}[T, N];",
    "// Offset term",
    onlyif(has_offset, "real offset_{y}[T, N];"),
    .indent = idt(c(0, 1))
  )
}

data_lines_exponential <- function(y, idt, has_missing, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(data_lines_default))]
  paste_rows(
    do.call(what = data_lines_default, args = args),
    "matrix<lower=0>[N, T] y_{y};",
    .indent = idt(c(0, 1))
  )
}

data_lines_gamma <- function(y, idt, has_missing, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(data_lines_default))]
  paste_rows(
    do.call(what = data_lines_default, args = args),
    "matrix<lower=0>[N, T] y_{y};",
    .indent = idt(c(0, 1))
  )
}

data_lines_beta <- function(y, idt, has_missing, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(data_lines_default))]
  paste_rows(
    do.call(what = data_lines_default, args = args),
    "matrix<lower=0, upper=1>[N, T] y_{y};",
    .indent = idt(c(0, 1))
  )
}

# Transformed data block --------------------------------------------------

transformed_data_lines_default <- function(y, idt) {
  list(
    declarations = "",
    statements = ""
  )
}

transformed_data_lines_categorical <- function(y, idt, write_alpha, write_beta,
                                               write_delta, write_tau,
                                               K, K_fixed, K_varying, S,
                                               alpha_prior_npars = 1L,
                                               alpha_prior_pars = "",
                                               beta_prior_npars = 1L,
                                               beta_prior_pars = "",
                                               delta_prior_npars = 1L,
                                               delta_prior_pars = "",
                                               tau_prior_npars = 1L,
                                               tau_prior_pars = "", ...) {
  mtext <- paste_rows(
    "vector[{K}] zeros_K_{y} = rep_vector(0, {K});",
    "vector[{S}] zeros_S_{y} = rep_vector(0, {S});",
    .indent = idt(1)
  )
  list(
    declarations = mtext,
    statements = ""
  )
}

transformed_data_lines_gaussian <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(transformed_data_lines_default))]
  do.call(what = transformed_data_lines_default, args = args)
}

transformed_data_lines_binomial <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(transformed_data_lines_default))]
  do.call(what = transformed_data_lines_default, args = args)
}

transformed_data_lines_bernoulli <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(transformed_data_lines_default))]
  do.call(what = transformed_data_lines_default, args = args)
}

transformed_data_lines_poisson <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(transformed_data_lines_default))]
  do.call(what = transformed_data_lines_default, args = args)
}

transformed_data_lines_negbin <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(transformed_data_lines_default))]
  do.call(what = transformed_data_lines_default, args = args)
}

transformed_data_lines_exponential <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(transformed_data_lines_default))]
  do.call(what = transformed_data_lines_default, args = args)
}

transformed_data_lines_gamma <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(transformed_data_lines_default))]
  do.call(what = transformed_data_lines_default, args = args)
}

transformed_data_lines_beta <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(transformed_data_lines_default))]
  do.call(what = transformed_data_lines_default, args = args)
}

# Parameters block --------------------------------------------------------

parameters_lines_default <- function(y, idt, noncentered, lb, has_fixed,
                                     has_varying, has_fixed_intercept,
                                     has_varying_intercept,
                                     has_lfactor, noncentered_psi,
                                     noncentered_lambda,
                                     nonzero_lambda) {
  oname <- ifelse_(noncentered, "omega_raw_", "omega_")
  allow_intercept <- !(has_lfactor && nonzero_lambda)
  if (!allow_intercept && (has_fixed_intercept || has_varying_intercept)) {
    warning_("Separate common intercept term of channel {y} was removed as
      channel predictors contain possibly nonzero latent factor.")
  }
  paste_rows(
    onlyif(has_fixed, "vector[K_fixed_{y}] beta_{y}; // Fixed coefficients"),
    onlyif(
      has_varying,
      "matrix[K_varying_{y}, D] {oname}{y}; // Spline coefficients"
    ),
    onlyif(
      has_varying,
      "vector<lower={lb}>[K_varying_{y}] tau_{y}; // SDs for the random walks"
    ),
    ifelse_(
      (has_fixed_intercept || has_varying_intercept) && allow_intercept,
      "real a_{y}; // Mean of the first time point",
      ""
    ),
    onlyif(
      has_varying_intercept && allow_intercept,
      "row_vector[D - 1] omega_raw_alpha_{y}; // Coefficients for alpha"
    ),
    onlyif(
      has_varying_intercept && allow_intercept,
      "real<lower={lb}> tau_alpha_{y}; // SD for the random walk"
    ),
    onlyif(
      has_lfactor && nonzero_lambda,
      "real<lower=0> tau_psi_{y}; // SD for for the random walk"
    ),
    onlyif(
      has_lfactor,
      "real<lower=0> sigma_lambda_{y}; // SD of factor loadings"
    ),
    onlyif(
      has_lfactor,
      "vector[N - 1] lambda_raw_{y}; // raw factor loadings"
    ),
    onlyif(
      has_lfactor && noncentered_psi,
      "real omega_raw_psi_1_{y}; // factor spline coef for first time point"
    ),
    .indent = idt(1)
  )
}

parameters_lines_categorical <- function(y, idt, noncentered, lb, has_fixed,
                                         has_varying, has_fixed_intercept,
                                         has_varying_intercept, ...) {
  oname <- ifelse_(noncentered, "omega_raw_", "omega_")
  paste_rows(
    onlyif(
      has_fixed,
      "matrix[K_fixed_{y}, S_{y} - 1] beta_{y}; // Fixed coefficients"
    ),
    onlyif(
      has_varying,
      "matrix[K_varying_{y}, D] {oname}{y}[S_{y} - 1]; // Spline coefficients"
    ),
    onlyif(
      has_varying,
      "vector<lower={lb}>[K_varying_{y}] tau_{y};  // SDs for the random walks"
    ),
    ifelse_(
      has_fixed_intercept || has_varying_intercept,
      "vector[S_{y} - 1] a_{y}; // Mean of the first time point",
      ""
    ),
    onlyif(
      has_varying_intercept,
      paste0(
        "row_vector[D - 1] omega_raw_alpha_{y}[S_{y} - 1]; ",
        "// Coefficients for alpha"
      )
    ),
    onlyif(
      has_varying_intercept,
      "real<lower={lb}> tau_alpha_{y}; // SD for the random walk"
    ),
    .indent = idt(1)
  )
}

parameters_lines_gaussian <- function(y, idt, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(parameters_lines_default))]
  paste_rows(
    do.call(what = parameters_lines_default, args = args),
    "real<lower=0> sigma_{y}; // SD of the normal distribution",
    .indent = idt(c(0, 1))
  )
}

parameters_lines_binomial <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(parameters_lines_default))]
  do.call(what = parameters_lines_default, args = args)
}

parameters_lines_bernoulli <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(parameters_lines_default))]
  do.call(what = parameters_lines_default, args = args)
}

parameters_lines_poisson <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(parameters_lines_default))]
  do.call(what = parameters_lines_default, args = args)
}

parameters_lines_negbin <- function(y, idt, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(parameters_lines_default))]
  paste_rows(
    do.call(what = parameters_lines_default, args = args),
    "real<lower=0> phi_{y}; // Dispersion parameter of the NB distribution",
    .indent = idt(c(0, 1))
  )
}

parameters_lines_exponential <- function(y, idt, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(parameters_lines_default))]
  do.call(what = parameters_lines_default, args = args)
}

parameters_lines_gamma <- function(y, idt, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(parameters_lines_default))]
  paste_rows(
    do.call(what = parameters_lines_default, args = args),
    "real<lower=0> phi_{y}; // Shape parameter of the Gamma distribution",
    .indent = idt(c(0, 1))
  )
}

parameters_lines_beta <- function(y, idt, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(parameters_lines_default))]
  paste_rows(
    do.call(what = parameters_lines_default, args = args),
    "real<lower=0> phi_{y}; // Precision parameter of the Beta distribution",
    .indent = idt(c(0, 1))
  )
}

# Transformed parameters block --------------------------------------------

transformed_parameters_lines_default <- function(y, idt, noncentered,
                                                 shrinkage,
                                                 has_fixed, has_varying,
                                                 has_fixed_intercept,
                                                 has_varying_intercept,
                                                 J, K, K_fixed, K_varying,
                                                 L_fixed, L_varying,
                                                 has_lfactor,
                                                 noncentered_psi,
                                                 noncentered_lambda,
                                                 nonzero_lambda) {
  if (noncentered) {
    xi_term <- ifelse_(shrinkage, " * xi[i - 1];", ";")
    declare_omega <- paste_rows(
      "matrix[K_varying_{y}, D] omega_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega <- paste_rows(
      "omega_{y}[, 1] = omega_raw_{y}[, 1];",
      "for (i in 2:D) {{",
      paste0(
        "omega_{y}[, i] = omega_{y}[, i - 1] + ",
        "omega_raw_{y}[, i] .* tau_{y}{xi_term}"
      ),
      "}}",
      .indent = idt(c(1, 1, 2, 1)),
      .parse = FALSE
    )
    declare_omega_alpha <- paste_rows(
      "row_vector[D] omega_alpha_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha <- paste_rows(
      "omega_alpha_{y}[1] = omega_alpha_1_{y};",
      "for (i in 2:D) {{",
      paste0(
        "omega_alpha_{y}[i] = omega_alpha_{y}[i - 1] + ",
        "omega_raw_alpha_{y}[i - 1] * tau_alpha_{y}{xi_term}"
      ),
      "}}",
      .indent = idt(c(1, 1, 2, 1)),
      .parse = FALSE
    )
  } else {
    declare_omega_alpha <- paste_rows(
      "row_vector[D] omega_alpha_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha <- paste_rows(
      "omega_alpha_{y}[1] = omega_alpha_1_{y};",
      "omega_alpha_{y}[2:D] = omega_raw_alpha_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
  }

  declare_delta <- paste_rows(
    "// Time-varying coefficients",
    "vector[K_varying_{y}] delta_{y}[T];",
    .indent = idt(1),
    .parse = FALSE
  )
  state_delta <- paste_rows(
    "for (t in 1:T) {{",
    "delta_{y}[t] = omega_{y} * Bs[, t];",
    "}}",
    .indent = idt(c(1, 2, 1)),
    .parse = FALSE
  )

  if (has_fixed || has_varying) {
    declare_omega_alpha_1 <- paste_rows(
      "// Time-varying intercept",
      "real alpha_{y}[T];",
      "// Spline coefficients",
      "real omega_alpha_1_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha_1 <- paste_rows(
      "// Define the first alpha using mean a_{y}",
      "{{",
      "vector[K_{y}] gamma__{y};",
      onlyif(has_fixed, "gamma__{y}[L_fixed_{y}] = beta_{y};"),
      onlyif(has_varying, "gamma__{y}[L_varying_{y}] = delta_{y}[1];"),
      "omega_alpha_1_{y} = a_{y} - X_m[J_{y}] * gamma__{y};",
      "}}",
      .indent = idt(c(1, 1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
    declare_fixed_intercept <- paste_rows(
      "// Time-invariant intercept",
      "real alpha_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_fixed_intercept <- paste_rows(
      "// Define the first alpha using mean a_{y}",
      "{{",
      "vector[K_{y}] gamma__{y};",
      onlyif(has_fixed, "gamma__{y}[L_fixed_{y}] = beta_{y};"),
      onlyif(has_varying, "gamma__{y}[L_varying_{y}] = delta_{y}[1];"),
      "alpha_{y} = a_{y} - X_m[J_{y}] * gamma__{y};",
      "}}",
      .indent = idt(c(1, 1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
  } else {
    declare_omega_alpha_1 <- paste_rows(
      "// Time-invariant intercept",
      "real alpha_{y}[T];",
      "real omega_alpha_1_{y} = a_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha_1 <- character(0L)
    declare_fixed_intercept <- "real alpha_{y} = a_{y};"
    state_fixed_intercept <- character(0L)
  }
  declare_varying_intercept <- paste_rows(
    declare_omega_alpha_1,
    declare_omega_alpha,
    .indent = idt(0),
    .parse = FALSE
  )
  state_varying_intercept <- paste_rows(
    state_omega_alpha_1,
    state_omega_alpha,
    "for (t in 1:T) {{",
    "alpha_{y}[t] = omega_alpha_{y} * Bs[, t];",
    "}}",
    .indent = idt(c(0, 0, 1, 2, 1)),
    .parse = FALSE
  )

  m <- ifelse(nonzero_lambda, "1 + ", "")
  if (noncentered_lambda) {
    declare_lambda <- paste_rows(
      "// hard sum-to-zero constraint",
      "vector[N] lambda_std_{y} = A_qr * lambda_raw_{y};",
      "vector[N] lambda_{y} = {m}sigma_lambda_{y} * lambda_std_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
  } else {
    declare_lambda <- paste_rows(
      "// hard sum-to-zero constraint",
      "vector[N] lambda_{y} = {m}sigma_lambda_{y} * A_qr * lambda_raw_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
  }
  if (noncentered_psi) {
    state_omega_psi <- paste_rows(
      "omega_psi_{y}[1] = omega_raw_psi_1_{y};",
      "omega_psi_{y} = cumulative_sum(omega_psi_{y});",
      .indent = idt(1),
      .parse = FALSE
    )
  }
  declare_psi <- paste_rows(
    "// Latent factor",
    "vector[T] psi_{y};",
    .indent = idt(1),
    .parse = FALSE
  )
  state_psi <- paste_rows(
    "for (t in 1:T) {{",
    "psi_{y}[t] = omega_psi_{y} * Bs[, t];",
    "}}",
    .indent = idt(c(1, 2, 1)),
    .parse = FALSE
  )
  allow_intercept <- !(has_lfactor && nonzero_lambda)
  list(
    declarations = paste_rows(
      onlyif(has_varying && noncentered, declare_omega),
      onlyif(has_varying, declare_delta),
      onlyif(has_fixed_intercept && allow_intercept, declare_fixed_intercept),
      onlyif(has_varying_intercept && allow_intercept, declare_varying_intercept),
      onlyif(has_lfactor, declare_psi),
      onlyif(has_lfactor, declare_lambda),
      .indent = idt(0)
    ),
    statements = paste_rows(
      onlyif(has_varying && noncentered, state_omega),
      onlyif(has_varying, state_delta),
      onlyif(has_fixed_intercept && allow_intercept, state_fixed_intercept),
      onlyif(has_varying_intercept && allow_intercept, state_varying_intercept),
      onlyif(has_lfactor && noncentered_psi, state_omega_psi),
      onlyif(has_lfactor, state_psi),
      .indent = idt(0)
    )
  )
}

transformed_parameters_lines_categorical <- function(y, idt, noncentered,
                                                     shrinkage, has_fixed,
                                                     has_varying,
                                                     has_fixed_intercept,
                                                     has_varying_intercept,
                                                     J, K, K_fixed, K_varying,
                                                     L_fixed, L_varying,
                                                     S, ...) {
  if (noncentered) {
    xi_term <- ifelse_(shrinkage, " * xi[i - 1];", ";")
    declare_omega <- paste_rows(
      "// Spline coefficients",
      "matrix[K_varying_{y}, D] omega_{y}[S_{y} - 1];",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega <- paste_rows(
      "for (s in 1:(S_{y} - 1)) {{",
      "omega_{y}[s, , 1] = omega_raw_{y}[s, , 1];",
      "for (i in 2:D) {{",
      paste0(
        "omega_{y}[s, , i] = omega_{y}[s, , i - 1] + ",
        "omega_raw_{y}[s, , i] .* tau_{y}{xi_term}"
      ),
      "}}",
      "}}",
      .indent = idt(c(1, 2, 2, 3, 2, 1)),
      .parse = FALSE
    )
  }
  declare_delta <- paste_rows(
    "// Varying coefficients",
    "matrix[K_varying_{y}, S_{y} - 1] delta_{y}[T];",
    .indent = idt(1),
    .parse = FALSE
  )
  state_delta <- paste_rows(
    "for (s in 1:(S_{y} - 1)) {{",
    "for (t in 1:T) {{",
    "delta_{y}[t, , s] = omega_{y}[s] * Bs[, t];",
    "}}",
    "}}",
    .indent = idt(c(1, 2, 3, 2, 1)),
    .parse = FALSE
  )

  if (has_fixed || has_varying) {
    declare_omega_alpha_1 <- paste_rows(
      "// Fixed intercept",
      "vector[S_{y} - 1] alpha_{y}[T];",
      "real omega_alpha_1_{y}[S_{y} - 1];",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha_1 <- paste_rows(
      "// Define the first alpha using mean a_{y}",
      "for (s in 1:(S_{y} - 1)) {{",
      "vector[K_{y}] gamma__{y};",
      onlyif(
        has_fixed,
        "gamma__{y}[L_fixed_{y}] = beta_{y}[, s];"
      ),
      onlyif(
        has_varying,
        "gamma__{y}[L_varying_{y}] = delta_{y}[1, , s];"
      ),
      "omega_alpha_1_{y}[s] = a_{y}[s] - X_m[J_{y}] * gamma__{y};",
      "}}",
      .indent = idt(c(1, 1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
    declare_fixed_intercept <- paste_rows(
      "// Time-invariant intercept",
      "vector[S_{y} - 1] alpha_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_fixed_intercept <- paste_rows(
      "for (s in 1:(S_{y} - 1)) {{",
      "vector[K_{y}] gamma__{y};",
      onlyif(
        has_fixed,
        "gamma__{y}[L_fixed_{y}] = beta_{y}[, s];"
      ),
      onlyif(
        has_varying,
        "gamma__{y}[L_varying_{y}] = delta_{y}[1, , s];"
      ),
      "alpha_{y}[s] = a_{y}[s] - X_m[J_{y}] * gamma__{y};",
      "}}",
      .indent = idt(c(1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
  } else {
    declare_omega_alpha_1 <- paste_rows(
      "// Fixed intercept",
      "vector[S_{y} - 1] alpha_{y}[T];",
      "real omega_alpha_1_{y}[S_{y} - 1];",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha_1 <- paste_rows(
      "for (s in 1:(S_{y} - 1)) {{",
      "omega_alpha_1_{y}[s] = a_{y}[s];",
      "}}",
      .indent = idt(c(1, 2, 1)),
      .parse = FALSE
    )
    declare_fixed_intercept <- paste_rows(
      "// Time-invariant intercept",
      "vector[S_{y} - 1] alpha_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_fixed_intercept <- paste_rows(
      "for (s in 1:(S_{y} - 1)) {{",
      "alpha_{y}[s] = a_{y}[s];",
      "}}",
      .indent = idt(c(1, 2, 1)),
      .parse = FALSE
    )
  }
  if (noncentered) {
    xi_term <- ifelse_(shrinkage, " * xi[i - 1];", ";")
    declare_omega_alpha <- paste_rows(
      "// Spline coefficients",
      "row_vector[D] omega_alpha_{y}[S_{y} - 1];",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha <- paste_rows(
      "for (s in 1:(S_y - 1)) {{",
      "omega_alpha_{y}[s, 1] = omega_alpha_1_{y}[s];",
      "for (i in 2:D) {{",
      paste0(
        "omega_alpha_{y}[s, i] = omega_alpha_{y}[s, i - 1] + ",
        "omega_raw_alpha_{y}[s, i - 1] * tau_alpha_{y}{xi_term}"
      ),
      "}}",
      "}}",
      .indent = idt(c(1, 2, 2, 3, 2, 1)),
      .parse = FALSE
    )
  } else {
    declare_omega_alpha <- paste_rows(
      "// Spline coefficients",
      "row_vector[D] omega_alpha_{y}[S_{y} - 1];",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha <- paste_rows(
      "for (s in 1:(S_{y} - 1)) {{",
      "omega_alpha_{y}[s, 1] = omega_alpha_1_{y}[s];",
      "omega_alpha_{y}[s, 2:D] = omega_raw_alpha_{y}[s];",
      "}}",
      .indent = idt(c(1, 2, 2, 1)),
      .parse = FALSE
    )
  }
  declare_varying_intercept <- paste_rows(
    declare_omega_alpha_1,
    declare_omega_alpha,
    .indent = idt(c(0, 1)),
    .parse = FALSE
  )
  state_varying_intercept <- paste_rows(
    state_omega_alpha_1,
    state_omega_alpha,
    "for (t in 1:T) {{",
    "for (s in 1:(S_{y} - 1)) {{",
    "alpha_{y}[t, s] = omega_alpha_{y}[s] * Bs[, t];",
    "}}",
    "}}",
    .indent = idt(c(0, 0, 1, 2, 3, 2, 1)),
    .parse = FALSE
  )

  list(
    declarations = paste_rows(
      onlyif(has_varying && noncentered, declare_omega),
      onlyif(has_varying, declare_delta),
      onlyif(has_fixed_intercept, declare_fixed_intercept),
      onlyif(has_varying_intercept, declare_varying_intercept),
      .indent = idt(c(1, 1, 0, 0))
    ),
    statements = paste_rows(
      onlyif(has_varying && noncentered, state_omega),
      onlyif(has_varying, state_delta),
      onlyif(has_fixed_intercept, state_fixed_intercept),
      onlyif(has_varying_intercept, state_varying_intercept),
      .indent = idt(0)
    )
  )
}

transformed_parameters_lines_gaussian <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in%
    names(formals(transformed_parameters_lines_default))]
  do.call(what = transformed_parameters_lines_default, args = args)
}

transformed_parameters_lines_binomial <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in%
    names(formals(transformed_parameters_lines_default))]
  do.call(what = transformed_parameters_lines_default, args = args)
}

transformed_parameters_lines_bernoulli <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in%
    names(formals(transformed_parameters_lines_default))]
  do.call(what = transformed_parameters_lines_default, args = args)
}

transformed_parameters_lines_poisson <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in%
    names(formals(transformed_parameters_lines_default))]
  do.call(what = transformed_parameters_lines_default, args = args)
}

transformed_parameters_lines_negbin <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in%
    names(formals(transformed_parameters_lines_default))]
  do.call(what = transformed_parameters_lines_default, args = args)
}

transformed_parameters_lines_exponential <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in%
    names(formals(transformed_parameters_lines_default))]
  do.call(what = transformed_parameters_lines_default, args = args)
}

transformed_parameters_lines_gamma <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in%
    names(formals(transformed_parameters_lines_default))]
  do.call(what = transformed_parameters_lines_default, args = args)
}

transformed_parameters_lines_beta <- function(...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in%
    names(formals(transformed_parameters_lines_default))]
  do.call(what = transformed_parameters_lines_default, args = args)
}

# Model block -------------------------------------------------------------

model_lines_default <- function(y, idt, obs, noncentered, shrinkage,
                                has_varying, has_fixed, has_random,
                                has_fixed_intercept,
                                has_varying_intercept, has_random_intercept,
                                has_lfactor, noncentered_psi,
                                noncentered_lambda, nonzero_lambda,
                                sigma_nu_prior_distr = "",
                                sigma_nu_prior_npars = 1L,
                                alpha_prior_distr = "",
                                tau_alpha_prior_distr = "",
                                lambda_prior_distr = "",
                                beta_prior_distr = "",
                                beta_prior_npars = 1L,
                                delta_prior_distr = "",
                                delta_prior_npars = 1L,
                                tau_prior_distr = "",
                                tau_prior_npars = 1L,
                                sigma_lambda_prior_distr = "",
                                psi_prior_distr = "",
                                tau_psi_prior_distr = "",
                                K_fixed, K_varying, K_random, J_random, ...) {
  allow_intercept <- !(has_lfactor && nonzero_lambda)

  mtext_alpha <- "a_{y} ~ {alpha_prior_distr};"

  if (vectorizable_prior(sigma_nu_prior_distr)) {
    dpars_sigma_nu <- ifelse_(sigma_nu_prior_npars > 0L,
      paste0(
        "sigma_nu_prior_pars_", y, "[, ", seq_len(sigma_nu_prior_npars), "]",
        collapse = ", "
      ),
      ""
    )
    mtext_sigma_nu <- "sigma_nu_{y} ~ {sigma_nu_prior_distr}({dpars_sigma_nu});"
  } else {
    mtext_sigma_nu <-
      "sigma_nu_{y}[{{{cs(1:K_random)}}}] ~ {sigma_nu_prior_distr};"
  }

  if (has_lfactor) {
    m <- ifelse_(nonzero_lambda, "1", "0")
    mtext_lambda <- paste_rows(
      ifelse_(
        noncentered_lambda,
        "lambda_std_{y} ~ normal(0, inv(sqrt(1 - inv(N))));",
        "lambda_{y} ~ normal({m}, sigma_lambda_{y} * inv(sqrt(1 - inv(N))));"
      ),
      "sigma_lambda_{y} ~ {sigma_lambda_prior_distr};",
      onlyif(nonzero_lambda, "tau_psi_{y} ~ {tau_psi_prior_distr};"),
      ifelse_(
        noncentered_psi,
        "omega_raw_psi_1_{y} ~ {psi_prior_distr};",
        "omega_psi_{y}[1] ~ {psi_prior_distr};"
      ),
      .indent = idt(c(0, 1, 1, 1)),
      .parse = TRUE
    )
  }

  if (noncentered) {
    mtext_omega <- "omega_raw_alpha_{y} ~ std_normal();"
    if (vectorizable_prior(delta_prior_distr)) {
      dpars_varying <- ifelse_(delta_prior_npars > 0L,
        paste0(
          "delta_prior_pars_", y, "[, ", seq_len(delta_prior_npars), "]",
          collapse = ", "
        ),
        ""
      )
      mtext_varying <-
        "omega_raw_{y}[, 1] ~ {delta_prior_distr}({dpars_varying});"
    } else {
      mtext_varying <-
        "omega_raw_{y}[{{{cs(1:K_varying)}}}, 1] ~ {delta_prior_distr};"
    }
    mtext_varying <- paste_rows(
      mtext_varying,
      "to_vector(omega_raw_{y}[, 2:D]) ~ std_normal();",
      .indent = idt(c(0, 1)),
      .parse = FALSE
    )
  } else {
    xi_term1 <- ifelse_(shrinkage, " * xi[1]", "")
    xi_term <- ifelse_(shrinkage, " * xi[i - 1]", "")
    mtext_omega <- paste_rows(
      paste0(
        "omega_raw_alpha_{y}[1] ~ normal(omega_alpha_1_{y}, ",
        "tau_alpha_{y}{xi_term1});"
      ),
      "for (i in 2:(D - 1)) {{",
      paste0(
        "omega_raw_alpha_{y}[i] ~ normal(omega_raw_alpha_{y}[i - 1], ",
        "tau_alpha_{y}{xi_term});"
      ),
      "}}",
      .indent = idt(c(0, 1, 2, 1)),
      .parse = FALSE
    )
    if (vectorizable_prior(delta_prior_distr)) {
      dpars_varying <- ifelse_(delta_prior_npars > 0L,
        paste0(
          "delta_prior_pars_", y, "[, ", seq_len(delta_prior_npars), "]",
          collapse = ", "
        ),
        ""
      )
      mtext_varying <-
        "omega_{y}[, 1] ~ {delta_prior_distr}({dpars_varying});"
    } else {
      mtext_varying <-
        "omega_{y}[{{{cs(1:K_varying)}}}, 1] ~ {delta_prior_distr};"
    }
    mtext_varying <- paste_rows(
      mtext_varying,
      "for (i in 2:D) {{",
      ifelse_(
        shrinkage,
        paste0(
          "omega_{y}[, i] ~ normal(omega_{y}[, i - 1], ",
          "xi[i - 1] * tau_{y});"
        ),
        "omega_{y}[, i] ~ normal(omega_{y}[, i - 1], tau_{y});"
      ),
      "}}",
      .indent = idt(c(0, 1, 2, 1)),
      .parse = FALSE
    )
  }

  mtext_varying_intercept <- paste_rows(
    mtext_alpha,
    mtext_omega,
    "tau_alpha_{y} ~ {tau_alpha_prior_distr};",
    .indent = idt(c(0, 1, 1)),
    .parse = FALSE
  )

  mtext_fixed_intercept <- mtext_alpha

  if (vectorizable_prior(beta_prior_distr)) {
    dpars_fixed <-  ifelse_(beta_prior_npars > 0L,
      paste0(
        "beta_prior_pars_", y, "[, ", seq_len(beta_prior_npars), "]",
        collapse = ", "
      ),
      ""
    )
    mtext_fixed <- "beta_{y} ~ {beta_prior_distr}({dpars_fixed});"
  } else {
    mtext_fixed <- "beta_{y}[{{{cs(1:K_fixed)}}}] ~ {beta_prior_distr};"
  }

  if (vectorizable_prior(tau_prior_distr)) {
    dpars_tau <-  ifelse_(tau_prior_npars > 0L,
      paste0(
        "tau_prior_pars_", y, "[, ", seq_len(tau_prior_npars), "]",
        collapse = ", "
      ),
      ""
    )
    mtext_tau <- "tau_{y} ~ {tau_prior_distr}({dpars_tau});"
  } else {
    mtext_tau <- "tau_{y}[{{{cs(1:K_varying)}}}] ~ {tau_prior_distr};"
  }

  intercept_alpha <- ifelse(
    allow_intercept,
    ifelse_(
      has_fixed_intercept,
      glue::glue("alpha_{y}"),
      ifelse_(
        has_varying_intercept,
        glue::glue("alpha_{y}[t]"),
        ""
      )
    ),
    ""
  )
  intercept_nu <- ifelse_(
    has_random_intercept,
    ifelse(
      nzchar(obs),
      glue::glue("nu_{y}[{obs}, 1]"),
      glue::glue("nu_{y}[, 1]")
    ),
    ""
  )
  lfactor <- ifelse_(
    has_lfactor,
    ifelse(
      nzchar(obs),
      glue::glue("lambda_{y}[{obs}] * psi_{y}[t]"),
      glue::glue("lambda_{y} * psi_{y}[t]")
    ),
    ""
  )
  plus <- ifelse_(nzchar(intercept_alpha) && has_random_intercept, " + ", "")
  intercept <- ifelse_(
    nzchar(intercept_alpha) || has_random_intercept,
    glue::glue("{intercept_alpha}{plus}{intercept_nu}"),
    "0"
  )
  plus <- ifelse_(nzchar(intercept) && has_lfactor, " + ", "")
  intercept <- ifelse_(
    nzchar(intercept) || has_lfactor,
    glue::glue("{intercept}{plus}{lfactor}"),
    ""
  )
  random <- ifelse_(
    has_random,
    ifelse_(
      has_random_intercept,
      glue::glue("rows_dot_product(X[t][{obs}, L_random_{y}], nu_{y}[{obs}, 2:K_random_{y}])"),
      glue::glue("rows_dot_product(X[t][{obs}, L_random_{y}], nu_{y}[{obs}, ])")
    ),
    ""
  )
  list(text = paste_rows(
    onlyif(has_lfactor, mtext_lambda),
    onlyif(has_random_intercept || has_random, mtext_sigma_nu),
    onlyif(has_fixed_intercept && allow_intercept, mtext_fixed_intercept),
    onlyif(has_varying_intercept && allow_intercept, mtext_varying_intercept),
    onlyif(has_fixed, mtext_fixed),
    onlyif(has_varying, mtext_varying),
    onlyif(has_varying, mtext_tau),
    .indent = idt(c(1, 1, 1, 1, 1, 1, 1))
  ), intercept = intercept, random = random)
}

model_lines_categorical <- function(y, idt, obs, noncentered, shrinkage,
                                    has_varying, has_fixed,
                                    has_fixed_intercept, has_varying_intercept,
                                    has_random_intercept,
                                    has_lfactor,
                                    alpha_prior_distr = "",
                                    alpha_prior_npars = 1L,
                                    tau_alpha_prior_distr = "",
                                    lambda_prior_distr = "",
                                    beta_prior_distr = "",
                                    beta_prior_npars = 1L,
                                    delta_prior_distr = "",
                                    delta_prior_npars = 1L,
                                    tau_prior_distr = "",
                                    tau_prior_npars = 1L,
                                    J, K, K_fixed, K_varying, J_random,
                                    K_random, L_fixed, L_varying, S, backend,
                                    ...) {
  if (vectorizable_prior(alpha_prior_distr)) {
    dpars_alpha <- ifelse_(alpha_prior_npars > 0L,
      paste0(
        "alpha_prior_pars_", y, "[, ", seq_len(alpha_prior_npars), "]",
        collapse = ", "
      ),
      ""
    )
    mtext_alpha <- "to_vector(a_{y}) ~ {alpha_prior_distr}({dpars_alpha});"
  } else {
    s <- seq_len(S - 1L)
    mtext_alpha <- "a_{y}[{s}] ~ {alpha_prior_distr};"
  }
  mtext_fixed_intercept <- mtext_alpha

  xi_term1 <- ifelse_(shrinkage, " * xi[1]", "")
  xi_term <- ifelse_(shrinkage, " * xi[i]", "")
  mtext_omega <- ifelse_(
    noncentered,
    paste_rows(
      "for (s in 1:(S_{y} - 1)) {{",
      "omega_raw_alpha_{y}[s] ~ std_normal();",
      "}}",
      .indent = idt(c(1, 2, 1)),
      .parse = FALSE
    ),
    paste_rows(
      "for (s in 1:(S_{y} - 1)) {{",
      paste0(
        "omega_raw_alpha_{y}[s, 1] ~ normal(omega_alpha_1_{y}[s], ",
        "tau_alpha_{y}{xi_term1});"
      ),
      "for (i in 2:(D - 1)) {{",
      paste0(
        "omega_raw_alpha_{y}[s, i] ~ ",
        "normal(omega_raw_alpha_{y}[s, i - 1], ",
        "tau_alpha_{y}{xi_term});"
      ),
      "}}",
      "}}",
      .indent = idt(c(1, 2, 2, 3, 2, 1)),
      .parse = FALSE
    )
  )
  mtext_varying_intercept <- paste_rows(
    mtext_alpha,
    mtext_omega,
    "tau_alpha_{y} ~ {tau_alpha_prior_distr};",
    .indent = idt(c(0, 0, 1)),
    .parse = FALSE
  )

  if (vectorizable_prior(beta_prior_distr)) {
    dpars_fixed <- ifelse_(beta_prior_npars > 0L,
      paste0(
        "beta_prior_pars_", y, "[, ", seq_len(beta_prior_npars), "]",
        collapse = ", "
      ),
      ""
    )
    mtext_fixed <- "to_vector(beta_{y}) ~ {beta_prior_distr}({dpars_fixed});"
  } else {
    k <- rep(seq_len(K_fixed), S - 1L)
    s <- rep(seq_len(S - 1L), each = K)
    mtext_fixed <- "beta_{y}[{k},{s}] ~ {beta_prior_distr};"
  }

  if (noncentered) {
    if (vectorizable_prior(delta_prior_distr)) {
      dpars_varying <- ifelse_(delta_prior_npars > 0L,
        paste0(
          "delta_prior_pars_", y, "[, ", seq_len(delta_prior_npars), "]",
          collapse = ", "
        ),
        ""
      )
      mtext_varying <-
        "to_vector(delta_{y}[1]) ~ {delta_prior_distr}({dpars_varying});"
    } else {
      k <- rep(seq_len(K_fixed), S - 1L)
      s <- rep(seq_len(S - 1L), each = K_fixed)
      mtext_varying <- "omega_raw_{y}[{s},{k},1] ~ {delta_prior_distr};"
    }
    mtext_varying <- paste_rows(
      mtext_varying,
      "for (s in 1:(S_{y} - 1)) {{",
      "to_vector(omega_raw_{y}[s, ,2:D]) ~ std_normal();",
      "}}",
      .indent = idt(c(0, 1, 2, 1)),
      .parse = FALSE
    )
  } else {
    if (vectorizable_prior(delta_prior_distr)) {
      dpars_varying <- ifelse_(delta_prior_npars > 0L,
        paste0(
          "delta_prior_pars_", y, "[, ", seq_len(delta_prior_npars), "]",
          collapse = ", "
        ),
        ""
      )
      mtext_varying <-
        "to_vector(delta_{y}[1]) ~ {delta_prior_distr}({dpars_varying});"
    } else {
      k <- rep(seq_len(K_fixed), S - 1L)
      s <- rep(seq_len(S - 1L), each = K_varying)
      mtext_varying <- "omega_{y}[{s}, {k}, 1] ~ {delta_prior_distr};"
    }
    mtext_varying <- paste_rows(
      mtext_varying,
      "for (s in 1:(S_{y} - 1)) {{",
      "for (i in 2:D) {{",
      ifelse_(
        shrinkage,
        paste0(
          "omega_{y}[s, , i] ~ normal(omega_{y}[s, , i - 1], ",
          "xi[i - 1] * tau_{y});"
        ),
        "omega_{y}[s, , i] ~ normal(omega_{y}[s, , i - 1], tau_{y});"
      ),
      "}}",
      "}}",
      .indent = idt(c(0, 1, 2, 3, 2, 1)),
      .parse = FALSE
    )
  }
  if (vectorizable_prior(tau_prior_distr)) {
    dpars_tau <- ifelse_(tau_prior_npars > 0L,
      paste0(
        "tau_prior_pars_", y, "[, ", seq_len(tau_prior_npars), "]",
        collapse = ", "
      ),
      "")
    mtext_tau <- "tau_{y} ~ {tau_prior_distr}({dpars_tau});"
  } else {
    mtext_tau <- "tau_{y}[{{{cs(1:K_varying)}}}] ~ {tau_prior_distr};"
  }

  intercept <- ifelse_(
    has_fixed_intercept,
    glue::glue("append_row(0, alpha_{y})"),
    ifelse_(
      has_varying_intercept,
      glue::glue("append_row(0, alpha_{y}[t])"),
      glue::glue("zeros_S_{y}")
    )
  )
  # categorical_logit_glm does not support id-varying effects
  # and categorical distribution is very slow without it
  stopifnot_(
    identical(K_random, 0L),
    "Categorical family does not yet support random effects."
  )
  stopifnot_(
    !has_lfactor,
    "Categorical family does not yet support latent factors."
  )
  likelihood_term <- ifelse_(
    stan_supports_categorical_logit_glm(backend),
    ifelse_(
      has_fixed || has_varying,
      paste0(
        "y_{y}[t, {obs}] ~ categorical_logit_glm(X[t][{obs}, J_{y}], ",
        "{intercept}, append_col(zeros_K_{y}, gamma__{y}));"
      ),
      "y_{y}[t, {obs}] ~ categorical_logit({intercept});"
    ),
    ifelse_(
      has_fixed || has_varying,
      paste_rows(
        ifelse_(nzchar(obs), "for (i in {obs}) {{", "for (i in 1:N) {{"),
        paste0(
          "y_{y}[t, i] ~ categorical_logit({intercept} + ",
          "to_vector(X[t][i, J_{y}] * ",
          "append_col(zeros_K_{y}, gamma__{y})));"
        ),
        "}}",
        .indent = idt(c(0, 1, 0)),
        .parse = FALSE
      ),
      "y_{y}[t, {obs}] ~ categorical_logit({intercept});"
    )
  )
  paste_rows(
    onlyif(has_fixed_intercept, mtext_fixed_intercept),
    onlyif(has_varying_intercept, mtext_varying_intercept),
    onlyif(has_fixed, mtext_fixed),
    onlyif(has_varying, mtext_varying),
    onlyif(has_varying, mtext_tau),
    "{{",
    onlyif(
      has_fixed || has_varying,
      "matrix[K_{y}, S_{y} - 1] gamma__{y};"
    ),
    onlyif(has_fixed, "gamma__{y}[L_fixed_{y}] = beta_{y};"),
    "for (t in 1:T) {{",
    onlyif(has_varying, "gamma__{y}[L_varying_{y}] = delta_{y}[t];"),
    likelihood_term,
    "}}",
    "}}",
    .indent = idt(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 2, 1))
  )
}

model_lines_gaussian <- function(y, idt, obs, has_fixed, has_varying,
                                 has_random, sigma_prior_distr = "", J, K,
                                 L_fixed, L_varying, K_random, J_random, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(model_lines_default))]
  args$glm <- TRUE
  mtext_def <- do.call(model_lines_default, args = args)
  intercept <- mtext_def$intercept
  random <- mtext_def$random
  plus_re <- ifelse_(
    has_random && (nchar(intercept)),
    " + ",
    ""
  )
  likelihood_term <- ifelse_(
    has_fixed || has_varying,
    paste0(
      "y_{y}[{obs}, t] ~ normal_id_glm(X[t][{obs}, J_{y}], ",
      "{intercept}{plus_re}{random}, gamma__{y}, sigma_{y});"
    ),
    "y_{y}[{obs}, t] ~ normal({intercept}{plus_re}{random}, sigma_{y});"
  )

  mtext <- paste_rows(
    "sigma_{y} ~ {sigma_prior_distr};",
    "{{",
    onlyif(
      has_fixed || has_varying,
      "vector[K_{y}] gamma__{y};"
    ),
    onlyif(has_fixed, "gamma__{y}[L_fixed_{y}] = beta_{y};"),
    "for (t in 1:T) {{",
    onlyif(has_varying, "gamma__{y}[L_varying_{y}] = delta_{y}[t];"),
    likelihood_term,
    "}}",
    "}}",
    .indent = idt(c(1, 1, 2, 2, 2, 3, 3, 2, 1))
  )
  paste_rows(mtext_def$text, mtext, .parse = FALSE)
}

model_lines_binomial <- function(y, idt, obs, has_varying, has_fixed,
                                 has_random, J_fixed, J_varying, J_random,
                                 K_random, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(model_lines_default))]
  mtext_def <- do.call(model_lines_default, args = args)
  intercept <- mtext_def$intercept
  random <- mtext_def$random
  plus_re <- ifelse_(
    has_random && (nchar(intercept) || has_fixed || has_varying),
    " + ",
    ""
  )
  fixed_term <- ifelse_(
    has_fixed,
    glue::glue("X[t][{obs}, J_fixed_{y}] * beta_{y}"),
    ""
  )
  varying_term <- ifelse_(
    has_varying,
    glue::glue("X[t][{obs}, J_varying_{y}] * delta_{y}[t]"),
    ""
  )

  plus_c <- ifelse_(has_fixed && has_varying, " + ", "")
  plus_ic <- ifelse_(
    nzchar(intercept) && (has_fixed || has_varying),
    " + ",
    ""
  )

  likelihood_term <- paste0(
    "y_{y}[t, {obs}] ~ binomial_logit(trials_{y}[t, {obs}], ",
    "{intercept}{plus_ic}{fixed_term}{plus_c}{varying_term}{plus_re}{random});"
  )

  mtext <- paste_rows(
    "for (t in 1:T) {{",
    likelihood_term,
    "}}",
    .indent = idt(c(1, 2, 1))
  )
  paste_rows(mtext_def$text, mtext, .parse = FALSE)
}

model_lines_bernoulli <- function(y, idt, obs, has_varying, has_fixed,
                                  has_random, J, K, J_random, K_random,
                                  L_fixed, L_varying, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(model_lines_default))]
  mtext_def <- do.call(model_lines_default, args = args)
  intercept <- mtext_def$intercept
  random <- mtext_def$random
  plus_re <- ifelse_(
    has_random && (nchar(intercept)),
    " + ",
    ""
  )
  likelihood_term <- ifelse_(
    has_fixed || has_varying,
    paste0(
      "y_{y}[t, {obs}] ~ bernoulli_logit_glm(X[t][{obs}, J_{y}],",
      "{intercept}{plus_re}{random}, gamma__{y});"
    ),
    "y_{y}[t, {obs}] ~ bernoulli_logit({intercept}{plus_re}{random});"
  )
  mtext <- paste_rows(
    "{{",
    onlyif(
      has_fixed || has_varying,
      "vector[K_{y}] gamma__{y};"
    ),
    onlyif(has_fixed, "gamma__{y}[L_fixed_{y}] = beta_{y};"),
    "for (t in 1:T) {{",
    onlyif(has_varying, "gamma__{y}[L_varying_{y}] = delta_{y}[t];"),
    likelihood_term,
    "}}",
    "}}",
    .indent = idt(c(1, 2, 2, 2, 3, 3, 2, 1))
  )
  paste_rows(mtext_def$text, mtext, .parse = FALSE)
}

model_lines_poisson <- function(y, idt, obs, has_varying, has_fixed, has_random,
                                has_offset, J, K, J_random, K_random, L_fixed,
                                L_varying, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(model_lines_default))]
  mtext_def <- do.call(model_lines_default, args = args)
  intercept <- mtext_def$intercept
  intercept <- ifelse_(
    has_offset,
    ifelse_(
      intercept == "0",
      glue::glue("to_vector(offset_{y}[t, {obs}])"),
      glue::glue("{intercept} + to_vector(offset_{y}[t, {obs}])")
    ),
    intercept
  )
  random <- mtext_def$random
  plus_re <- ifelse_(
    has_random && (nchar(intercept)),
    " + ",
    ""
  )
  likelihood_term <- ifelse_(
    has_fixed || has_varying,
    paste0(
      "y_{y}[t, {obs}] ~ poisson_log_glm(X[t][{obs}, J_{y}], ",
      "{intercept}{plus_re}{random}, gamma__{y});"
    ),
    "y_{y}[t, {obs}] ~ poisson_log({intercept}{plus_re}{random});"
  )
  mtext <- paste_rows(
    "{{",
    onlyif(
      has_fixed || has_varying,
      "vector[K_{y}] gamma__{y};"
    ),
    onlyif(has_fixed, "gamma__{y}[L_fixed_{y}] = beta_{y};"),
    "for (t in 1:T) {{",
    onlyif(has_varying, "gamma__{y}[L_varying_{y}] = delta_{y}[t];"),
    likelihood_term,
    "}}",
    "}}",
    .indent = idt(c(1, 2, 2, 2, 3, 3, 2, 1))
  )
  paste_rows(mtext_def$text, mtext, .parse = FALSE)
}

model_lines_negbin <- function(y, idt, obs, has_varying, has_fixed, has_random,
                               has_offset, phi_prior_distr, J, K, J_random,
                               K_random, L_fixed, L_varying, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(model_lines_default))]
  mtext_def <- do.call(model_lines_default, args = args)
  intercept <- mtext_def$intercept
  intercept <- ifelse_(
    has_offset,
    ifelse_(
      intercept == "0",
      glue::glue("to_vector(offset_{y}[t, {obs}])"),
      glue::glue("{intercept} + to_vector(offset_{y}[t, {obs}])")
    ),
    intercept
  )
  random <- mtext_def$random
  plus_re <- ifelse_(
    has_random && (nchar(intercept)),
    " + ",
    ""
  )
  likelihood_term <- ifelse_(
    has_fixed || has_varying,
    paste0(
      "y_{y}[t, {obs}] ~ neg_binomial_2_log_glm(X[t][{obs}, J_{y}], ",
      "{intercept}{plus_re}{random}, gamma__{y}, phi_{y});"
    ),
    paste0(
      "y_{y}[t, {obs}] ~ neg_binomial_2_log(",
      "{intercept}{plus_re}{random}, phi_{y});"
    )
  )
  mtext <- paste_rows(
    "phi_{y} ~ {phi_prior_distr};",
    "{{",
    onlyif(
      has_fixed || has_varying,
      "vector[K_{y}] gamma__{y};"
    ),
    onlyif(has_fixed, "gamma__{y}[L_fixed_{y}] = beta_{y};"),
    "for (t in 1:T) {{",
    onlyif(has_varying, "gamma__{y}[L_varying_{y}] = delta_{y}[t];"),
    likelihood_term,
    "}}",
    "}}",
    .indent = idt(c(1, 1, 2, 2, 2, 3, 3, 2, 1))
  )
  paste_rows(mtext_def$text, mtext, .parse = FALSE)
}

model_lines_exponential <- function(y, idt, obs, has_varying, has_fixed,
                                    has_random, J_fixed, J_varying, J_random,
                                    K_random, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(model_lines_default))]
  mtext_def <- do.call(model_lines_default, args = args)
  intercept <- mtext_def$intercept
  random <- mtext_def$random
  plus_re <- ifelse_(
    has_random && (nchar(intercept) || has_fixed || has_varying),
    " + ",
    ""
  )
  fixed_term <- ifelse_(
    has_fixed,
    glue::glue("X[t][{obs}, J_fixed_{y}] * beta_{y}"),
    ""
  )
  varying_term <- ifelse_(
    has_varying,
    glue::glue("X[t][{obs}, J_varying_{y}] * delta_{y}[t]"),
    ""
  )
  plus_c <- ifelse_(has_fixed && has_varying, " + ", "")
  plus_ic <- ifelse_(
    nzchar(intercept) && (has_fixed || has_varying),
    " + ",
    ""
  )
  likelihood_term <- paste0(
    "y_{y}[{obs}, t] ~ exponential(exp(-({intercept}{plus_ic}",
    "{fixed_term}{plus_c}{varying_term}{plus_re}{random})));"
  )
  mtext <- paste_rows(
    "for (t in 1:T) {{",
    likelihood_term,
    "}}",
    .indent = idt(c(1, 2, 1))
  )
  paste_rows(mtext_def$text, mtext, .parse = FALSE)
}

model_lines_gamma <- function(y, idt, obs, has_varying, has_fixed, has_random,
                              phi_prior_distr, J_fixed, J_varying, J_random,
                              K_random, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(model_lines_default))]
  mtext_def <- do.call(model_lines_default, args = args)
  intercept <- mtext_def$intercept
  random <- mtext_def$random
  plus_re <- ifelse_(
    has_random && (nchar(intercept) || has_fixed || has_varying),
    " + ",
    ""
  )
  phi_term <- "phi_{y} ~ {phi_prior_distr};"
  fixed_term <- ifelse_(
    has_fixed,
    glue::glue("X[t][{obs}, J_fixed_{y}] * beta_{y}"),
    ""
  )
  varying_term <- ifelse_(
    has_varying,
    glue::glue("X[t][{obs}, J_varying_{y}] * delta_{y}[t]"),
    ""
  )

  plus_c <- ifelse_(has_fixed && has_varying, " + ", "")
  plus_ic <- ifelse_(
    nzchar(intercept) && (has_fixed || has_varying),
    " + ",
    ""
  )
  likelihood_term <- paste0(
    "y_{y}[{obs}, t] ~ gamma(phi_{y}, phi_{y} * ",
    "exp(-({intercept}{plus_ic}{fixed_term}{plus_c}{varying_term}{plus_re}{random})));"
  )
  mtext <- paste_rows(
    phi_term,
    "for (t in 1:T) {{",
    likelihood_term,
    "}}",
    .indent = idt(c(1, 1, 2, 1))
  )
  paste_rows(mtext_def$text, mtext, .parse = FALSE)
}

model_lines_beta <- function(y, idt, obs, has_varying, has_fixed, has_random,
                             phi_prior_distr, J_fixed, J_varying, J_random,
                             K_random, ...) {
  args <- as.list(match.call()[-1L])
  args <- args[names(args) %in% names(formals(model_lines_default))]
  mtext_def <- do.call(model_lines_default, args = args)
  intercept <- mtext_def$intercept
  random <- mtext_def$random
  plus_re <- ifelse_(
    has_random && (nchar(intercept) || has_fixed || has_varying),
    " + ",
    ""
  )
  fixed_term <- ifelse_(
    has_fixed,
    glue::glue("X[t][{obs}, J_fixed_{y}] * beta_{y}"),
    ""
  )
  varying_term <- ifelse_(
    has_varying,
    glue::glue("X[t][{obs}, J_varying_{y}] * delta_{y}[t]"),
    ""
  )
  plus_c <- ifelse_(has_fixed && has_varying, " + ", "")
  plus_ic <- ifelse_(
    nzchar(intercept) && (has_fixed || has_varying),
    " + ",
    ""
  )
  likelihood_term <- paste0(
    "y_{y}[{obs}, t] ~ beta_proportion(",
    "inv_logit(({intercept}{plus_ic}{fixed_term}{plus_c}{varying_term}{plus_re}{random})), ",
    "phi_{y});"
  )
  mtext <- paste_rows(
    "phi_{y} ~ {phi_prior_distr};",
    "for (t in 1:T) {{",
    likelihood_term,
    "}}",
    .indent = idt(c(1, 1, 2, 1))
  )
  paste_rows(mtext_def$text, mtext, .parse = FALSE)
}

# Generated quantities block ----------------------------------------------

# generated_quantities_lines_default <- function() {
#   ""
# }
#
# generated_quantities_lines_categorical <- function() {
#   generated_quantities_lines_default()
# }
#
# generated_quantities_lines_gaussian <- function() {
#   generated_quantities_lines_default()
# }
#
# generated_quantities_lines_binomial <- function() {
#   generated_quantities_lines_default()
# }
#
# generated_quantities_lines_bernoulli <- function() {
#   generated_quantities_lines_default()
# }
#
# generated_quantities_lines_poisson <- function() {
#   generated_quantities_lines_default()
# }
#
# generated_quantities_lines_negbin <- function() {
#   generated_quantities_lines_default()
# }
#
# generated_quantities_lines_exponential <- function() {
#   generated_quantities_lines_default()
# }
#
# generated_quantities_lines_gamma <- function() {
#   generated_quantities_lines_default()
# }
#
# generated_quantities_lines_beta <- function() {
#   generated_quantities_lines_default()
# }

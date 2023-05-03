#' Wrapper to Parse Stan Model Blocks Based on the Family of the Response
#'
#' @param prefix \[`character(1)`]\cr Stan model block name, e.g., "model".
#' @param family \[`dynamitefamily`, `character(1)`]\cr A family object, a
#'   supported family name, or `"default"`.
#' @param args Channel specific component of `channel_vars`.
#' @param idt An indenter function.
#' @param backend The Stan backend.
#'   (see [create_blocks()])
#' @noRd
lines_wrap <- function(prefix, family, idt, backend, args) {
  args$idt <- idt
  args$backend <- backend
  suffix <- ifelse_(
    is.dynamitefamily(family),
    family$name,
    family
  )
  do.call(what = paste0(prefix, "_lines_", suffix), args = args)
}

# Block parameters --------------------------------------------------------

#' Parameters for Stan Code Generation
#'
#' Each of the following functions of the form `(stan_block)_lines_(family)`
#' uses a subset from a  common set of arguments that are documented here.
#'
#' @param backend \[`character(1)`]\cr `"rstan"` or `"cmdstanr"`.
#' @param y \[`character(1)`]\cr The name of the channel or vector of
#'   channel names for multivariate channels.
#' @param t \[`character(1)`]\cr  The time index. Either `"t"` if there
#'   are observations at every time point, otherwise a channel-specific
#'   vector of indices with at least one observation.
#' @param obs \[`character(1)`]\cr A vector of indices indicating non-missing
#'   values of the channel.
#' @param t_obs \[`character(1)`]\cr Either `"T"` if the are observations
#'   at every time point, otherwise a channel specific integer that represents
#'   the number of time points with at least one observation for the channel.
#' @param default \[`character(1)`]\cr Default channel specifications for the
#'   block.
#' @param intercept \[`character(1)`]\cr Model block intercept definitions
#' @param priors \[`character(1)`]\cr Model block prior definitions
#' @param y_cg \[`character(1)`]\cr Name of the channel group for multivariate
#'   channels.
#' @param idt \[`function`]\cr An indentation function, see [indenter_()].
#' @param has_missing \[`logical(1)`]\cr Does the channel contain missing
#'   observations?
#' @param has_fully_missing \[`logical(1)`]\cr Does the channel have time
#'   points with zero observations?
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
#' @param vectorized_beta \[`logical(1)`]\cr If `TRUE`, use a vectorized prior
#'   for time-invariant predictors `beta`.
#' @param vectorized_delta \[`logical(1)`]\cr If `TRUE`, use a vectorized prior
#'   for time-varying predictors `delta`.
#' @param vectorized_tau \[`logical(1)`]\cr If `TRUE`, use a vectorized prior
#'   for `tau` parameters.
#' @param vectorized_sigma_nu \[`logical(1)`]\cr If `TRUE`, use a vectorized
#'   prior for `sigma_nu` parameters.
#' @param sigma_prior_distr \[`character(1)`]\cr `sigma` parameter prior
#'   specification.
#' @param sigma_nu_prior_distr \[`character(1)`]\cr `sigma_nu` parameter prior
#'   specification.
#' @param alpha_prior_distr \[`character(1)`]\cr `alpha` parameter prior
#'   specification.
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

# # Functions block -----------------------------------------------------------

intercept_lines <- function(y, obs, family, has_varying, has_fixed, has_random,
                            has_fixed_intercept,
                            has_varying_intercept, has_random_intercept,
                            has_lfactor, has_offset, backend, ydim = y, ...) {

  intercept_alpha <- ifelse_(
    has_fixed_intercept,
    glue::glue("alpha_{y}"),
    ifelse_(
      has_varying_intercept,
      glue::glue("alpha_{y}[t]"),
      "0"
    )
  )
  offset <- ifelse_(
    has_offset,
    glue::glue(" + offset_{y}[{obs}, t]"),
    ""
  )
  intercept_nu <- ifelse_(
    has_random_intercept,
    ifelse(
      nzchar(obs),
      glue::glue(" + nu_{y}[{obs}, 1]"),
      glue::glue(" + nu_{y}[, 1]")
    ),
    ""
  )
  random <- ifelse_(
    has_random,
    ifelse_(
      has_random_intercept,
      paste0(
        glue::glue(" + rows_dot_product(X[t][{obs}, L_random_{ydim}], "),
        glue::glue("nu_{y}[{obs}, 2:K_random_{ydim}])")
      ),
      paste0(
        glue::glue(" + rows_dot_product(X[t][{obs}, L_random_{ydim}], "),
        glue::glue("nu_{y}[{obs}, ])")
      )
    ),
    ""
  )
  lfactor <- ifelse_(
    has_lfactor,
    ifelse(
      nzchar(obs),
      glue::glue(" + lambda_{y}[{obs}] * psi_{y}[t]"),
      glue::glue(" + lambda_{y} * psi_{y}[t]")
    ),
    ""
  )
  fixed <- ifelse_(
    has_fixed,
    glue::glue(" + X[t][{obs}, J_fixed_{ydim}] * beta_{y}"),
    ""
  )
  varying <- ifelse_(
    has_varying,
    glue::glue(" + X[t][{obs}, J_varying_{ydim}] * delta_{y}[t]"),
    ""
  )
  common_intercept <- !has_random && !has_random_intercept && !has_lfactor
  glm <- stan_supports_glm_likelihood(family, backend, common_intercept) &&
    (has_fixed || has_varying)
  intercept <- ifelse_(
    glm,
    glue::glue("{intercept_alpha}{offset}{intercept_nu}{random}{lfactor}"),
    glue::glue("{intercept_alpha}{offset}{intercept_nu}{random}{lfactor}{fixed}{varying}")
  )
  intercept <- sub(
    "^0 \\+",
    "",
    intercept
  )
  attr(intercept, "glm") <- glm
  intercept
}
functions_lines_default <- function(y, idt, obs, family, has_missing,
                                    has_varying, has_fixed, has_random,
                                    has_fixed_intercept, has_varying_intercept,
                                    has_random_intercept, has_lfactor,
                                    has_offset, backend, threading, ydim = y,
                                    K, ...) {

  y_type <- ifelse_(
    family$name %in% c("gaussian", "exponential", "gamma", "beta", "student"),
    "matrix",
    ifelse_(
      stan_supports_array_keyword(backend),
      "array[,] int",
      "int[,]"
    )
  )

  has_X <- has_fixed || has_varying || has_random

  extra_pars <- onlyif(
    family$name %in% c("gaussian", "student"),
    "real sigma_{y}"
  )
  extra_pars <- ifelse_(
    family$name %in% c("negbin", "gamma", "beta", "student"),
    "real phi_{y}",
    extra_pars
  )
  intercept <- intercept_lines(
    y, obs, family, has_varying, has_fixed, has_random, has_fixed_intercept,
    has_varying_intercept, has_random_intercept, has_lfactor, has_offset,
    backend, ydim
  )
  glm <- attr(intercept, "glm")
  scalar_intercept <- !has_offset && !has_random && !has_random_intercept &&
    !has_lfactor && (glm || !has_X)
  n_obs <- ifelse_(
    nchar(obs),
    paste0("n_obs_", y, "[t]"),
    "N"
  )
  intercept_line <- ifelse_(
    scalar_intercept,
    "real intercept_{y} = {intercept};",
    "vector[{n_obs}] intercept_{y} = {intercept};"
  )
  fun_name <- paste0(family$name, "_loglik_", y, "_lpmf")

  LJ <- ifelse_(
    glm,
    "L",
    "J"
  )
  fun_args <- glue::glue(cs(c(
    stan_array_arg(backend, "int", "t_obs_{y}", 0L, TRUE),
    onlyif(threading, c("int start", "int end")),
    ifelse_(
      has_missing,
      c(
        stan_array_arg(backend, "int", "obs_{y}", 1L, TRUE),
        stan_array_arg(backend, "int", "n_obs_{y}", 0L, TRUE)
      ),
      "data int N"
    ),
    paste("data", y_type, "y_{y}"),
    onlyif(has_fixed_intercept, "real alpha_{y}"),
    onlyif(
      has_varying_intercept,
      stan_array_arg(backend, "real", "alpha_{y}", 0L)
    ),
    onlyif(
      has_random,
      c(
        stan_array_arg(backend, "int", "L_random_{y}", 0L, TRUE),
        "int K_random_{y}"
      )
    ),
    onlyif(has_random || has_random_intercept, "matrix nu_{y}"),
    onlyif(has_lfactor, c("vector lambda_{y}", "vector psi_{y}")),
    onlyif(
      has_fixed,
      c(
        stan_array_arg(backend, "int", "{LJ}_fixed_{y}", 0L, TRUE),
        "vector beta_{y}"
      )
    ),
    onlyif(
      has_varying,
      c(
        stan_array_arg(backend, "int", "{LJ}_varying_{y}", 0L, TRUE),
        stan_array_arg(backend, "vector", "delta_{y}")
      )
    ),
    onlyif(
      has_fixed || has_varying,
      c(
        stan_array_arg(backend, "int", "J_{y}", 0L, TRUE),
        "data int K_{y}"
      )
    ),
    onlyif(has_X, "data array[] matrix X"),
    onlyif(has_offset, "data matrix offset_{y}"),
    onlyif(
      is_binomial(family),
      stan_array_arg(backend, "int", "trials_{y}", 1L, TRUE)
    ),
    extra_pars
  ))
  )

  fun_body <- paste_rows(
    "real ll = 0.0;",
    onlyif(glm, "vector[K_{y}] gamma__{y};"),
    onlyif(glm && has_fixed, "gamma__{y}[L_fixed_{y}] = beta_{y};"),
    "for (t in t_obs_{y}) {{",
    intercept_line,
    onlyif(glm && has_varying, "gamma__{y}[L_varying_{y}] = delta_{y}[t];"),
    "__likelihood__",
    "}}",
    "return ll;",
    .indent = idt(c(2, 2, 2, 2, 3, 3, 3, 2, 2))
  )
  list(
    fun_name = fun_name,
    fun_args = fun_args,
    fun_body = fun_body,
    use_glm = attr(intercept, "glm"),
    threading = threading)
}

functions_lines_gaussian <- function(y, obs, idt, default, ...) {

  likelihood <- ifelse_(
    default$use_glm,
    glue::glue(
      "ll += normal_id_glm_lupdf(y_{y}[{obs}, t] | X[t][{obs}, J_{y}], ",
      "intercept_{y}, gamma__{y}, sigma_{y});"
    ),
    glue::glue(
      "ll += normal_lupdf(y_{y}[{obs}, t] | intercept_{y}, sigma_{y});"
    )
  )
  fun_body <- default$fun_body
  fun_body <- gsub("__likelihood__", likelihood, fun_body)
  paste_rows(
    "real {default$fun_name}({default$fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}

functions_lines_binomial <- function(y, obs, idt, default, ...) {
  likelihood <- glue::glue(
    "ll += binomial_logit_lupmf(y_{y}[t, {obs}] | ",
    "trials_{y}[t, {obs}], intercept_{y});"
  )
  fun_body <- default$fun_body
  fun_body <- gsub("__likelihood__", likelihood, fun_body)
  paste_rows(
    "real {default$fun_name}({default$fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}

functions_lines_bernoulli <- function(y, obs, idt, default, ...) {
  likelihood <- ifelse_(
    default$use_glm,
    glue::glue(
      "ll += bernoulli_logit_glm_lupmf(y_{y}[t, {obs}] | X[t][{obs}, J_{y}], ",
      "intercept_{y}, gamma__{y});"
    ),
    glue::glue("ll += bernoulli_logit_lupmf(y_{y}[t, {obs}] | intercept_{y});")
  )
  fun_body <- default$fun_body
  fun_body <- gsub("__likelihood__", likelihood, fun_body)
  paste_rows(
    "real {default$fun_name}({default$fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}

functions_lines_poisson <- function(y, obs, idt, default, ...) {
  likelihood <- ifelse_(
    default$use_glm,
    glue::glue(
      "ll += poisson_log_glm_lupmf(y_{y}[t, {obs}] | X[t][{obs}, J_{y}], ",
      "intercept_{y}, gamma__{y});"
    ),
    glue::glue("ll += poisson_log_lupmf(y_{y}[t, {obs}] | intercept_{y});")
  )
  fun_body <- default$fun_body
  fun_body <- gsub("__likelihood__", likelihood, fun_body)
  paste_rows(
    "real {default$fun_name}({default$fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}

functions_lines_negbin <- function(y, obs, idt, default, ...) {
  likelihood <- ifelse_(
    default$use_glm,
    glue::glue(
      "ll += neg_binomial_2_log_glm_lupmf(y_{y}[t, {obs}] | ",
      "X[t][{obs}, J_{y}], intercept_{y}, gamma__{y}, phi_{y});"
    ),
    glue::glue(
      "ll += neg_binomial_2_log_lupmf(y_{y}[t, {obs}] | ",
      "intercept_{y}, phi_{y});"
    )
  )
  fun_body <- default$fun_body
  fun_body <- gsub("__likelihood__", likelihood, fun_body)
  paste_rows(
    "real {default$fun_name}({default$fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}

functions_lines_exponential <- function(y, obs, idt, default, ...) {
  likelihood <- glue::glue(
    "ll += exponential_lupdf(y_{y}[{obs}, t] | inv(exp(intercept_{y})));"
  )
  fun_body <- default$fun_body
  fun_body <- gsub("__likelihood__", likelihood, fun_body)
  paste_rows(
    "real {default$fun_name}({default$fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}

functions_lines_gamma <- function(y, obs, idt, default, ...) {
  likelihood <- glue::glue(
    "ll += gamma_lupdf(y_{y}[{obs}, t] | phi_{y}, phi_{y} * ",
    "inv(exp(intercept_{y})));"
  )
  fun_body <- default$fun_body
  fun_body <- gsub("__likelihood__", likelihood, fun_body)
  paste_rows(
    "real {default$fun_name}({default$fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}

functions_lines_beta <- function(y, obs, idt, default, ...) {
  likelihood <- glue::glue(
    "ll += beta_proportion_lupdf(y_{y}[{obs}, t] | ",
    "inv_logit(intercept_{y}), phi_{y});"
  )
  fun_body <- default$fun_body
  fun_body <- gsub("__likelihood__", likelihood, fun_body)
  paste_rows(
    "real {default$fun_name}({default$fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}
functions_lines_student <- function(y, obs, idt, default, ...) {
  likelihood <- glue::glue(
    "ll += student_t_lupdf(y_{y}[{obs}, t] | phi_{y}, intercept_{y}, ",
    "sigma_{y});"
  )
  fun_body <- default$fun_body
  fun_body <- gsub("__likelihood__", likelihood, fun_body)
  paste_rows(
    "real {default$fun_name}({default$fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}
functions_lines_categorical <- function(y, idt, obs, family, has_missing,
                                        has_varying, has_fixed, has_random,
                                        has_fixed_intercept,
                                        has_varying_intercept,
                                        has_random_intercept, has_lfactor,
                                        has_offset, backend, threading,
                                        ydim = y, K, categories,
                                        multinomial = FALSE, ...) {
  distr <- ifelse_(multinomial, "multinomial", "categorical")
  S <- length(categories)
  cats <- categories[seq.int(2L, S)]
  fun_args <- character(length(y))
  intercept <- vector("list", length(y))
  for (i in seq_along(cats)) {
    yi <- ifelse_(
      multinomial,
      cats[i],
      paste0(y, "_", cats[i])
    )
    intercept[[i]] <- intercept_lines(
      yi, obs, family, has_varying, has_fixed, has_random, has_fixed_intercept,
      has_varying_intercept, has_random_intercept, has_lfactor, has_offset,
      backend, ydim = y
    )
    fun_args[i] <- glue::glue(cs(c(
      onlyif(has_fixed_intercept, "real alpha_{yi}"),
      onlyif(
        has_varying_intercept,
        stan_array_arg(backend, "real", "alpha_{yi}", 0L)
      ),
      onlyif(has_random || has_random_intercept, "matrix nu_{yi}"),
      onlyif(has_lfactor, c("vector lambda_{yi}", "vector psi_{yi}")),
      onlyif(has_fixed, "vector beta_{yi}"),
      onlyif(has_varying, stan_array_arg(backend, "vector", "delta_{yi}"))
    )))
  }

  has_X <- has_fixed || has_varying || has_random
  glm <- attr(intercept[[1]], "glm")
  LJ <- ifelse_(glm, "L", "J")
  fun_args <- glue::glue(cs(c(
    stan_array_arg(backend, "int", "t_obs_{y}", 0L, TRUE),
    onlyif(threading, c("int start", "int end")),
    ifelse_(
      has_missing,
      c(stan_array_arg(backend, "int", "obs_{y}", 1L, TRUE),
        stan_array_arg(backend, "int", "n_obs_{y}", 0L, TRUE)
      ),
      "data int N"
    ),
    stan_array_arg(backend, "int", "y_{y}", 1 + multinomial, TRUE),
    "data int S_{y}",
    fun_args,
    onlyif(
      has_random,
      c(
        stan_array_arg(backend, "int", "L_random_{y}", 0L, TRUE),
        "int K_random_{y}"
      )
    ),
    onlyif(
      has_fixed,
      stan_array_arg(backend, "int", "{LJ}_fixed_{y}", 0L, TRUE)
    ),
    onlyif(
      has_varying,
      stan_array_arg(backend, "int", "{LJ}_varying_{y}", 0L, TRUE)
    ),
    onlyif(
      has_fixed || has_varying,
      c(
        stan_array_arg(backend, "int", "J_{y}", 0L, TRUE),
        "data int K_{y}"
      )
    ),
    onlyif(has_X, "data array[] matrix X")
  )))

  n_obs <- ifelse_(
    has_missing,
    glue::glue("n_obs_{y}[t]"),
    "N"
  )
  if (glm && !multinomial) {
    # combine intercepts and gammas
    icpt_y <- c("0", paste0("intercept_", cats))
    icpt <- paste_rows(
      "real intercept_{cats} = {unlist(intercept)};",
      "vector[S_{y}] intercept_{y} = [{cs(icpt_y)}]';",
      .indent = idt(1)
    )
    gamma <- onlyif(
      has_fixed || has_varying,
      paste_rows(
        "matrix[K_{y}, S_{y}] gamma__{y};",
        "gamma__{y}[, 1] = zeros_K_{y};",
        .indent = idt(1)
      )
    )
    beta <- onlyif(
      has_fixed,
      paste0(
        "gamma__{y}[L_fixed_{y}, {seq.int(2L, S)}] = ",
        "beta_{y}_{cats};"
      )
    )
    delta <- onlyif(
      has_varying,
      paste0(
        "gamma__{y}[L_varying_{y}, {seq.int(2L, S)}] = ",
        "delta_{y}_{cats}[t];"
      )
    )
    likelihood_term <- paste0(
      "ll += categorical_logit_glm_lupmf(y_{y}[t, {obs}] | X[t][{obs},",
      " J_{y}], intercept_{y}, gamma__{y});"
    )
    likelihood <- paste_rows(
      gamma,
      beta,
      "for (t in t_obs_{y}) {{",
      icpt,
      delta,
      likelihood_term,
      "}}",
      .indent = idt(c(0, 0, 0, 1, 2, 2, 1))
    )
  } else {
    category_dim <- ifelse_(multinomial, ", ", "")
    icpt_y <- c("0", paste0("intercept_", cats, "[i]"))
    icpt <- paste_rows(
      ifelse_(
        nzchar(obs),
        "vector[n_obs_{y}[t]] intercept_{cats} = {unlist(intercept)};",
        "vector[N] intercept_{cats} = {unlist(intercept)};"
      ),
      .indent = idt(2)
    )
    likelihood_term <- ifelse_(
      nzchar(obs),
      paste0("ll += {distr}_logit_lupmf(y_{y}[t, obs_{y}[i, t]{category_dim}]",
             "| intercept_{y});"),
      "ll += {distr}_logit_lupmf(y_{y}[t, i{category_dim}] | intercept_{y});"
    )
    likelihood <- paste_rows(
      "vector[S_{y}] intercept_{y};",
      "for (t in t_obs_{y}) {{",
      icpt,
      ifelse_(
        nzchar(obs),
        "for (i in 1:n_obs_{y}[t]) {{",
        "for (i in 1:N) {{"
      ),
      "intercept_{y} = [{cs(icpt_y)}]';",
      likelihood_term,
      "}}",
      "}}",
      .indent = idt(c(1, 1, 0, 2, 3, 3, 2, 1))
    )
  }
  fun_body <- paste_rows(
    onlyif(
      K > 0L,
      glue::glue("vector[K_{y}] zeros_K_{y} = rep_vector(0, K_{y});")
    ),
    "real ll = 0.0;",
    likelihood,
    "return ll;",
    .parse = FALSE,
    .indent = idt(c(1, 1, 0, 1))
  )
  paste_rows(
    "real {distr}_loglik_{y}_lpmf({fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}
functions_lines_multinomial <- function(idt, cvars, cgvars, backend,
                                        threading, ...) {
  cgvars$categories <- cgvars$y
  cgvars$y <- cgvars$y_cg
  cgvars$multinomial <- TRUE
  cgvars$threading <- threading
  cgvars$backend <- backend
  do.call(functions_lines_categorical, args = c(cgvars, idt = idt))
}
functions_lines_mvgaussian <- function(idt, cvars, cgvars, backend,
                                       threading, ...) {

  y <- cgvars$y
  y_cg <- cgvars$y_cg
  obs <- cgvars$obs
  mu <- fun_args <- character(length(y))
  has_X <- FALSE
  for (i in seq_along(y)) {
    yi <- y[i]
    args <- cvars[[i]]
    args$obs <- obs
    args$backend <- backend
    mu[i] <- do.call(intercept_lines, args = args)
    fun_args[i] <- glue::glue(cs(c(
      onlyif(cvars[[i]]$has_fixed_intercept, "real alpha_{yi}"),
      onlyif(
        cvars[[i]]$has_varying_intercept,
        stan_array_arg(backend, "real", "alpha_{yi}")
      ),
      onlyif(
        cvars[[i]]$has_random,
        c(
          stan_array_arg(backend, "int", "L_random_{yi}", 0L, TRUE),
          "int K_random_{yi}"
        )
      ),
      onlyif(
        cvars[[i]]$has_random || cvars[[i]]$has_random_intercept,
        "matrix nu_{yi}"
      ),
      onlyif(
        cvars[[i]]$has_lfactor,
        c("vector lambda_{yi}", "vector psi_{yi}")
      ),
      onlyif(
        cvars[[i]]$has_fixed,
        c(
          stan_array_arg(backend, "int", "J_fixed_{yi}", 0L, TRUE),
          "vector beta_{yi}"
        )
      ),
      onlyif(
        cvars[[i]]$has_varying,
        c(
          stan_array_arg(backend, "int", "J_varying_{yi}", 0L, TRUE),
          stan_array_arg(backend, "vector", "delta_{yi}")
        )
      ),
      "real sigma_{yi}"
    )))
    has_X <- has_X ||
      cvars[[i]]$has_fixed || cvars[[i]]$has_varying || cvars[[i]]$has_random
  }

  fun_args <- glue::glue(cs(c(
    stan_array_arg(backend, "int", "t_obs_{y_cg}", 0L, TRUE),
    onlyif(threading, c("int start", "int end")),
    ifelse_(
      cgvars$has_missing,
      c(stan_array_arg(backend, "int", "obs_{y_cg}", 1L, TRUE),
        stan_array_arg(backend, "int", "n_obs_{y_cg}", 0L, TRUE)
      ),
      "data int N"
    ),
    paste0("data vector[,] y_{y_cg}"),
    "int O_{y_cg}",
    fun_args,
    onlyif(has_X, "data array[] matrix X"),
    "matrix L_{y_cg}"
  )))
  n_obs <- ifelse_(
    cgvars$has_missing,
    glue::glue("n_obs_{y_cg}[t]"),
    "N"
  )
  sd_y <- paste0("sigma_", y)
  mu_y <- paste0("mu_", y, "[i]")
  likelihood <- paste_rows(
    "vector[O_{y_cg}] sigma_{y_cg} = [{cs(sd_y)}]';",
    paste0(
      "matrix[O_{y_cg}, O_{y_cg}] Lsigma = ",
      "diag_pre_multiply(sigma_{y_cg}, L_{y_cg});"
    ),
    "for (t in t_obs_{y_cg}) {{",
    "vector[O_{y_cg}] mu[{n_obs}];",
    "vector[{n_obs}] mu_{y} = {mu};",
    "for (i in 1:{n_obs}) {{",
    "mu[i] = [{cs(mu_y)}]';",
    "}}",
    "ll += multi_normal_cholesky_lupdf(y_{y_cg}[t, {obs}] | mu, Lsigma);",
    "}}",
    .indent = idt(c(1, 1, 1, 2, 2, 2, 3, 2, 2, 1))
  )
  fun_body <- paste_rows(
    "real ll = 0.0;",
    likelihood,
    "return ll;",
    .parse = FALSE,
    .indent = idt(c(1, 0, 1))
  )
  paste_rows(
    "real mvgaussian_loglik_{y_cg}_lpmf({fun_args}) {{",
    "{fun_body}",
    "}}",
    .indent = idt(c(1, 0, 1))
  )
}

# Data block --------------------------------------------------------------

missing_data_lines <- function(y, idt,
                               has_missing, has_fully_missing, backend) {
  mis <- character(0L)
  full_mis <- character(0L)
  if (has_missing) {
    mis <- paste_rows(
      "// Missing data indicators",
      stan_array(backend, "int", "obs_{y}", "N, T", "lower=0"),
      stan_array(backend, "int", "n_obs_{y}", "T", "lower=0"),
      .indent = idt(1),
      .parse = TRUE
    )
  }
  if (has_fully_missing) {
    full_mis <- paste_rows(
      "int<lower=0> T_obs_{y};",
      stan_array(backend, "int", "t_obs_{y}", "T_obs_{y}", "lower=0"),
      .indent = idt(1),
      .parse = TRUE
    )
  }
  paste_rows(
    mis,
    full_mis,
    .indent = idt(0),
    .parse = FALSE
  )
}

prior_data_lines <- function(y, idt, prior_distr,
                             K_fixed, K_varying, K_random, category = "",
                             multinomial = FALSE, ...) {
  ycat <- ifelse(
    nzchar(category),
    ifelse_(
      multinomial,
      category,
      paste0(y, "_", category)
    ),
    y
  )
  vectorize_beta <-
    K_fixed > 0L &&
    prior_distr$vectorized_beta &&
    prior_distr$beta_prior_npars > 0L
  vectorize_delta <-
    K_varying > 0L &&
    prior_distr$vectorized_delta &&
    prior_distr$delta_prior_npars > 0L
  vectorize_tau <-
    K_varying > 0L &&
    prior_distr$vectorized_tau &&
    prior_distr$tau_prior_npars > 0L
  vectorize_sigma_nu <-
    K_random > 0L &&
    prior_distr$vectorized_sigma_nu &&
    prior_distr$sigma_nu_prior_npars > 0L

  any_vectorized <-
    vectorize_beta ||
    vectorize_delta ||
    vectorize_tau ||
    vectorize_sigma_nu
  paste_rows(
    onlyif(
      any_vectorized,
      "// Parameters of vectorized priors"
    ),
    onlyif(
      vectorize_beta,
      paste0(
        "matrix[K_fixed_{y}, {prior_distr$beta_prior_npars}] ",
        "beta_prior_pars_{ycat};"
      )
    ),
    onlyif(
      vectorize_delta,
      paste0(
        "matrix[K_varying_{y}, {prior_distr$delta_prior_npars}] ",
        "delta_prior_pars_{ycat};"
      )
    ),
    onlyif(
      vectorize_tau,
      paste0(
        "matrix[K_varying_{y}, {prior_distr$tau_prior_npars}] ",
        "tau_prior_pars_{ycat};"
      )
    ),
    onlyif(
      vectorize_sigma_nu,
      paste0(
        "matrix[K_random_{y}, {prior_distr$sigma_nu_prior_npars}] ",
        "sigma_nu_prior_pars_{ycat};"
      )
    ),
    .indent = idt(1)
  )
}

data_lines_default <- function(y, idt, has_random_intercept,
                               K_fixed, K_varying, K_random, backend, ...) {
  icpt <- ifelse_(has_random_intercept, " - 1", "")
  re_icpt <- ifelse_(has_random_intercept, 1L, 0L)
  paste_rows(
    "// number of fixed, varying and random coefficients, and related indices",
    onlyif(K_fixed > 0L, "int<lower=0> K_fixed_{y};"),
    onlyif(K_varying > 0L, "int<lower=0> K_varying_{y};"),
    onlyif(
      K_random > 0L,
      "int<lower=0> K_random_{y}; // includes the potential random intercept"
    ),
    onlyif(
      K_fixed + K_varying > 0L,
      "int<lower=0> K_{y}; // K_fixed + K_varying"
    ),
    onlyif(
      K_fixed > 0L,
      stan_array(backend, "int", "J_fixed_{y}", "K_fixed_{y}")
    ),
    onlyif(
      K_varying > 0L,
      stan_array(backend, "int", "J_varying_{y}", "K_varying_{y}")
    ),
    onlyif(
      K_random > re_icpt,
      stan_array(
        backend,
        "int",
        "J_random_{y}",
        "K_random_{y}{icpt}",
        comment = "no intercept"
      )
    ),
    onlyif(
      K_fixed + K_varying > 0L,
      stan_array(
        backend,
        "int",
        "J_{y}",
        "K_{y}",
        comment = "fixed and varying"
      )
    ),
    onlyif(
      K_fixed > 0L,
      stan_array(backend, "int", "L_fixed_{y}", "K_fixed_{y}")
    ),
    onlyif(
      K_varying > 0L,
      stan_array(backend, "int", "L_varying_{y}", "K_varying_{y}")
    ),
    onlyif(
      K_random > re_icpt,
      stan_array(backend, "int", "L_random_{y}", "K_random_{y}{icpt}")
    ),
    .indent = idt(1)
  )
}

data_lines_categorical <- function(y, idt, default, has_missing,
                                   has_fully_missing, prior_distr,
                                   K_fixed, K_varying, K_random,
                                   categories, backend, ...) {
  prior_data <- lapply(
    categories[-1L],
    function(s) {
      prior_data_lines(
        y, idt, prior_distr[[s]], K_fixed, K_varying, K_random, category = s
      )
    }
  )
  paste_rows(
    "int<lower=0> S_{y}; // number of categories",
    missing_data_lines(y, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data,
    "// Response",
    stan_array(backend, "int", "y_{y}", "T, N", "lower=0"),
    .indent = idt(c(1, 0, 0, 0, 1, 1))
  )
}

data_lines_multinomial <- function(y_cg, idt, default, has_missing,
                                   has_fully_missing, prior_distr,
                                   K_fixed, K_varying, K_random,
                                   y, has_random_intercept, backend, ...) {

  default <- data_lines_default(
    y_cg, idt, has_random_intercept, K_fixed, K_varying, K_random, backend
  )
  prior_data <- ulapply(
    y[-1L],
    function(s) {
      prior_data_lines(
        y_cg, idt, prior_distr[[s]],
        K_fixed, K_varying, K_random,
        category = s, multinomial = TRUE
      )
    }
  )
  paste_rows(
    "int<lower=0> S_{y_cg}; // number of categories",
    missing_data_lines(y_cg, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data,
    "// Response",
    stan_array(backend, "int", "y_{y_cg}", "T, N, S_{y_cg}", "lower=0"),
    .indent = idt(c(1, 0, 0, 0, 1, 1))
  )
}

data_lines_gaussian <- function(y, idt, default, has_missing,
                                has_fully_missing, prior_distr,
                                K_fixed, K_varying, K_random, backend, ...) {
  paste_rows(
    missing_data_lines(y, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data_lines(y, idt, prior_distr, K_fixed, K_varying, K_random),
    "matrix[N, T] y_{y};",
    .indent = idt(c(0, 0, 0, 1))
  )
}

data_lines_mvgaussian <- function(y_cg, idt, default, has_missing,
                                  has_fully_missing, prior_distr,
                                  K_fixed, K_varying, K_random, backend, ...) {
  paste_rows(
    missing_data_lines(y_cg, idt, has_missing, has_fully_missing, backend),
    default,
    "int<lower=0> O_{y_cg};",
    stan_array(backend, "vector", "y_{y_cg}", "T, N", dims = "O_{y_cg}"),
    .indent = idt(c(0, 0, 1, 1))
  )
}

data_lines_binomial <- function(y, idt, default, has_missing,
                                has_fully_missing, prior_distr,
                                K_fixed, K_varying, K_random, backend, ...) {
  paste_rows(
    missing_data_lines(y, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data_lines(y, idt, prior_distr, K_fixed, K_varying, K_random),
    stan_array(backend, "int", "y_{y}", "T, N", "lower=0"),
    "// Trials for binomial response {y}",
    stan_array(backend, "int", "trials_{y}", "T, N", "lower=1"),
    .indent = idt(c(0, 0, 0, 1, 1, 1))
  )
}

data_lines_bernoulli <- function(y, idt, default, has_missing,
                                 has_fully_missing, prior_distr,
                                 K_fixed, K_varying, K_random,  backend, ...) {
  paste_rows(
    missing_data_lines(y, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data_lines(y, idt, prior_distr, K_fixed, K_varying, K_random),
    stan_array(backend, "int", "y_{y}", "T, N", "lower=0,upper=1"),
    .indent = idt(c(0, 0, 0, 1))
  )
}

data_lines_poisson <- function(y, idt, default, has_missing,
                               has_fully_missing, has_offset, prior_distr,
                               K_fixed, K_varying, K_random, backend, ...) {
  paste_rows(
    missing_data_lines(y, idt, has_missing,has_fully_missing, backend),
    default,
    prior_data_lines(y, idt, prior_distr, K_fixed, K_varying, K_random),
    stan_array(backend, "int", "y_{y}", "T, N", "lower=0"),
    onlyif(
      has_offset,
      "// Offset term"
    ),
    onlyif(
      has_offset,
      "matrix[N, T] offset_{y}"
    ),
    .indent = idt(c(0, 0, 0, 1, 1, 1))
  )
}

data_lines_negbin <- function(y, idt, default, has_missing,
                              has_fully_missing, has_offset, prior_distr,
                              K_fixed, K_varying, K_random, backend, ...) {
  paste_rows(
    missing_data_lines(y, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data_lines(y, idt, prior_distr, K_fixed, K_varying, K_random),
    stan_array(backend, "int", "y_{y}", "T, N", "lower=0"),
    onlyif(
      has_offset,
      "// Offset term"
    ),
    onlyif(
      has_offset,
      "matrix[N, T] offset_{y}"
    ),
    .indent = idt(c(0, 0, 0, 1, 1, 1))
  )
}

data_lines_exponential <- function(y, idt, default, has_missing,
                                   has_fully_missing, prior_distr,
                                   K_fixed, K_varying, K_random,
                                   backend, ...) {
  paste_rows(
    missing_data_lines(y, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data_lines(y, idt, prior_distr, K_fixed, K_varying, K_random),
    "matrix<lower=0>[N, T] y_{y};",
    .indent = idt(c(0, 0, 0, 1))
  )
}

data_lines_gamma <- function(y, idt, default, has_missing,
                             has_fully_missing, prior_distr,
                             K_fixed, K_varying, K_random, backend, ...) {
  paste_rows(
    missing_data_lines(y, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data_lines(y, idt, prior_distr, K_fixed, K_varying, K_random),
    "matrix<lower=0>[N, T] y_{y};",
    .indent = idt(c(0, 0, 0, 1))
  )
}

data_lines_beta <- function(y, idt, default, has_missing,
                            has_fully_missing, prior_distr,
                            K_fixed, K_varying, K_random, backend, ...) {
  paste_rows(
    missing_data_lines(y, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data_lines(y, idt, prior_distr, K_fixed, K_varying, K_random),
    "matrix<lower=0, upper=1>[N, T] y_{y};",
    .indent = idt(c(0, 0, 0, 1))
  )
}

data_lines_student <- function(y, idt, default, has_missing,
                               has_fully_missing, prior_distr,
                               K_fixed, K_varying, K_random, backend, ...) {
  paste_rows(
    missing_data_lines(y, idt, has_missing, has_fully_missing, backend),
    default,
    prior_data_lines(y, idt, prior_distr, K_fixed, K_varying, K_random),
    "matrix[N, T] y_{y};",
    .indent = idt(c(0, 0, 0, 1))
  )
}

# Transformed data block --------------------------------------------------

transformed_data_lines_default <- function(y, idt, ...) {
  list(
    declarations = "",
    statements = ""
  )
}

transformed_data_lines_categorical <- function(y, idt, K, S, ...) {
  list(
    declarations = "",
    statements = ""
  )
}

transformed_data_lines_multinomial <- function(y_cg, idt, K, S, ...) {
  list(
    declarations = "",
    statements = ""
  )
}

transformed_data_lines_gaussian <- function(default, ...) {
  default
}

transformed_data_lines_mvgaussian <- function(default, ...) {
  default[[1L]]
}

transformed_data_lines_binomial <- function(default, ...) {
  default
}

transformed_data_lines_bernoulli <- function(default, ...) {
  default
}

transformed_data_lines_poisson <- function(default, ...) {
  default
}

transformed_data_lines_negbin <- function(default, ...) {
  default
}

transformed_data_lines_exponential <- function(default, ...) {
  default
}

transformed_data_lines_gamma <- function(default, ...) {
  default
}

transformed_data_lines_beta <- function(default, ...) {
  default
}

transformed_data_lines_student <- function(default, ...) {
  default
}

# Parameters block --------------------------------------------------------

parameters_lines_default <- function(y, idt, noncentered, lb, has_fixed,
                                     has_varying, has_fixed_intercept,
                                     has_varying_intercept,
                                     has_lfactor, noncentered_psi,
                                     nonzero_lambda, ydim = y, ...) {

  oname <- ifelse_(noncentered, "omega_raw_", "omega_")

  # positivity constraint to deal with label-switching
  tr <- ifelse(has_lfactor && !nonzero_lambda, "<lower=0>", "")
  paste_rows(
    onlyif(
      has_fixed,
      "vector[K_fixed_{ydim}] beta_{y}; // Fixed coefficients"
    ),
    onlyif(
      has_varying,
      "matrix[K_varying_{ydim}, D] {oname}{y}; // Spline coefficients"
    ),
    onlyif(
      has_varying,
      paste0(
        "vector<lower={lb}>[K_varying_{ydim}] tau_{y}; ",
        "// SDs for the random walks"
      )
    ),
    ifelse_(
      has_fixed_intercept || has_varying_intercept,
      "real a_{y}; // Mean of the first time point",
      ""
    ),
    onlyif(
      has_varying_intercept,
      "row_vector[D - 1] omega_raw_alpha_{y}; // Coefficients for alpha"
    ),
    onlyif(
      has_varying_intercept,
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
      has_lfactor,
      paste0(
        "real{tr} omega_raw_psi_1_{y}; ",
        "// factor spline coef for first time point"
      )
    ),
    .indent = idt(1)
  )
}

parameters_lines_categorical <- function(y, idt, default, ...) {
  paste_rows(
    default,
    .parse = FALSE,
    .indent = idt(0)
  )
}

parameters_lines_multinomial <- function(y_cg, idt, univariate, ...) {
  paste_rows(
    univariate,
    .indent = idt(0)
  )
}


parameters_lines_gaussian <- function(y, idt, default, ...) {
  paste_rows(
    default,
    "real<lower=0> sigma_{y}; // SD of the normal distribution",
    .indent = idt(c(0, 1))
  )
}

parameters_lines_mvgaussian <- function(y_cg, idt, univariate, ...) {
  paste_rows(
    univariate,
    "cholesky_factor_corr[O_{y_cg}] L_{y_cg}; // Cholesky for gaussian",
    .indent = idt(c(0, 1))
  )
}

parameters_lines_binomial <- function(default, ...) {
  default
}

parameters_lines_bernoulli <- function(default, ...) {
  default
}

parameters_lines_poisson <- function(default, ...) {
  default
}

parameters_lines_negbin <- function(y, idt, default, ...) {
  paste_rows(
    default,
    "real<lower=0> phi_{y}; // Dispersion parameter of the NB distribution",
    .indent = idt(c(0, 1))
  )
}

parameters_lines_exponential <- function(y, idt, default, ...) {
  default
}

parameters_lines_gamma <- function(y, idt, default, ...) {
  paste_rows(
    default,
    "real<lower=0> phi_{y}; // Shape parameter of the Gamma distribution",
    .indent = idt(c(0, 1))
  )
}

parameters_lines_beta <- function(y, idt, default, ...) {
  paste_rows(
    default,
    "real<lower=0> phi_{y}; // Precision parameter of the Beta distribution",
    .indent = idt(c(0, 1))
  )
}

parameters_lines_student <- function(y, idt, default, ...) {
  paste_rows(
    default,
    "real<lower=0> sigma_{y}; // scale of the t-distribution",
    "real<lower=1> phi_{y}; // Degrees of freedom of the t-distribution",
    .indent = idt(c(0, 1, 1))
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
                                                 nonzero_lambda,
                                                 backend, ydim = y, ...) {
  if (noncentered) {
    xi_term <- ifelse_(shrinkage, " * xi[i - 1];", ";")
    declare_omega <- paste_rows(
      "matrix[K_varying_{ydim}, D] omega_{y};",
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
    stan_array(backend, "vector", "delta_{y}", "T", dims = "K_varying_{ydim}"),
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
  psi <- ifelse_(
    has_lfactor && nonzero_lambda,
    glue::glue(" - psi_{y}[1]"),
    ""
  )
  declare_fixed_intercept <- paste_rows(
    "// Time-invariant intercept",
    "real alpha_{y};",
    .indent = idt(1),
    .parse = FALSE
  )
  if (has_fixed || has_varying) {
    declare_omega_alpha_1 <- paste_rows(
      "// Time-varying intercept",
      stan_array(backend, "real", "alpha_{y}", "T"),
      "// Spline coefficients",
      "real omega_alpha_1_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha_1 <- paste_rows(
      "// Define the first alpha using mean a_{y}",
      "{{",
      "vector[K_{ydim}] gamma__{y};",
      onlyif(has_fixed, "gamma__{y}[L_fixed_{ydim}] = beta_{y};"),
      onlyif(has_varying, "gamma__{y}[L_varying_{ydim}] = delta_{y}[1];"),
      "omega_alpha_1_{y} = a_{y} - X_m[J_{ydim}] * gamma__{y};",
      "}}",
      .indent = idt(c(1, 1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
    state_fixed_intercept <- paste_rows(
      "// Define the first alpha using mean a_{y}",
      "{{",
      "vector[K_{ydim}] gamma__{y};",
      onlyif(has_fixed, "gamma__{y}[L_fixed_{ydim}] = beta_{y};"),
      onlyif(has_varying, "gamma__{y}[L_varying_{ydim}] = delta_{y}[1];"),
      "alpha_{y} = a_{y} - X_m[J_{ydim}] * gamma__{y}{psi};",
      "}}",
      .indent = idt(c(1, 1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
  } else {
    declare_omega_alpha_1 <- paste_rows(
      "// Time-invariant intercept",
      stan_array(backend, "real", "alpha_{y}", "T"),
      "real omega_alpha_1_{y} = a_{y};",
      .indent = idt(1),
      .parse = FALSE
    )
    state_omega_alpha_1 <- character(0L)
    state_fixed_intercept <- paste_rows(
      "alpha_{y} = a_{y}{psi};",
      .indent = idt(1),
      .parse = FALSE
    )
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
  declare_lambda <- paste_rows(
    "// hard sum constraint",
    paste0(
      "vector[N] lambda_{y} = {m}sigma_lambda_{y} * ",
      "sum_to_zero(lambda_raw_{y}, QR_Q);"
    ),
    .indent = idt(1),
    .parse = FALSE
  )
  if (noncentered_psi) {
    state_omega_psi <- paste_rows(
      "omega_psi_{y}[1] = omega_raw_psi_1_{y};",
      "omega_psi_{y} = cumulative_sum(omega_psi_{y});",
      .indent = idt(1),
      .parse = FALSE
    )
  } else {
    state_omega_psi <- paste_rows(
      "omega_psi_{y}[1] = omega_raw_psi_1_{y};",
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
  list(
    declarations = paste_rows(
      onlyif(has_lfactor, declare_psi),
      onlyif(has_lfactor, declare_lambda),
      onlyif(has_varying && noncentered, declare_omega),
      onlyif(has_varying, declare_delta),
      onlyif(has_fixed_intercept, declare_fixed_intercept),
      onlyif(has_varying_intercept, declare_varying_intercept),
      .indent = idt(0)
    ),
    statements = paste_rows(
      onlyif(has_lfactor, state_omega_psi),
      onlyif(has_lfactor, state_psi),
      onlyif(has_varying && noncentered, state_omega),
      onlyif(has_varying, state_delta),
      onlyif(has_fixed_intercept, state_fixed_intercept),
      onlyif(has_varying_intercept, state_varying_intercept),
      .indent = idt(0)
    )
  )
}

transformed_parameters_lines_categorical <- function(default, idt, ...) {
  list(
    declarations = paste_rows(
      ulapply(default, "[[", "declarations"),
      .parse = FALSE,
      .indent = idt(0)
    ),
    statements = paste_rows(
      ulapply(default, "[[", "statements"),
      .parse = FALSE,
      .indent = idt(0)
    )
  )
}

transformed_parameters_lines_multinomial <- function(default, idt, ...) {
  list(
    declarations = paste_rows(
      ulapply(default, "[[", "declarations"),
      .parse = FALSE,
      .indent = idt(0)
    ),
    statements = paste_rows(
      ulapply(default, "[[", "statements"),
      .parse = FALSE,
      .indent = idt(0)
    )
  )
}

transformed_parameters_lines_gaussian <- function(default, ...) {
  default
}

transformed_parameters_lines_mvgaussian <- function(default, idt, ...) {
  list(
    declarations = paste_rows(
      ulapply(default, "[[", "declarations"),
      .parse = FALSE,
      .indent = idt(0)
    ),
    statements = paste_rows(
      ulapply(default, "[[", "statements"),
      .parse = FALSE,
      .indent = idt(0)
    )
  )
}

transformed_parameters_lines_binomial <- function(default, ...) {
  default
}

transformed_parameters_lines_bernoulli <- function(default, ...) {
  default
}

transformed_parameters_lines_poisson <- function(default, ...) {
  default
}

transformed_parameters_lines_negbin <- function(default, ...) {
  default
}

transformed_parameters_lines_exponential <- function(default, ...) {
  default
}

transformed_parameters_lines_gamma <- function(default, ...) {
  default
}

transformed_parameters_lines_beta <- function(default, ...) {
  default
}

transformed_parameters_lines_student <- function(default, ...) {
  default
}

# Model block -------------------------------------------------------------

# priors which do not depend on family
prior_lines <- function(y, idt, noncentered, shrinkage,
                        has_varying, has_fixed, has_random,
                        has_fixed_intercept,
                        has_varying_intercept, has_random_intercept,
                        has_lfactor, noncentered_psi,
                        nonzero_lambda,
                        K_fixed, K_varying, K_random, prior_distr, ...) {
  if (prior_distr$vectorized_sigma_nu) {
    dpars_sigma_nu <- ifelse_(
      prior_distr$sigma_nu_prior_npars > 0L,
      paste0(
        "sigma_nu_prior_pars_", y, "[, ",
        seq_len(prior_distr$sigma_nu_prior_npars), "]",
        collapse = ", "
      ),
      ""
    )
    mtext_sigma_nu <- paste0(
      "sigma_nu_{y} ~ {prior_distr$sigma_nu_prior_distr}",
      "({dpars_sigma_nu});"
    )
  } else {
    mtext_sigma_nu <-
      "sigma_nu_{y}[{1:K_random}] ~ {prior_distr$sigma_nu_prior_distr};"
  }

  if (has_lfactor) {
    m <- ifelse_(nonzero_lambda, "1", "0")
    mtext_lambda <- paste_rows(
      "lambda_raw_{y} ~ normal(0, inv(sqrt(1 - inv(N))));",
      "sigma_lambda_{y} ~ {prior_distr$sigma_lambda_prior_distr};",
      onlyif(
        nonzero_lambda,
        "tau_psi_{y} ~ {prior_distr$tau_psi_prior_distr};"
      ),
      "omega_raw_psi_1_{y} ~ {prior_distr$psi_prior_distr};",
      .indent = idt(c(0, 1, 1, 1)),
      .parse = TRUE
    )
  }

  if (noncentered) {
    mtext_omega <- "omega_raw_alpha_{y} ~ std_normal();"
    if (prior_distr$vectorized_delta) {
      dpars_varying <- ifelse_(
        prior_distr$delta_prior_npars > 0L,
        paste0(
          "delta_prior_pars_", y, "[, ",
          seq_len(prior_distr$delta_prior_npars), "]",
          collapse = ", "
        ),
        ""
      )
      mtext_varying <- paste0(
        "omega_raw_{y}[, 1] ~ {prior_distr$delta_prior_distr}",
        "({dpars_varying});"
      )
    } else {
      mtext_varying <-
        "omega_raw_{y}[{1:K_varying}, 1] ~ {prior_distr$delta_prior_distr};"
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
    if (prior_distr$vectorized_delta) {
      dpars_varying <- ifelse_(
        prior_distr$delta_prior_npars > 0L,
        paste0(
          "delta_prior_pars_", y, "[, ",
          seq_len(prior_distr$delta_prior_npars), "]",
          collapse = ", "
        ),
        ""
      )
      mtext_varying <- paste0(
        "omega_{y}[, 1] ~ {prior_distr$delta_prior_distr}",
        "({dpars_varying});"
      )
    } else {
      mtext_varying <-
        "omega_{y}[{1:K_varying}, 1] ~ {prior_distr$delta_prior_distr};"
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

  mtext_alpha <- "a_{y} ~ {prior_distr$alpha_prior_distr};"
  mtext_varying_intercept <- paste_rows(
    mtext_alpha,
    mtext_omega,
    "tau_alpha_{y} ~ {prior_distr$tau_alpha_prior_distr};",
    .indent = idt(c(0, 1, 1)),
    .parse = FALSE
  )
  mtext_fixed_intercept <- mtext_alpha

  if (prior_distr$vectorized_beta) {
    dpars_fixed <-  ifelse_(
      prior_distr$beta_prior_npars > 0L,
      paste0(
        "beta_prior_pars_", y, "[, ",
        seq_len(prior_distr$beta_prior_npars), "]",
        collapse = ", "
      ),
      ""
    )
    mtext_fixed <- "beta_{y} ~ {prior_distr$beta_prior_distr}({dpars_fixed});"
  } else {
    mtext_fixed <- "beta_{y}[{1:K_fixed}] ~ {prior_distr$beta_prior_distr};"
  }

  if (prior_distr$vectorized_tau) {
    dpars_tau <- ifelse_(
      prior_distr$tau_prior_npars > 0L,
      paste0(
        "tau_prior_pars_", y, "[, ", seq_len(prior_distr$tau_prior_npars), "]",
        collapse = ", "
      ),
      ""
    )
    mtext_tau <- "tau_{y} ~ {prior_distr$tau_prior_distr}({dpars_tau});"
  } else {
    mtext_tau <- "tau_{y}[{1:K_varying}] ~ {prior_distr$tau_prior_distr};"
  }
  paste_rows(
    onlyif(has_lfactor, mtext_lambda),
    onlyif(has_random_intercept || has_random, mtext_sigma_nu),
    onlyif(has_fixed_intercept, mtext_fixed_intercept),
    onlyif(has_varying_intercept, mtext_varying_intercept),
    onlyif(has_fixed, mtext_fixed),
    onlyif(has_varying, mtext_varying),
    onlyif(has_varying, mtext_tau),
    .indent = idt(c(1, 1, 1, 1, 1, 1, 1))
  )
}

loglik_fun_args <- function(y, has_fixed, has_varying, has_missing,
                            has_fixed_intercept, has_varying_intercept,
                            has_random_intercept, has_random, has_lfactor,
                            glm) {
  has_X <- has_fixed || has_varying || has_random
  LJ <- ifelse_(
    glm,
    "L",
    "J"
  )
  glue::glue(cs(c(
    ifelse_(
      has_missing,
      c("obs_{y}", "n_obs_{y}"),
      "N"
    ),
    "y_{y}",
    onlyif(has_fixed_intercept || has_varying_intercept, "alpha_{y}"),
    onlyif(has_random, c("L_random_{y}", "K_random_{y}")),
    onlyif(has_random || has_random_intercept, "nu_{y}"),
    onlyif(has_lfactor, c("lambda_{y}", "psi_{y}")),
    onlyif(has_fixed, c("{LJ}_fixed_{y}", "beta_{y}")),
    onlyif(has_varying, c("{LJ}_varying_{y}", "delta_{y}")),
    onlyif(has_fixed || has_varying, c("J_{y}", "K_{y}")),
    onlyif(has_X, "X")
  )))
}

model_lines_categorical <- function(y, obs, idt, priors,
                                    has_missing, has_fully_missing,
                                    has_fixed_intercept, has_varying_intercept,
                                    has_random_intercept,
                                    has_fixed, has_varying, has_random,
                                    has_lfactor, threading,
                                    categories, multinomial = FALSE,
                                    backend, ...) {
  distr <- ifelse_(multinomial, "multinomial", "categorical")
  S <- length(categories)
  cats <- categories[seq.int(2L, S)]
  fun_args <- character(length(y))
  for (i in seq_along(cats)) {
    yi <- ifelse_(
      multinomial,
      cats[i],
      paste0(y, "_", cats[i])
    )
    fun_args[i] <- glue::glue(cs(c(
      onlyif(has_fixed_intercept || has_varying_intercept, "alpha_{yi}"),
      onlyif(has_random || has_random_intercept, "nu_{yi}"),
      onlyif(has_lfactor, c("lambda_{y}_{yi}", "psi_{yi}")),
      onlyif(has_fixed, "beta_{yi}"),
      onlyif(has_varying, "delta_{yi}")
    )))
  }
  common_intercept <- !has_random && !has_random_intercept && !has_lfactor
  glm <- stan_supports_categorical_logit_glm(backend, common_intercept) &&
    (has_fixed || has_varying) && !multinomial
  LJ <- ifelse_(
    glm,
    "L",
    "J"
  )
  has_X <- has_fixed || has_varying || has_random
  fun_args <- glue::glue(cs(c(
    ifelse_(
      has_missing,
      c("obs_{y}", "n_obs_{y}"),
      "N"
    ),
    "y_{y}",
    "S_{y}",
    fun_args,
    onlyif(has_random, c("L_random_{y}", "K_random_{y}")),
    onlyif(has_fixed, "{LJ}_fixed_{y}"),
    onlyif(has_varying, "{LJ}_varying_{y}"),
    onlyif(has_fixed || has_varying, c("J_{y}", "K_{y}")),
    onlyif(has_X, "X")
  )))
  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum({distr}_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += {distr}_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)

}

model_lines_multinomial <- function(cvars, cgvars, idt, threading, ...) {
  cgvars$priors <- lapply(
    cgvars$y[-1L],
    function(s) {
      cvars[[s]]$y <- s
      cvars[[s]]$prior_distr <- cgvars$prior_distr[[s]]
      do.call(prior_lines, c(cvars[[s]], idt = idt))
    }
  )
  cgvars$categories <- cgvars$y
  cgvars$y <- cgvars$y_cg
  cgvars$multinomial <- TRUE
  do.call(
    model_lines_categorical,
    args = c(cgvars, idt = idt, threading = threading)
  )
}

model_lines_gaussian <- function(y, obs, idt, priors,
                                 has_missing, has_fully_missing,
                                 has_fixed_intercept, has_varying_intercept,
                                 has_random_intercept,
                                 has_fixed, has_varying, has_random,
                                 has_lfactor, threading, prior_distr, ...) {

  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  fun_args <- paste0(
    c(
      loglik_fun_args(
        y, has_fixed, has_varying, has_missing, has_fixed_intercept,
        has_varying_intercept, has_random_intercept, has_random, has_lfactor,
        glm = TRUE
      ),
      glue::glue("sigma_{y}")
    ), collapse = ", "
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(gaussian_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += gaussian_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    "sigma_{y} ~ {prior_distr$sigma_prior_distr};",
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

model_lines_mvgaussian <- function(cvars, cgvars, idt, threading, ...) {
  y <- cgvars$y
  y_cg <- cgvars$y_cg
  priors <- fun_args <- character(length(y))
  has_X <- FALSE
  for (i in seq_along(y)) {
    yi <- y[i]
    args <- c(cvars[[i]], idt = idt)
    priors[i] <- do.call(prior_lines, args = args)
    priors[i] <- paste_rows(
      priors[i],
      glue::glue("sigma_{y[i]} ~ {cvars[[i]]$sigma_prior_distr};"),
      .indent = idt(c(0, 1)),
      .parse = FALSE
    )
    has_X <- has_X ||
      cvars[[i]]$has_fixed || cvars[[i]]$has_varying || cvars[[i]]$has_random
    fun_args[i] <- glue::glue(cs(c(
      onlyif(
        cvars[[i]]$has_fixed_intercept || cvars[[i]]$has_varying_intercept,
        "alpha_{yi}"
      ),
      onlyif(cvars[[i]]$has_random, c("L_random_{yi}", "K_random_{yi}")),
      onlyif(
        cvars[[i]]$has_random || cvars[[i]]$has_random_intercept,
        "nu_{yi}"
      ),
      onlyif(cvars[[i]]$has_lfactor, c("lambda_{yi}", "psi_{yi}")),
      onlyif(cvars[[i]]$has_fixed, c("J_fixed_{yi}", "beta_{yi}")),
      onlyif(cvars[[i]]$has_varying, c("J_varying_{yi}", "delta_{yi}")),
      "sigma_{yi}"
    )))
  }

  fun_args <- glue::glue(cs(c(
    ifelse_(
      cgvars$has_missing,
      c("obs_{y_cg}", "n_obs_{y_cg}"),
      "N"
    ),
    "y_{y_cg}",
    "O_{y_cg}",
    fun_args,
    onlyif(has_X, "X"),
    "L_{y_cg}"
  )))

  seq1T <- ifelse_(
    cgvars$has_fully_missing,
    "t_obs_{y_cg}",
    "seq1T"
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(mvgaussian_loglik_{y_cg}_lpmf, {seq1T}, grainsize,",
      " {fun_args});"
    ),
    glue::glue(
      "target += mvgaussian_loglik_{y_cg}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    "L_{y_cg} ~ {cgvars$prior_distr$L_prior_distr};",
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

model_lines_bernoulli <- function(y, obs, idt, priors,
                                  has_missing, has_fully_missing,
                                  has_fixed_intercept, has_varying_intercept,
                                  has_random_intercept,
                                  has_fixed, has_varying, has_random,
                                  has_lfactor, threading, ...) {


  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  fun_args <- loglik_fun_args(
    y, has_fixed, has_varying, has_missing, has_fixed_intercept,
    has_varying_intercept, has_random_intercept, has_random, has_lfactor,
    glm = TRUE
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(bernoulli_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += bernoulli_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

model_lines_poisson <- function(y, obs, idt, priors,
                                has_missing, has_fully_missing,
                                has_fixed_intercept, has_varying_intercept,
                                has_random_intercept,
                                has_fixed, has_varying, has_random,
                                has_lfactor, threading, has_offset, ...) {
  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  fun_args <- paste0(c(
    loglik_fun_args(
      y, has_fixed, has_varying, has_missing, has_fixed_intercept,
      has_varying_intercept, has_random_intercept, has_random, has_lfactor,
      glm = TRUE
    ),
    onlyif(has_offset, glue::glue("offset_{y}"))),
    collapse = ", "
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(poisson_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += poisson_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

model_lines_negbin <- function(y, obs, idt, priors,
                               has_missing, has_fully_missing,
                               has_fixed_intercept, has_varying_intercept,
                               has_random_intercept,
                               has_fixed, has_varying, has_random,
                               has_lfactor, threading, prior_distr,
                               has_offset, ...) {
  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  fun_args <- paste0(
    c(
      loglik_fun_args(
        y, has_fixed, has_varying, has_missing, has_fixed_intercept,
        has_varying_intercept, has_random_intercept, has_random, has_lfactor,
        glm = TRUE
      ),
      onlyif(has_offset, glue::glue("offset_{y}")),
      "phi_{y}"
    ),
    collapse = ","
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(negbin_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += negbin_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    "phi_{y} ~ {prior_distr$phi_prior_distr};",
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

model_lines_binomial <- function(y, obs, idt, priors,
                                 has_missing, has_fully_missing,
                                 has_fixed_intercept, has_varying_intercept,
                                 has_random_intercept,
                                 has_fixed, has_varying, has_random,
                                 has_lfactor, threading, ...) {
  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  fun_args <- cs(
    c(
      loglik_fun_args(
        y, has_fixed, has_varying, has_missing, has_fixed_intercept,
        has_varying_intercept, has_random_intercept, has_random, has_lfactor,
        glm = FALSE
      ),
      glue::glue("trials_{y}")
    )
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(binomial_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += binomial_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

model_lines_exponential <- function(y, obs, idt, priors,
                                    has_missing, has_fully_missing,
                                    has_fixed_intercept, has_varying_intercept,
                                    has_random_intercept,
                                    has_fixed, has_varying, has_random,
                                    has_lfactor, threading, ...) {
  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  fun_args <- loglik_fun_args(
    y, has_fixed, has_varying, has_missing, has_fixed_intercept,
    has_varying_intercept, has_random_intercept, has_random, has_lfactor,
    glm = FALSE
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(exponential_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += exponential_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

model_lines_gamma <- function(y, obs, idt, priors,
                              has_missing, has_fully_missing,
                              has_fixed_intercept, has_varying_intercept,
                              has_random_intercept,
                              has_fixed, has_varying, has_random,
                              has_lfactor, threading, prior_distr, ...) {
  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  fun_args <- paste0(
    c(
      loglik_fun_args(
        y, has_fixed, has_varying, has_missing, has_fixed_intercept,
        has_varying_intercept, has_random_intercept, has_random, has_lfactor,
        glm = FALSE
      ),
      glue::glue("phi_{y}")
    ),
    collapse = ", "
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(gamma_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += gamma_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    "phi_{y} ~ {prior_distr$phi_prior_distr};",
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

model_lines_beta <- function(y, obs, idt, priors,
                             has_missing, has_fully_missing,
                             has_fixed_intercept, has_varying_intercept,
                             has_random_intercept,
                             has_fixed, has_varying, has_random,
                             has_lfactor, threading, prior_distr, ...) {
  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  fun_args <- paste0(
    c(
      loglik_fun_args(
        y, has_fixed, has_varying, has_missing, has_fixed_intercept,
        has_varying_intercept, has_random_intercept, has_random, has_lfactor,
        glm = FALSE
      ),
      glue::glue("phi_{y}")
    ),
    collapse = ", "
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(beta_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += beta_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    "phi_{y} ~ {prior_distr$phi_prior_distr};",
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

model_lines_student <- function(y, obs, idt, priors,
                                has_missing, has_fully_missing,
                                has_fixed_intercept, has_varying_intercept,
                                has_random_intercept,
                                has_fixed, has_varying, has_random,
                                has_lfactor, threading, prior_distr, ...) {
  seq1T <- ifelse_(
    has_fully_missing,
    "t_obs_{y}",
    "seq1T"
  )
  fun_args <- paste0(
    c(
      loglik_fun_args(
        y, has_fixed, has_varying, has_missing, has_fixed_intercept,
        has_varying_intercept, has_random_intercept, has_random, has_lfactor,
        glm = FALSE
      ),
      glue::glue("sigma_y"),
      glue::glue("phi_{y}")
    ),
    collapse = ", "
  )
  likelihood <- ifelse_(
    threading,
    glue::glue(
      "target += reduce_sum(student_loglik_{y}_lpmf, {seq1T}, grainsize, ",
      "{fun_args});"
    ),
    glue::glue(
      "target += student_loglik_{y}_lpmf({seq1T} | {fun_args});"
    )
  )
  model_text <- paste_rows(
    "sigma_{y} ~ {prior_distr$sigma_prior_distr};",
    "phi_{y} ~ {prior_distr$phi_prior_distr};",
    likelihood,
    .indent = idt(1)
  )
  paste_rows(priors, model_text, .parse = FALSE)
}

# Generated quantities block ----------------------------------------------

generated_quantities_lines_default <- function() {
  ""
}

generated_quantities_lines_categorical <- function(...) {
  ""
}

generated_quantities_lines_multinomial <- function(...) {
  ""
}

generated_quantities_lines_gaussian <- function(...) {
  ""
}

generated_quantities_lines_mvgaussian <- function(y, y_cg, idt, ...) {
  O <- length(y)
  paste_rows(
    "matrix[O_{y_cg},O_{y_cg}] corr_matrix_{y_cg} = ",
    "multiply_lower_tri_self_transpose(L_{y_cg});",
    "vector[{(O * (O - 1L)) %/% 2L}] corr_{y_cg};",
    "for (k in 1:O_{y_cg}) {{",
    "for (j in 1:(k - 1)) {{",
    "corr_{y_cg}[choose(k - 1, 2) + j] = corr_matrix_{y_cg}[j, k];",
    "}}",
    "}}",
    .indent = idt(c(1, 2, 1, 1, 2, 3, 2, 1))
  )
}

generated_quantities_lines_binomial <- function(...) {
  ""
}

generated_quantities_lines_bernoulli <- function(...) {
  ""
}

generated_quantities_lines_poisson <- function(...) {
  ""
}

generated_quantities_lines_negbin <- function(...) {
  ""
}

generated_quantities_lines_exponential <- function(...) {
  ""
}

generated_quantities_lines_gamma <- function(...) {
  ""
}

generated_quantities_lines_beta <- function(...) {
  ""
}

generated_quantities_lines_student <- function(...) {
  ""
}

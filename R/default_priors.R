#' Create Default Priors for Non-categorical Data
#'
#' @param y \[`character(1)`]\cr Name of the response variable of the channel.
#' @param channel \[`list()`]\cr Channel-specific helper variables.
#' @param mean_gamma Prior mean betas and deltas (at time `fixed + 1`).
#' @param sd_gamma Prior SD betas and deltas (at time `fixed + 1`).
#' @param mean_y Mean of the response variable at time `fixed + 1`.
#' @param sd_y Standard deviation of the response variable at time `fixed + 1`.
#' @noRd
default_priors <- function(y, channel, mean_gamma, sd_gamma, mean_y, sd_y) {
  mean_y <- signif(mean_y, 2)
  sd_y <- signif(sd_y, 2)
  mean_gamma <- signif(mean_gamma, 2)
  sd_gamma <- signif(sd_gamma, 2)
  priors <- list()
  if (channel$has_random_intercept || channel$has_random) {
    icpt <- ifelse_(
      channel$has_random_intercept,
      "alpha",
      NULL
    )
    ns <- c(icpt, names(sd_gamma[channel$J_random]))
    channel$sigma_nu_prior_distr <- "normal"
    channel$sigma_nu_prior_npars <- 2L
    channel$sigma_nu_prior_pars <- cbind(0, rep(1, channel$K_random))
    priors$sigma_nu <- data.frame(
      parameter = paste0("sigma_nu_", y, "_", ns),
      response = y,
      prior = "normal(0, 1)",
      type = "sigma_nu",
      category = ""
    )
  }
  if (channel$has_lfactor) {
    channel$sigma_lambda_prior_distr <- "normal(0, 1)"
    priors$sigma_lambda <- data.frame(
      parameter = paste0("sigma_lambda_", y),
      response = y,
      prior = channel$sigma_lambda_prior_distr,
      type = "sigma_lambda",
      category = ""
    )
    if (channel$nonzero_lambda) {
      channel$tau_psi_prior_distr <- "normal(0, 1)"
      priors$tau_psi <- data.frame(
        parameter = paste0("tau_psi_", y),
        response = y,
        prior = channel$tau_psi_prior_distr,
        type = "tau_psi",
        category = ""
      )
    }
    channel$psi_prior_distr <- "normal(0, 1)"
    priors$psi <- data.frame(
      parameter = paste0("psi_", y),
      response = y,
      prior = channel$psi_prior_distr,
      type = "psi",
      category = ""
    )
  }
  if (channel$has_fixed_intercept || channel$has_varying_intercept) {
    channel$alpha_prior_distr <- paste0("normal(", mean_y, ", ", 2 * sd_y, ")")
    priors$alpha <- data.frame(
      parameter = paste0("alpha_", y),
      response = y,
      prior = channel$alpha_prior_distr,
      type = "alpha",
      category = ""
    )
    if (channel$has_varying_intercept) {
      channel$tau_alpha_prior_distr <- "normal(0, 1)"
      priors$tau_alpha <- data.frame(
        parameter = paste0("tau_alpha_", y),
        response = y,
        prior = "normal(0, 1)",
        type = "tau_alpha",
        category = ""
      )
    }
  }
  if (channel$has_fixed) {
    m <- mean_gamma[channel$J_fixed]
    s <- sd_gamma[channel$J_fixed]
    channel$beta_prior_distr <- "normal"
    channel$beta_prior_npars <- 2L
    channel$beta_prior_pars <- unname(cbind(m, s))
    priors$beta <- data.frame(
      parameter = paste0("beta_", y, "_", names(s)),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "beta",
      category = ""
    )
  }
  if (channel$has_varying) {
    m <- mean_gamma[channel$J_varying]
    s <- sd_gamma[channel$J_varying]
    channel$delta_prior_distr <- "normal"
    channel$delta_prior_npars <- 2L
    channel$delta_prior_pars <- unname(cbind(m, s))
    priors$delta <- data.frame(
      parameter = paste0("delta_", y, "_", names(s)),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "delta",
      category = ""
    )
    channel$tau_prior_distr <- "normal"
    channel$tau_prior_npars <- 2
    channel$tau_prior_pars <- cbind(0, rep(1, channel$K_varying))
    priors$tau <- data.frame(
      parameter = paste0("tau_", y, "_", names(s)),
      response = y,
      prior = "normal(0, 1)",
      type = "tau",
      category = ""
    )
  }
  list(
    channel = channel,
    priors = data.table::setDF(data.table::rbindlist(priors))
  )
}

#' Create Default Priors for Categorical Data
#'
#' @param y \[`character(1)`]\cr Name of the response variable of the channel.
#' @param channel \[`list()`]\cr Channel-specific helper variables.
#' @param sd_x \[`numeric(1)`]\cr
#'   Standard deviation of the explanatory variables at time `fixed + 1`.
#' @param resp_class \[`character(1)`]\cr Class of the response variable.
#' @noRd
default_priors_categorical <- function(y, channel, sd_x, resp_class) {
  S_y <- length(attr(resp_class, "levels"))
  # remove the first level which acts as reference
  resp_levels <- attr(resp_class, "levels")[-1]
  sd_gamma <- signif(2 / sd_x, 2)
  priors <- list()
  if (channel$has_fixed_intercept || channel$has_varying_intercept) {
    m <- rep(0.0, S_y - 1L)
    s <- rep(2.0, S_y - 1L)
    channel$alpha_prior_distr <- "normal"
    channel$alpha_prior_npars <- 2L
    channel$alpha_prior_pars <- unname(cbind(m, s))
    priors$alpha <- data.frame(
      parameter = paste0("alpha_", y),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "alpha",
      category = resp_levels
    )
    if (channel$has_varying_intercept) {
      channel$tau_alpha_prior_distr <- "normal(0, 1)"
      priors$tau_alpha <- data.frame(
        parameter = paste0("tau_alpha_", y),
        response = y,
        prior = "normal(0, 1)",
        type = "tau_alpha",
        category = ""
      )
    }
  }
  if (channel$has_fixed) {
    m <- rep(0.0, channel$K_fixed * (S_y - 1L))
    s <- rep(sd_gamma[channel$J_fixed], S_y - 1L)
    channel$beta_prior_distr <- "normal"
    channel$beta_prior_npars <- 2L
    channel$beta_prior_pars <- unname(cbind(m, s))
    priors$beta <- data.frame(
      parameter = paste0("beta_", y, "_", names(s)),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "beta",
      category = rep(resp_levels, each = channel$K_fixed)
    )
  }
  if (channel$has_varying) {
    m <- rep(0.0, channel$K_varying * (S_y - 1L))
    s <- rep(sd_gamma[channel$J_varying], S_y - 1L)
    channel$delta_prior_distr <- "normal"
    channel$delta_prior_npars <- 2L
    channel$delta_prior_pars <- unname(cbind(m, s))
    priors$delta <- data.frame(
      parameter = paste0("delta_", y, "_", names(s)),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "delta",
      category = rep(resp_levels, each = channel$K_varying)
    )
    channel$tau_prior_distr <- "normal"
    channel$tau_prior_npars <- 2L
    channel$tau_prior_pars <- cbind(0.0, rep(1.0, channel$K_varying))
    priors$tau <- data.frame(
      parameter = paste0("tau_", y, "_", names(sd_gamma[channel$J_varying])),
      response = y,
      prior = "normal(0, 1)",
      type = "tau",
      category = ""
    )
  }
  list(
    channel = channel,
    priors = data.table::setDF(data.table::rbindlist(priors))
  )
}

#' Check and Correct the User-defined Priors
#'
#' This function makes a crude check that the user-supplied prior distributions
#' are valid. The actual syntax is later tested automatically during
#' compilation, this is mainly for checking the support of the distribution,
#' so that the users don't supply constrained priors for coefficients, and for
#' ordering the corresponding data frame coherently.
#'
#' @param priors A data frame of prior definitions.
#' @param defaults A data frame of default prior definitions.
#' @noRd
check_priors <- function(priors, defaults) {
  not_found <-
    defaults$parameter[which(!(defaults$parameter %in% priors$parameter))]
  not_found_len <- length(not_found)
  stopifnot_(
    identical(not_found_len, 0L),
    c(
      "Argument {.var priors} must contain all relevant parameters:",
      `x` = "{cli::qty(not_found_len)} Prior{?s} for parameter{?s}
             {.var {not_found}} {?is/are} not defined."
    )
  )
  extras <-
    priors$parameter[which(!(priors$parameter %in% defaults$parameter))]
  extras_len <- length(extras)
  stopifnot_(
    identical(extras_len, 0L),
    c(
      "Argument {.var priors} must contain only relevant parameters:",
      `x` = "{cli::qty(extras)} Found {?a/} prior{?s} for parameter{?s}
             {.var {extras}} but the model does not contain such
             {?a/} parameter{?s}."
    )
  )
  # order to match the code generation
  prior_order <- match(priors$parameter, defaults$parameter)
  priors <- priors[order(prior_order), ]
  unconstrained_dists <- c(
    "normal", "student_t", "double_exponential", "cauchy", "exp_mod_normal",
    "skew_normal", "logistic", "gumbel", "skew_double_exponential"
  )
  positive_dists <- c(
    "gamma", "exponential", "lognormal", "chi_square", "inv_chi_square",
    "scaled_inv_chi_square", "inv_gamma", "weibull", "frechet", "rayleigh"
  )
  all_dists <- c(unconstrained_dists, positive_dists)
  dists <- sub("\\(.*", "", priors$prior)
  unsupported <- unique(dists[which(!(dists %in% all_dists))])
  unsupported_len <- length(unsupported)
  stopifnot_(
    identical(unsupported_len, 0L),
    c(
      "{cli::qty(unsupported_len)} Found {?an/} unsupported prior
       distribution{?s} in {.var priors}:",
      `x` = "Distribution{?s} {.var {unsupported}} {?is/are} not available."
    )
  )
  unsupported <- which(
    priors$type %in% c("alpha", "beta", "delta") &
      !(dists %in% unconstrained_dists)
  )
  unsupported_len <- length(unsupported)
  pars <- priors$parameter[unsupported]
  dists <- dists[unsupported]
  stopifnot_(
    identical(unsupported_len, 0L),
    c(
      "Priors for parameters alpha, beta, and delta should have unconstrained
      support:",
      `x` = "{cli::qty(unsupported_len)} Found {?an/} unconstrained
             distribution{?s} {.var {dists}} for parameter{?s} {.var {pars}}."
    )
  )
  priors
}

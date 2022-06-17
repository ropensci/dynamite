#' Create Default Priors for Non-categorical Data
#'
#' @param y \[`character(1)`]\cr Name of the response variable of the channel.
#' @param channel \[`list()`]\cr Channel-specific helper variables.
#' @param mean_gamma Prior mean betas and deltas (at time fixed + 1).
#' @param sd_gamma Prior SD betas and deltas (at time fixed + 1).
#' @param mean_y Mean of the response variable at time fixed + 1.
#' @param sd_y SD of the response variable at time fixed + 1.
#' @noRd
default_priors <- function(y, channel, mean_gamma, sd_gamma,
                           mean_y, sd_y) {

  priors <- list()
  if (channel$has_random_intercept) {
    channel$sigma_nu_prior_distr <-  paste0("normal(", 0, ", ", sd_y, ")")
    priors$sigma_nu <- data.frame(
      parameter = paste0("sigma_nu_", y),
      response = y,
      prior = channel$sigma_nu_prior_distr,
      type = "sigma_nu",
      category = ""
    )
  }
  if (channel$has_fixed_intercept || channel$has_varying_intercept) {
    channel$alpha_prior_distr <-  paste0("normal(", mean_y, ", ", 2 * sd_y, ")")
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
    channel$beta_prior_npars <- 2
    channel$beta_prior_pars <- cbind(m, s, deparse.level = 0)
    channel$beta_prior_distr <- "normal"
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
    channel$delta_prior_npars <- 2
    channel$delta_prior_pars <- cbind(m, s, deparse.level = 0)
    channel$delta_prior_distr <- "normal"
    priors$delta <- data.frame(
      parameter = paste0("delta_", y, "_", names(s)),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "delta",
      category = ""
    )
    channel$tau_prior_npars <- 2
    channel$tau_prior_pars <- cbind(0, rep(1, channel$K_varying))
    channel$tau_prior_distr <- "normal"
    priors$tau <- data.frame(
      parameter = paste0("tau_", y, "_", names(s)),
      response = y,
      prior = "normal(0, 1)",
      type = "tau",
      category = ""
    )
  }
  list(channel = channel, priors = dplyr::bind_rows(priors))
}
#' Create Default Priors for Categorical Data
#'
#' @param y \[`character(1)`]\cr Name of the response variable of the channel.
#' @param channel \[`list()`]\cr Channel-specific helper variables.
#' @param sd_x Standard deviation of the explanatory variables at time
#'   fixed + 1.
#' @param resp_class Class of the response variable.
#' @noRd
default_priors_categorical <- function(y, channel, sd_x, resp_class) {
  S_y <- length(attr(resp_class, "levels"))
  # remove the first level which acts as reference
  resp_levels <- attr(resp_class, "levels")[-1]
  sd_gamma <- 2 / sd_x
  priors <- list()

  if (channel$has_fixed_intercept || channel$has_varying_intercept) {
    m <- rep(0, S_y - 1)
    s <- rep(2, S_y - 1)
    channel$alpha_prior_npars <- 2
    channel$alpha_prior_pars <- cbind(m, s, deparse.level = 0)
    channel$alpha_prior_distr <- "normal"
    priors$alpha <- data.frame(
      parameter = paste0("alpha_", y),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "alpha",
      category = resp_levels
    )
    if (channel$has_varying_intercept) {
      channel$tau_alpha_prior_distr <- "normal(0, 1);"
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
    m <- rep(0, channel$K_fixed * (S_y - 1))
    s <- rep(sd_gamma[channel$J_fixed], S_y - 1)
    channel$beta_prior_npars <- 2
    channel$beta_prior_distr <- "normal"
    channel$beta_prior_pars <- cbind(m, s, deparse.level = 0)
    priors$beta <- data.frame(
      parameter = paste0("beta_", y, "_", names(s)),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "beta",
      category = rep(resp_levels, each = channel$K_fixed)
    )
  }
  if (channel$has_varying) {
    m <- rep(0, channel$K_varying * (S_y - 1))
    s <- rep(sd_gamma[channel$J_varying], S_y - 1)
    channel$delta_prior_npars <- 2
    channel$delta_prior_pars <- cbind(m, s, deparse.level = 0)
    channel$delta_prior_distr <- "normal"
    priors$delta <- data.frame(
      parameter = paste0("delta_", y, "_", names(s)),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "delta",
      category = rep(resp_levels, each = channel$K_varying)
    )
    channel$tau_prior_npars <- 2
    channel$tau_prior_pars <- cbind(0, rep(1, channel$K_varying))
    channel$tau_prior_distr <- "normal"
    priors$tau <- data.frame(
      parameter = paste0("tau_", y, "_", names(s)),
      response = y,
      prior = "normal(0, 1)",
      type = "tau",
      category = ""
    )
  }
  list(channel = channel, priors = dplyr::bind_rows(priors))
}
#' Check and Correct the User-defined Priors
#'
#' This function makes a crude check that the user-supplied prior distributions
#' are valid. The actual syntax is later tested automatically during
#' compilation, this is mainly for checking the support of the distribution,
#' so that the users don't supply constrained priors for coefficients, and for
#' ordering the corresponding data frame coherently.
#'
#' @noRd
check_priors <- function(priors, defaults) {

  not_found <-
    defaults$parameter[which(!(defaults$parameter %in% priors$parameter))]
  if (length(not_found) > 0) {
    stop_(c(
      "Argument `priors` should contain all relevant parameters.",
      `x` ="Prior{?s} for parameter{?s} {not_found} {?is/are} not defined."))
  }
  extras <-
    priors$parameter[which(!(priors$parameter %in% defaults$parameter))]
  if (length(extras) > 0) {
    stop_(c(
      "Argument `priors` should contain only relevant parameters.",
      `x` = "Found {?a/} prior{?s} for parameter{?s} {.var {extras}}, but
        the model does not contain such {?a/} parameter{?s}."))
  }
  # order to match the code generation
  priors <- priors |>
    dplyr::arrange(match(.data$parameter, defaults$parameter))

  unconstrained_dists <- c(
    "normal", "student_t", "double_exponential", "cauchy", "exp_mod_normal",
    "skew_normal", "logistic","gumbel", "skew_double_exponential"
  )
  positive_dists <- c(
    "gamma", "exponential", "lognormal", "chi_square", "inv_chi_square",
    "scaled_inv_chi_square", "inv_gamma", "weibull", "frechet", "rayleigh"
  )
  all_dists <- c(unconstrained_dists, positive_dists)
  dists <- sub("\\(.*", "", priors$prior)
  unsupported <- unique(dists[which(!(dists %in% all_dists))])
  if (length(unsupported) > 0) {
    stop_(c("Found unsupported prior distributions in `priors:",
            `x` = "Distribution{?s} {.var {unsupported}}, are not available."))
  }
  unsupported <- which(
    priors$type %in% c("alpha", "beta", "delta") &
      !(dists %in% unconstrained_dists))
  if (length(unsupported) > 0) {
    pars <- priors$parameter[unsupported]
    dists <- dists[unsupported]
    stop_(c(
      "Priors for parameters alpha, beta, and delta should have unconstrained
      support.",
      `x` = "Found unconstrained distribution {.var {dists}} for parameter{?s}
      {.var {pars}}."))
  }

  priors
}

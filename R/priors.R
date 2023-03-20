#' Construct a Vectorized Priors
#'
#' @param ptype \[character(1L)]\cr Type of the parameter.
#' @param priors \[`data.frame`]\cr Prior definitions.
#' @param channel \[`list()`]\cr Channel-specific variables for Stan sampling
#' @param category \[character(1L)]\cr Category of the categorical response.
#' @noRd
create_vectorized_prior <- function(ptype, priors, channel, category = "") {
  ycat <- ifelse(
    nchar(category) > 0L,
    paste0("_", category),
    ""
  )
  pdef <- priors[priors$type == ptype, ]
  vectorized_prior <- list()
  vectorized_prior[[paste0(ptype, "_prior_distr", ycat)]] <- pdef$prior
  dists <- sub("\\(.*", "", pdef$prior)
  if (nrow(pdef) > 0L && identical(n_unique(dists), 1L)) {
    pars <- strsplit(sub(".*\\((.*)\\).*", "\\1", pdef$prior), ",")
    pars <- do.call("rbind", lapply(pars, as.numeric))
    vectorized_prior[[paste0(ptype, "_prior_distr", ycat)]] <- dists[1L]
    vectorized_prior[[paste0(ptype, "_prior_npars", ycat)]] <- ncol(pars)
    vectorized_prior[[paste0(ptype, "_prior_pars", ycat)]] <- pars
    vectorized_prior[[paste0("vectorized_", ptype, ycat)]] <- TRUE
  } else {
    vectorized_prior[[paste0("vectorized_", ptype, ycat)]] <- FALSE
  }
  vectorized_prior
}
#' Find Vectorizable Priors
#' @param priors \[`data.frame`]\cr Prior definitions.
#' @param y \[`character(1)`]\cr Name of the response variable.
#' @param category \[character(1L)]\cr Category of the categorical response.
#' @noRd
extract_vectorizable_priors <- function(priors, y, category = "") {
  ycat <- ifelse(
    nchar(category) > 0L,
    paste0(".", category),
    ""
  )
  priors <- priors[which(
    names(priors) %in% paste0(
      c("alpha", "beta", "delta", "tau", "sigma_nu"), "_prior_pars", ycat
    )
  )]
  ifelse_(
    length(priors) > 0,
    setNames(priors, paste0(names(priors), "_", y, ycat)),
    NULL
  )
}
#' Construct Common Priors among Channels
#'
#' @inheritParams splines
#' @param priors Custom prior definitions or `NULL` if not specified.
#' @param M Number of random effects.
#' @param shrinkage Does the model contain shrinkage parameter?
#' @param correlated_nu Does the model contain correlated random effects?
#' @param P Number of channels with latent factor.
#' @param correlated_lf Does the model contain correlated latent factors?
#' @noRd
prepare_common_priors <- function(priors, M, shrinkage, P,
                                  correlated_nu, correlated_lf) {
  common_priors <- NULL
  if (shrinkage) {
    common_priors <- ifelse_(
      is.null(priors),
      data.frame(
        parameter = "xi",
        response = "",
        prior = "normal(0, 1)",
        type = "xi",
        category = ""
      ),
      priors[priors$type == "xi", ]
    )
  }
  if (M > 1L && correlated_nu) {
    common_priors <- ifelse_(
      is.null(priors),
      rbind(
        common_priors,
        data.frame(
          parameter = "L_nu",
          response = "",
          prior = "lkj_corr_cholesky(1)",
          type = "L",
          category = ""
        )
      ),
      rbind(
        common_priors,
        priors[priors$parameter == "L_nu", ]
      )
    )
  }
  if (P > 1L && correlated_lf) {
    common_priors <- ifelse_(
      is.null(priors),
      rbind(
        common_priors,
        data.frame(
          parameter = "L_lf",
          response = "",
          prior = "lkj_corr_cholesky(1)",
          type = "L",
          category = ""
        )
      ),
      rbind(
        common_priors,
        priors[priors$parameter == "L_lf", ]
      )
    )
  }
  common_priors
}

#' Create Default Priors
#'
#' @param y \[`character(1)`]\cr Name of the response variable of the channel.
#' @param channel \[`list()`]\cr Channel-specific helper variables.
#' @param mean_gamma Prior mean betas and deltas (at time `fixed + 1`).
#' @param sd_gamma Prior SD betas and deltas (at time `fixed + 1`).
#' @param mean_y Mean of the response variable at time `fixed + 1`.
#' @param sd_y Standard deviation of the response variable at time `fixed + 1`.
#' @param category Name of the category for categorical response variable.
#' @noRd
default_priors <- function(y, channel, mean_gamma, sd_gamma, mean_y, sd_y,
                           category = "") {
  ycat <- ifelse(
    nchar(category) > 0L,
    paste0(".", category),
    ""
  )
  mean_y <- signif(mean_y, 2)
  sd_y <- signif(2 * max(1, sd_y), 2)
  mean_gamma <- signif(mean_gamma, 2)
  sd_gamma <- signif(pmax(sd_gamma, 1), 2)
  priors <- list()
  prior_distributions <- list()
  if (channel$has_random_intercept || channel$has_random) {
    icpt <- ifelse_(
      channel$has_random_intercept,
      c(alpha = sd_y),
      NULL
    )
    s <- c(
      icpt,
      ifelse_(
        channel$has_random,
        sd_gamma[channel$J_random],
        NULL
      )
    )
    ns <- names(s)
    prior_distributions$sigma_nu_prior_distr <- "normal"
    prior_distributions$sigma_nu_prior_npars <- 2L
    prior_distributions$sigma_nu_prior_pars <- cbind(0, s)
    prior_distributions$vectorized_sigma_nu <- TRUE
    priors$sigma_nu <- data.frame(
      parameter = paste0("sigma_nu_", y, ycat, "_", ns),
      response = y,
      prior = paste0("normal(0, ", s, ")"),
      type = "sigma_nu",
      category = category
    )
  } else {
    prior_distributions$vectorized_sigma_nu <- FALSE
  }
  if (channel$has_lfactor) {
    prior_distributions$sigma_lambda_prior_distr <- "normal(0, 1)"
    priors$sigma_lambda <- data.frame(
      parameter = paste0("sigma_lambda_", y, ycat),
      response = y,
      prior = prior_distributions$sigma_lambda_prior_distr,
      type = "sigma_lambda",
      category = category
    )
    if (channel$nonzero_lambda) {
      prior_distributions$tau_psi_prior_distr <- paste0("normal(0, ", sd_y, ")")
      priors$tau_psi <- data.frame(
        parameter = paste0("tau_psi_", y, ycat),
        response = y,
        prior = prior_distributions$tau_psi_prior_distr,
        type = "tau_psi",
        category = category
      )
    }
    prior_distributions$psi_prior_distr <- "normal(0, 1)"
    priors$psi <- data.frame(
      parameter = paste0("psi_", y, ycat),
      response = y,
      prior = prior_distributions$psi_prior_distr,
      type = "psi",
      category = category
    )
  }
  if (channel$has_fixed_intercept || channel$has_varying_intercept) {
    prior_distributions$alpha_prior_distr <-
      paste0("normal(", mean_y, ", ", sd_y, ")")
    priors$alpha <- data.frame(
      parameter = paste0("alpha_", y, ycat),
      response = y,
      prior = prior_distributions$alpha_prior_distr,
      type = "alpha",
      category = category
    )
    if (channel$has_varying_intercept) {
      prior_distributions$tau_alpha_prior_distr <- paste0("normal(0, ", sd_y, ")")
      priors$tau_alpha <- data.frame(
        parameter = paste0("tau_alpha_", y, ycat),
        response = y,
        prior = prior_distributions$tau_alpha_prior_distr,
        type = "tau_alpha",
        category = category
      )
    }
  }
  if (channel$has_fixed) {
    m <- mean_gamma[channel$J_fixed]
    s <- sd_gamma[channel$J_fixed]
    prior_distributions$beta_prior_distr <- "normal"
    prior_distributions$beta_prior_npars <- 2L
    prior_distributions$beta_prior_pars <- unname(cbind(m, s))
    prior_distributions$vectorized_beta <- TRUE
    priors$beta <- data.frame(
      parameter = paste0("beta_", y, ycat, "_", names(s)),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "beta",
      category = category
    )
  } else {
    prior_distributions$vectorized_beta <- FALSE
  }
  if (channel$has_varying) {
    m <- mean_gamma[channel$J_varying]
    s <- sd_gamma[channel$J_varying]
    prior_distributions$delta_prior_distr <- "normal"
    prior_distributions$delta_prior_npars <- 2L
    prior_distributions$delta_prior_pars <- unname(cbind(m, s))
    prior_distributions$vectorized_delta <- TRUE
    priors$delta <- data.frame(
      parameter = paste0("delta_", y, ycat, "_", names(s)),
      response = y,
      prior = paste0("normal(", m, ", ", s, ")"),
      type = "delta",
      category = category
    )
    prior_distributions$tau_prior_distr <- "normal"
    prior_distributions$tau_prior_npars <- 2
    prior_distributions$tau_prior_pars <- cbind(0, s)
    prior_distributions$vectorized_tau <- TRUE
    priors$tau <- data.frame(
      parameter = paste0("tau_", y, ycat, "_", names(s)),
      response = y,
      prior = paste0("normal(0, ", s, ")"),
      type = "tau",
      category = category
    )
  } else {
    prior_distributions$vectorized_delta <- FALSE
    prior_distributions$vectorized_tau <- FALSE
  }
  list(
    prior_distributions = prior_distributions,
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

  dupl <- duplicated(priors$parameter)
  stopifnot_(
    all(!dupl),
    c(
      "Argument {.var priors} contains multiple priors for same parameter.",
      `x` = "{cli::qty(sum(dupl))} Found multiple prior{?s} for parameter{?s}
             {.var {priors$parameter[dupl]}}."
    )
  )
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
    "skew_normal", "logistic", "gumbel", "skew_double_exponential", "std_normal"
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
}
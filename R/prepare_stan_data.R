#' Prepare data for Stan
#'
#' Prepares data for Stan sampling, Stan model code construction and
#' default/user-modifiable prior definitions.
#'
#' @param data \[`data.table`]\cr The data used for model fitting
#' @param dformula \[`dynamiteformula`]\cr The model formula of stochastic
#'   channels
#' @param group_var \[`character(1)`]\cr The grouping variable name
#' @param time_var \[`character(1)`]\cr The time index variable name
#' @param priors TODO
#' @param fixed \[`integer(1)`]\cr Number of fixed time points
#' @srrstats {RE2.3} *Where applicable, Regression Software should enable data to be centred (for example, through converting to zero-mean equivalent values; or to z-scores) or offset (for example, to zero-intercept equivalent values) via additional parameters, with the effects of any such parameters clearly documented and tested.*
#' @srrstats {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstats {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @srrstats {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#' @srrstats {G2.4c} *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
#' @srrstats {G2.4d} *explicit conversion to factor via `as.factor()`*
#' @srrstats {G2.4e} *explicit conversion from factor via `as...()` functions*
#' @srrstats {G2.16} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.*
#' @srrstats {BS2.1} *Bayesian Software should implement pre-processing routines to ensure all input data is dimensionally commensurate, for example by ensuring commensurate lengths of vectors or numbers of rows of tabular inputs.*
#' @srrstats {BS2.1a} *The effects of such routines should be tested.*
#' @srrstats {BS2.2} *Ensure that all appropriate validation and pre-processing of distributional parameters are implemented as distinct pre-processing steps prior to submitting to analytic routines, and especially prior to submitting to multiple parallel computational chains.*
#' @srrstats {BS2.3} *Ensure that lengths of vectors of distributional parameters are checked, with no excess values silently discarded (unless such output is explicitly suppressed, as detailed below).*
#' @srrstats {BS2.4} *Ensure that lengths of vectors of distributional parameters are commensurate with expected model input (see example immediately below)*
#' @srrstats {BS2.5} *Where possible, implement pre-processing checks to validate appropriateness of numeric values submitted for distributional parameters; for example, by ensuring that distributional parameters defining second-order moments such as distributional variance or shape parameters, or any parameters which are logarithmically transformed, are non-negative.*
# TODO prior validation, conversion warnings?
# TODO missing values in prior computations mean_x, sd etc
# TODO X and y starting from fixed + 1, modify T_full
#' @noRd
prepare_stan_data <- function(data, dformula, group_var, time_var,
                              priors = NULL, fixed) {

  resp_names <- get_responses(dformula)
  responses <- as.data.frame(data[, .SD, .SDcols = resp_names])
  # Needs sapply/lapply instead of apply to keep factors as factors
  attr(responses, "resp_class") <- lapply(responses, function(x) {
    cl <- class(x)
    attr(cl, "levels") <- levels(x)
    cl
  })
  specials <- evaluate_specials(dformula, data)
  model_matrix <- full_model.matrix(dformula, data)
  #resp_names <- colnames(responses)
  n_channels <- length(resp_names)
  # A list of variables for stan sampling without grouping by channel
  sampling_vars <- list()
  empty_list <- setNames(vector(mode = "list", length = n_channels), resp_names)
  # A list containing a list for each channel consisting of
  # variables used to construct the stan model code
  model_vars <- empty_list
  # A list for getting current prior definitions
  prior_list <- empty_list

  time <- sort(unique(data[[time_var]]))
  T_full <- length(time)
  group <- NULL
  groups <- !is.null(group_var)
  if (groups) {
    group <- data[[group_var]]
  }

  spline_defs <- attr(dformula, "splines")
  has_splines <- !is.null(spline_defs)
  lb <- numeric(n_channels)
  shrinkage <- FALSE
  if (has_splines) {
    lb <- spline_defs$lb_tau
    shrinkage <- spline_defs$shrinkage
    bs_opts <- spline_defs$bs_opts
    bs_opts$x <- time # TODO change this
    if (is.null(bs_opts$Boundary.knots)) {
      bs_opts$Boundary.knots <- range(bs_opts$x)
    }
    Bs <- t(do.call(splines::bs, args = bs_opts))
    D <- nrow(Bs)
    noncentered <- spline_defs$noncentered
    if (length(noncentered) %in% c(1L, n_channels)) {
      noncentered <- rep(noncentered, length = n_channels)
    } else {
      warning_(
        "Length of the {.var noncentered} argument of {.fun splines} function
        is not equal to 1 or the number of the channels.",
        `i` = "Recycling."
      )
    }
    sampling_vars$D <- D
    sampling_vars$Bs <- Bs
  }
  N <- 1L
  if (groups) {
    N <- length(unique(group))
  }
  K <- ncol(model_matrix)
  X <- aperm(
    array(as.numeric(unlist(split(model_matrix, gl(T_full, 1, N * T_full)))),
      dim = c(N, K, T_full)
    ),
    c(3, 1, 2)
  )
  sd_x <- apply(X, 3, sd, na.rm = TRUE)
  # needed for default priors, 0.5 is pretty arbitrary
  sd_x <- setNames(pmax(0.5, sd_x, na.rm = TRUE),
                   colnames(model_matrix))
  # TODO use fixed + 1
  x_means <- apply(X[1, , , drop = FALSE], 3, mean, na.rm = TRUE)
  # For missing lagged covariates etc
  x_means[is.na(x_means)] <- 0
  X_na <- is.na(X)
  # Placeholder for NAs in Stan
  X[X_na] <- 0
  assigned <- attr(model_matrix, "assign")
  fixed_pars <- attr(model_matrix, "fixed")
  varying_pars <- attr(model_matrix, "varying")
  resp_classes <- attr(responses, "resp_class")
  warn_nosplines <- FALSE
  for (i in seq_len(n_channels)) {
    channel <- list()
    resp <- resp_names[i]
    if (groups) {
      resp_split <- split(responses[, resp], group)
    } else {
      resp_split <- responses[, resp]
    }
    Y <- array(as.numeric(unlist(resp_split)), dim = c(T_full, N))
    Y_na <- is.na(Y)
    # Separate copy of Y for Stan, so that added zeros do not influence channel
    # preparation nor influence other checks related to response variables.
    Y_out <- Y
    # Placeholder for NAs in Stan
    Y_out[Y_na] <- 0
    form_specials <- specials[[i]]
    channel$resp <- resp
    channel$L_fixed <- as.array(match(fixed_pars[[i]], assigned[[i]]))
    channel$L_varying <- as.array(match(varying_pars[[i]], assigned[[i]]))
    channel$J <- as.array(assigned[[i]])
    channel$J_fixed <- as.array(fixed_pars[[i]]) #as.array(assigned[[i]][channel$L_fixed])
    channel$J_varying <- as.array(varying_pars[[i]]) #as.array(assigned[[i]][channel$L_varying])
    channel$K <- length(assigned[[i]])
    channel$K_fixed <- length(fixed_pars[[i]])
    channel$K_varying <- length(varying_pars[[i]])
    obs_idx <- array(0, dim = c(N, T_full))
    obs_len <- integer(T_full)
    for (j in seq_len(T_full)) {
      x_na <- X_na[j, , channel$J, drop = FALSE]
      dim(x_na) <- c(N, channel$K)
      y_na <- Y_na[j, ]
      obs_XY <- which(apply(x_na, 1, function(z) all(!z)) & !y_na)
      obs_XY_len <- length(obs_XY)
      obs_idx[, j] <- c(obs_XY, rep(0, N - obs_XY_len))
      obs_len[j] <- obs_XY_len
    }
    channel$has_missing <- any(obs_len < N)
    if (channel$has_missing) {
      sampling_vars[[paste0("obs_", resp)]] <- obs_idx
      sampling_vars[[paste0("n_obs_", resp)]] <- obs_len
      channel$obs <- glue::glue("obs_{resp}[1:n_obs_{resp}[t], t]")
    } else {
      channel$obs <- ""
    }
    channel$has_fixed_intercept <- dformula[[i]]$has_fixed_intercept
    channel$has_varying_intercept <- dformula[[i]]$has_varying_intercept
    channel$has_random_intercept <- dformula[[i]]$has_random_intercept
    channel$has_fixed <- channel$K_fixed > 0
    channel$has_varying <- channel$K_varying > 0
    channel$lb <- lb[i]
    channel$shrinkage <- shrinkage
    if (channel$has_varying || channel$has_varying_intercept) {
      if (!has_splines) {
        stop_("Model for response variable {.var {resp}}
              contains time-varying definitions
              but splines have not been defined.")
        # TODO switch back to warning after defining default splines?
        #warn_nosplines <- TRUE
      }
      if (warn_nosplines) {
        # channel$has_varying <- FALSE
        # channel$noncentered <- FALSE
      } else {
        channel$noncentered <- noncentered[i]
      }
    } else {
      channel$noncentered <- FALSE
    }
    for (spec in formula_special_funs) {
      if (!is.null(form_specials[[spec]])) {
        if (groups) {
          spec_split <- split(form_specials[[spec]], group)
        } else {
          spec_split <- form_specials[[spec]]
        }
        sampling_vars[[paste0(spec, "_", resp)]] <-
          array(as.numeric(unlist(spec_split)), dim = c(T_full, N))
        channel[[paste0("has_", spec)]] <- TRUE
      } else {
        channel[[paste0("has_", spec)]] <- FALSE
      }
    }
    family <- dformula[[i]]$family
    if (family %in% c("gaussian", "gamma", "exponential") ) {
      sampling_vars[[resp]] <- t(Y_out)
    } else {
      sampling_vars[[resp]] <- Y_out
    }
    prep <- do.call(
      paste0("prepare_channel_", family),
      list(
        y = resp,
        Y = Y,
        channel,
        sd_x = sd_x,
        resp_class = resp_classes[[resp]],
        priors = priors
      )
    )
    prior_list[[resp]] <- prep$priors
    model_vars[[resp]] <- prep$channel
    sampling_vars <- c(sampling_vars, prep$sampling_vars)
  }
  # TODO do we have default splines?
  #if (warn_nosplines) {
  #  warning_("All channels will now default to time-constant ",
  #           "coefficients for all predictors.")
  #}
  sampling_vars$N <- N
  sampling_vars$K <- K
  sampling_vars$X <- X
  sampling_vars$T <- T_full
  sampling_vars$X_m <- as.array(x_means)
  list(
    model_vars = model_vars,
    sampling_vars = sampling_vars,
    priors = prior_list,
    responses = responses,
    u_names = colnames(model_matrix)
  )
}

#' Default channel preparation
#'
#' Computes default channel-specific variables for Stan sampling,
#' Stan model code construction, and prior definitions.
#'
#' @param y \[`character(1)`]\cr Name of the response variable of the channel.
#' @param Y \[`matrix()`]\cr A matrix of values of the response variable.
#' @param channel \[`list()`]\cr Channel-specific helper variables.
#' @param sd_gamma TODO
#' @param resp_class \[`character()`]\cr Class(es) of the response `Y`.
#' @param priors TODO
#'
#' @srrstats {RE1.2} *Regression Software should document expected format (types or classes) for inputting predictor variables, including descriptions of types or classes which are not accepted.*
#' @noRd
prepare_channel_default <- function(y, Y, channel, mean_gamma, sd_gamma,
                                    mean_y, sd_y, resp_class, priors) {
  if (is.null(priors)) {
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
    priors <- dplyr::bind_rows(priors)
  } else {
    # TODO add a warning to documentation that the only the 'prior' column
    # of the priors data.frame should be altered (i.e. there's no checks for names or reordering of rows)
    # Or arrange...
    priors <- priors |> dplyr::filter(.data$response == y)
    for (ptype in c("alpha", "tau_alpha", "sigma_nu")) {
      pdef <- priors |> dplyr::filter(.data$type == ptype)
      if (nrow(pdef) > 0) {
        channel[[paste0(ptype, "_prior_distr")]] <- pdef$prior
      }
    }
    for (ptype in c("beta", "delta", "tau")) {
      pdef <- priors |> dplyr::filter(.data$type == ptype)
      if (nrow(pdef) > 0) {
        dists <- sub("\\(.*", "", pdef$prior)
        vectorized <- length(unique(dists)) == 1
        if (vectorized) {
          pars <- strsplit(sub(".*\\((.*)\\).*", "\\1", pdef$prior), ",")
          pars <- do.call("rbind", lapply(pars, as.numeric))
          channel[[paste0(ptype, "_prior_npars")]] <- ncol(pars)
          channel[[paste0(ptype, "_prior_pars")]] <- pars
          channel[[paste0(ptype, "_prior_distr")]] <- dists[1]
        } else {
          channel[[paste0(ptype, "_prior_distr")]] <- pdef$prior
        }
      }
    }
  }
  channel$write_beta <- channel$has_fixed &&
    length(channel$beta_prior_distr) == 1
  channel$write_delta <- channel$has_varying &&
    length(channel$delta_prior_distr) == 1
  channel$write_tau <- channel$has_varying &&
    length(channel$tau_prior_distr) == 1
  list(channel = channel, priors = priors)
}

#' @describeIn prepare_channel_default Prepare a categorical channel
#' @noRd
prepare_channel_categorical <- function(y, Y, channel, sd_x, resp_class,
                                        priors) {
  if (!("factor" %in% resp_class)) {
    stop_(c(
      "Response variable {.var {y}} is invalid:",
      `x` = "Categorical family supports only <factor> variables."
    ))
  }
  S_y <- length(unique(na.exclude(as.vector(Y))))
  channel$S <- S_y
  if (is.null(priors)) {
    # remove the first level which acts as reference
    resp_levels <- attr(resp_class, "levels")[-1]
    priors <- list()
    sd_gamma <- 2 / sd_x

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
    priors <- dplyr::bind_rows(priors)
  } else {
    # TODO add a warning to documentation that the only the 'prior' column
    # of the priors data.frame should be altered (i.e. there's no checks for names or reordering of rows)
    # Or arrange...
    priors <- priors |> dplyr::filter(.data$response == y)
    for (ptype in c("alpha", "tau_alpha", "beta", "delta", "tau")) {
      pdef <- priors |> dplyr::filter(.data$type == ptype)
      if (nrow(pdef) > 0) {
        dists <- sub("\\(.*", "", pdef$prior)
        vectorized <- length(unique(dists)) == 1
        if (vectorized) {
          pars <- strsplit(sub(".*\\((.*)\\).*", "\\1", pdef$prior), ",")
          pars <- do.call("rbind", lapply(pars, as.numeric))
          channel[[paste0(ptype, "_prior_npars")]] <- ncol(pars)
          channel[[paste0(ptype, "_prior_pars")]] <- pars
          channel[[paste0(ptype, "_prior_distr")]] <- dists[1]
        } else {
          channel[[paste0(ptype, "_prior_distr")]] <- pdef$prior # write separate priors
        }
      }
    }
  }
  channel$write_alpha <-
    (channel$has_fixed_intercept || channel$has_varying_intercept) &&
    length(channel$alpha_prior_distr) == 1
  channel$write_beta <- channel$has_fixed &&
    length(channel$beta_prior_distr) == 1
  channel$write_delta <- channel$has_varying &&
    length(channel$delta_prior_distr) == 1
  channel$write_tau <- channel$has_varying &&
    length(channel$tau_prior_distr) == 1
  list(channel = channel, priors = priors)
}

#' @describeIn prepare_channel_default Prepare a gaussian channel
#' @noRd
prepare_channel_gaussian <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Gaussian", call = rlang::caller_env())
  }
  if (ncol(Y) > 1) {
    sd_y <- mean(apply(Y, 1, sd, na.rm = TRUE))
    mean_y <- mean(Y[1, ], na.rm = TRUE)
  } else {
    sd_y <- sd(Y, na.rm = TRUE)
    mean_y <- Y[1]
  }
  if(is.na(sd_y) || sd_y == 0) {
    sd_y <- 1
  }
  if(is.na(mean_y)) {
    mean_y <- 0
  }
  sd_gamma <- 2 * sd_y / sd_x
  mean_gamma <- rep(0, length(sd_gamma))

  out <- prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma,
                                 mean_y, sd_y, resp_class, priors)
  if (is.null(priors)) {

    out$channel$sigma_prior_distr <- paste0("exponential(", 1 / sd_y, ")")
    out$priors <- dplyr::bind_rows(
      out$priors,
      data.frame(
        parameter = paste0("sigma_", y),
        response = y,
        prior = out$channel$sigma_prior_distr,
        type = "sigma",
        category = ""
      )
    )
  } else {
    pdef <- priors |>
      dplyr::filter(.data$response == y & .data$type == "sigma")
    if (nrow(pdef) == 1) {
      out$channel$sigma_prior_distr <- pdef$prior
    }
  }
  out
}

#' @describeIn prepare_channel_default Prepare a binomial channel
#' @noRd
prepare_channel_binomial <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Binomial", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(y, "Binomial", type = "integers", call = rlang::caller_env())
  }
  # TODO could be adjusted
  sd_y <- 0.5
  mean_y <- 0
  sd_gamma <- 2 / sd_x
  mean_gamma <- rep(0, length(sd_gamma))
  prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma, mean_y, sd_y,
                          resp_class, priors)
}

#' @describeIn prepare_channel_default Prepare a bernoulli channel
#' @noRd
prepare_channel_bernoulli <- function(y, Y, channel, sd_x, resp_class,
                                      priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Bernoulli", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (!all(Y_obs %in% 0:1)) {
    stop_(c(
      "Response variable {.var {y}} is invalid:",
      `x` = "Bernoulli family supports only 0/1 integers."
    ))
  }
  prepare_channel_binomial(y, Y, channel, sd_x, resp_class, priors)
}

#' @describeIn prepare_channel_default Prepare a Poisson channel
#' @noRd
prepare_channel_poisson <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Poisson", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(y, "Poisson", type = "integers", call = rlang::caller_env())
  }
  # TODO could be adjusted
  sd_y <- 1
  if (ncol(Y) > 1) {
    mean_y <- log(mean(Y[1, ], na.rm = TRUE))
  } else {
    mean_y <- log(Y[1])
  }
  if(is.na(mean_y)) {
    mean_y <- 0
  }
  sd_gamma <- 2 / sd_x
  mean_gamma <- rep(0, length(sd_gamma))
  prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma,  mean_y, sd_y,
                          resp_class, priors)
}

#' @describeIn prepare_channel_default Prepare a negative binomial channel
#' @noRd
prepare_channel_negbin <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Negative binomial", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(y, "Negative binomial", type = "integers",
                   call = rlang::caller_env())
  }
  # TODO could be adjusted
  sd_y <- 1
  if (ncol(Y) > 1) {
    mean_y <- log(mean(Y[1, ], na.rm = TRUE))
  } else {
    mean_y <- log(Y[1])
  }
  if(is.na(mean_y)) {
    mean_y <- 0
  }
  sd_gamma <- 2 / sd_x
  mean_gamma <- rep(0, length(sd_gamma))
  out <- prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma,
                                 mean_y, sd_y, resp_class, priors)
  if (is.null(priors)) {
    out$channel$phi_prior_distr <- "exponential(1)"
    out$priors <- dplyr::bind_rows(
      out$priors,
      data.frame(
        parameter = paste0("phi_", y),
        response = y,
        prior = out$channel$phi_prior_distr,
        type = "phi",
        category = ""
      )
    )
  } else {
    pdef <- priors |> dplyr::filter(.data$response == y & .data$type == "phi")
    if (nrow(pdef) == 1) {
      out$channel$phi_prior_distr <- pdef$prior
    }
  }
  out
}

#' @describeIn prepare_channel_default Prepare a gamma channel
#' @noRd
prepare_channel_exponential <- function(y, Y, channel, sd_x, resp_class,
                                        priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Exponential", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs <= 0)) {
    abort_negative(y, "Exponential", type = "values",
                   call = rlang::caller_env())
  }
  # TODO could be adjusted
  sd_y <- 1
  if (ncol(Y) > 1) {
    mean_y <- log(mean(Y[1, ], na.rm = TRUE))
  } else {
    mean_y <- log(Y[1])
  }
  if(is.na(mean_y)) {
    mean_y <- 0
  }
  sd_gamma <- 2 / sd_x
  mean_gamma <- rep(0, length(sd_gamma))
  prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma,
                                 mean_y, sd_y, resp_class, priors)

}

#' @describeIn prepare_channel_default Prepare a gamma channel
#' @noRd
prepare_channel_gamma <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Gamma", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs <= 0)) {
    abort_negative(y, "Gamma", type = "values", call = rlang::caller_env())
  }
  # TODO could be adjusted
  sd_y <- 1
  if (ncol(Y) > 1) {
    mean_y <- log(mean(Y[1, ], na.rm = TRUE))
  } else {
    mean_y <- log(Y[1])
  }
  if(is.na(mean_y)) {
    mean_y <- 0
  }
  sd_gamma <- 2 / sd_x
  mean_gamma <- rep(0, length(sd_gamma))
  out <- prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma,
                                 mean_y, sd_y, resp_class, priors)
  if (is.null(priors)) {
    out$channel$phi_prior_distr <- "exponential(1)"
    out$priors <- dplyr::bind_rows(
      out$priors,
      data.frame(
        parameter = paste0("phi_", y),
        response = y,
        prior = out$channel$phi_prior_distr,
        type = "phi",
        category = ""
      )
    )
  } else {
    pdef <- priors |> dplyr::filter(.data$response == y & .data$type == "phi")
    if (nrow(pdef) == 1) {
      out$channel$phi_prior_distr <- pdef$prior
    }
  }
  out
}

#' Give a warning about nonfinite standard deviation
#'
#' @param y Response varible the warning is related to
#'
#' @noRd
warn_nonfinite <- function(y) {
  warning_(c(
    "Found nonfinite prior standard deviation when using default priors
    for regression coeffients for response {.var {y}}
    indicating constant covariate:",
    `i` = "Switching to N(0, 0.01) prior."
  ))
}

#' Raise error if factor type is not supported by a family
#'
#' @param y Response variable the error ir related to
#' @param family Family as character vector
#' @param call call to be passed to [stop_()]
#'
#' @noRd
abort_factor <- function(y, family, call) {
  stop_(c(
    "Response variable {.var {y}} is invalid:",
    `x` = "{family} family is not supported for {.cls factor} variables."
  ), call = call)
}

#' Raise error if negative values are not supported by a family
#'
#' @param y Response variable the error ir related to
#' @param family Family as character vector
#' @param type Value type of the family
#' @param call call to be passed to [stop_()]
#'
#' @noRd
abort_negative <- function(y, family, type, call) {
  stop_(c(
    "Response variable {.var {y}} is invalid:",
    `x` = "{family} family supports only non-negative {type}."
  ), call = call)
}

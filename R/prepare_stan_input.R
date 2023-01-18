#' Prepare Input for Stan
#'
#' Prepares data for Stan sampling, variables for Stan model code construction
#' and default/user-modifiable prior definitions.
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula of stochastic
#'   channels.
#' @param data \[`data.table`]\cr The data used to fit the model.
#' @param group_var \[`character(1)`]\cr The grouping variable name.
#' @param time_var \[`character(1)`]\cr The time index variable name.
#' @param priors \[`data.frame`]\cr A data frame containing the prior
#'   definitions, or `NULL`, in which case default priors are used.
#' @param fixed \[`integer(1)`]\cr Number of fixed time points.
#' @param verbose \[`logical(1)`]\cr If `TRUE`, outputs warnings.
#'
#' @srrstats {G2.13, G2.14, G2.14a, G2.14b, G2.14c, G2.15}
#'   Missing data is appropriately considered.
#' @srrstats {G2.4, G2.4a, G2.4b, G2.4c, G2.4d, G2.4e}
#'   Data is appropriately converted for Stan.
#' @srrstats {G2.16} Non-finite values are not supported.
#' @srrstats {BS2.1, BS2.1a} Proper dimensionality is ensured for Stan.
#' @srrstats {BS2.2, BS2.3, BS2.4, BS2.5}
#'   Distributional parameters are checked.
#' @noRd
prepare_stan_input <- function(dformula, data, group_var, time_var,
                               priors = NULL, fixed, verbose) {
  resp_names <- get_responses(dformula)
  missing_resp <- !(resp_names %in% names(data))
  stopifnot_(
    all(!missing_resp),
    "Can't find variable{?s} {.var {resp_names[missing_resp]}} in {.arg data}."
  )
  responses <- as.data.frame(data[, .SD, .SDcols = resp_names])
  # Needs lapply instead of apply to keep factors as factors
  attr(responses, "resp_class") <- lapply(responses, function(x) {
    cl <- class(x)
    attr(cl, "levels") <- levels(x)
    cl
  })
  specials <- evaluate_specials(dformula, data)
  model_matrix <- full_model.matrix(dformula, data, verbose)
  n_channels <- length(resp_names)
  # A list of variables for Stan sampling without grouping by channel
  sampling_vars <- list()
  empty_list <- setNames(vector(mode = "list", length = n_channels), resp_names)
  # A list containing a list for each channel consisting of
  # variables used to construct the Stan model code
  model_vars <- empty_list
  # A list for getting current prior definitions
  prior_list <- empty_list
  time <- sort(unique(data[[time_var]]))
  T_full <- length(time)
  T_idx <- seq.int(fixed + 1L, T_full)
  has_groups <- !is.null(group_var)
  group <- data[[group_var]]
  spline_defs <- attr(dformula, "splines")
  random_defs <- attr(dformula, "random_spec")
  lfactor_defs <- attr(dformula, "lfactor")
  has_splines <- spline_defs$has_splines
  sampling_vars$D <- spline_defs$D
  sampling_vars$Bs <- spline_defs$Bs
  has_lfactor <- lfactor_defs$has_lfactor
  stopifnot_(
    !has_lfactor || has_splines,
    "Model contains time-varying latent factor(s) but splines have not been
    defined."
  )
  N <- n_unique(group)
  K <- ncol(model_matrix)
  X <- aperm(
    array(as.numeric(unlist(split(model_matrix, gl(T_full, 1, N * T_full)))),
      dim = c(N, K, T_full)
    ),
    c(3L, 1L, 2L)
  )[T_idx, , , drop = FALSE]
  sd_x <- apply(X[1L, , , drop = FALSE], 3L, sd, na.rm = TRUE)
  # needed for default priors, 0.5 is pretty arbitrary
  sd_x <- setNames(pmax(0.5, sd_x, na.rm = TRUE), colnames(model_matrix))
  x_means <- apply(X[1L, , , drop = FALSE], 3L, mean, na.rm = TRUE)
  # For totally missing covariates
  x_means[is.na(x_means)] <- 0.0
  X_na <- is.na(X)
  # Placeholder for NAs in Stan
  X[X_na] <- 0.0
  assigned <- attr(model_matrix, "assign")
  fixed_pars <- attr(model_matrix, "fixed")
  varying_pars <- attr(model_matrix, "varying")
  random_pars <- attr(model_matrix, "random")
  resp_classes <- attr(responses, "resp_class")
  random_defs$M <- sum(lengths(random_pars)) +
    sum(unlist(lapply(dformula, "[[", "has_random_intercept")))
  lfactor_defs$P <- length(lfactor_defs$responses)
  for (i in seq_len(n_channels)) {
    channel <- list()
    resp <- resp_names[i]
    resp_split <- split(responses[, resp], group)
    Y <- array(as.numeric(unlist(resp_split)), dim = c(T_full, N))
    Y <- Y[T_idx, , drop = FALSE]
    Y_na <- is.na(Y)
    # Separate copy of Y for Stan, so that added zeros do not influence channel
    # preparation nor influence other checks related to response variables.
    Y_out <- Y
    # Placeholder for NAs in Stan
    Y_out[Y_na] <- 0.0
    form_specials <- specials[[i]]
    channel$resp <- resp
    channel$L_fixed <- as.array(match(fixed_pars[[i]], assigned[[i]]))
    channel$L_varying <- as.array(match(varying_pars[[i]], assigned[[i]]))
    channel$J <- as.array(assigned[[i]])
    channel$J_fixed <- as.array(fixed_pars[[i]])
    channel$J_varying <- as.array(varying_pars[[i]])
    channel$J_random <- as.array(random_pars[[i]])
    channel$K <- length(assigned[[i]])
    channel$K_fixed <- length(fixed_pars[[i]])
    channel$K_varying <- length(varying_pars[[i]])
    # note! Random intercept is counted to K_random but not to J_random...
    channel$K_random <- length(random_pars[[i]]) +
      dformula[[i]]$has_random_intercept
    obs_idx <- array(0L, dim = c(N, T_full - fixed))
    obs_len <- integer(T_full - fixed)
    for (j in seq_len(T_full - fixed)) {
      x_na <- X_na[j, , channel$J, drop = FALSE]
      dim(x_na) <- c(N, channel$K)
      y_na <- Y_na[j, ]
      obs_XY <- which(apply(x_na, 1L, function(z) all(!z)) & !y_na)
      obs_XY_len <- length(obs_XY)
      obs_idx[, j] <- c(obs_XY, rep(0L, N - obs_XY_len))
      obs_len[j] <- obs_XY_len
    }
    channel$has_missing <- any(obs_len < N)
    sampling_vars[[paste0("obs_", resp)]] <- obs_idx
    sampling_vars[[paste0("n_obs_", resp)]] <- obs_len
    # obs selects complete cases if there are missing observations
    channel$obs <- ifelse_(
      channel$has_missing,
      glue::glue("obs_{resp}[1:n_obs_{resp}[t], t]"),
      ""
    )
    channel$has_fixed_intercept <- dformula[[i]]$has_fixed_intercept
    channel$has_varying_intercept <- dformula[[i]]$has_varying_intercept
    channel$has_random_intercept <-  dformula[[i]]$has_random_intercept
    channel$has_fixed <- channel$K_fixed > 0L
    channel$has_varying <- channel$K_varying > 0L
    # note! Random intercept is counted to K_random above, while has_random is
    # for checking non-intercept terms....
    channel$has_random <- channel$K_random > 1L
    channel$lb <- spline_defs$lb[i]
    channel$shrinkage <- spline_defs$shrinkage
    channel$noncentered <- spline_defs$noncentered[i]
    channel$has_lfactor <- resp %in% lfactor_defs$responses
    channel$noncentered_psi <- lfactor_defs$noncentered_psi
    channel$noncentered_lambda <- lfactor_defs$noncentered_lambda[i]
    channel$nonzero_lambda <- lfactor_defs$nonzero_lambda[i]
    stopifnot_(
      has_splines || !(channel$has_varying || channel$has_varying_intercept),
      "Model for response variable {.var {resp}} contains time-varying
       definitions but splines have not been defined."
    )
    # evaluate specials such as offset and trials
    for (spec in formula_special_funs) {
      if (!is.null(form_specials[[spec]])) {
        spec_split <- split(form_specials[[spec]], group)
        spec_array <- array(as.numeric(unlist(spec_split)), dim = c(T_full, N))
        sampling_vars[[paste0(spec, "_", resp)]] <-
          spec_array[seq.int(fixed + 1L, T_full), , drop = FALSE]
        channel[[paste0("has_", spec)]] <- TRUE
      } else {
        channel[[paste0("has_", spec)]] <- FALSE
      }
    }
    family <- dformula[[i]]$family

    stopifnot_(!channel$has_random || family != "categorical",
      "Random effects are not (yet) supported for categorical responses.")
    sampling_vars[[paste0("y_", resp)]] <- ifelse_(
      family %in% c("gaussian", "gamma", "exponential", "beta"),
      t(Y_out),
      Y_out
    )
    prep <- do.call(
      paste0("prepare_channel_", family),
      list(
        y = resp,
        Y = Y,
        channel = channel,
        sd_x = sd_x,
        resp_class = resp_classes[[resp]],
        priors = priors
      )
    )
    prior_list[[resp]] <- prep$priors
    model_vars[[resp]] <- prep$channel
    sampling_vars <- c(sampling_vars, prep$sampling_vars)
  }
  sampling_vars$N <- N
  sampling_vars$K <- K
  sampling_vars$X <- X
  sampling_vars$M <- random_defs$M
  sampling_vars$P <- lfactor_defs$P
  # avoid goodpractice warning, T is a Stan variable, not an R variable
  sampling_vars[["T"]] <- T_full - fixed
  sampling_vars$X_m <- as.array(x_means)
  prior_list[["common_priors"]] <- prepare_common_priors(
    priors = priors,
    M = sampling_vars$M,
    shrinkage = spline_defs$shrinkage,
    correlated_nu = random_defs$correlated,
    P = sampling_vars$P,
    correlated_lf = lfactor_defs$correlated
  )
  # for stanblocks
  attr(model_vars, "common_priors") <- prior_list[["common_priors"]]
  attr(model_vars, "spline_defs") <- spline_defs
  attr(model_vars, "random_defs") <- random_defs
  attr(model_vars, "lfactor_defs") <- lfactor_defs
  list(
    model_vars = model_vars,
    sampling_vars = sampling_vars,
    priors = prior_list,
    responses = responses,
    u_names = colnames(model_matrix),
    fixed = fixed
  )
}

#' Construct a Prior Definition for a Regression Parameters
#'
#' @param ptype \[character(1L)]\cr Type of the parameter.
#' @param priors \[`data.frame`]\cr Prior definitions.
#' @param channel \[`list()`]\cr Channel-specific variables for Stan sampling
#' @noRd
prepare_prior <- function(ptype, priors, channel) {
  pdef <- priors[priors$type == ptype, ]
  channel[[paste0(ptype, "_prior_distr")]] <- pdef$prior
  dists <- sub("\\(.*", "", pdef$prior)
  if (nrow(pdef) > 0L && identical(n_unique(dists), 1L)) {
    pars <- strsplit(sub(".*\\((.*)\\).*", "\\1", pdef$prior), ",")
    pars <- do.call("rbind", lapply(pars, as.numeric))
    channel[[paste0(ptype, "_prior_distr")]] <- dists[1L]
    channel[[paste0(ptype, "_prior_npars")]] <- ncol(pars)
    channel[[paste0(ptype, "_prior_pars")]] <- pars
  }
  channel
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
prepare_common_priors <- function(priors, M, shrinkage, P, correlated_nu,
  correlated_lf) {
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

#' Default channel preparation
#'
#' Computes default channel-specific variables for Stan sampling,
#' Stan model code construction, and prior definitions.
#'
#' @param y \[`character(1)`]\cr Name of the response variable of the channel.
#' @param Y \[`matrix()`]\cr A matrix of values of the response variable.
#' @param channel \[`list()`]\cr Channel-specific helper variables.
#' @param mean_gamma \[`numeric(1)`]\cr Prior mean betas and deltas
#'   (at time `fixed + 1`).
#' @param sd_gamma \[`numeric(1)`]\cr  Prior SD betas and deltas
#'   (at time `fixed + 1`).
#' @param mean_y \[`numeric(1)`]\cr  Mean of the response variable
#'   (at time `fixed + 1`).
#' @param sd_y \[`numeric(1)`]\cr  SD of the response variable
#'   (at time `fixed + 1`).
#' @param resp_class \[`character()`]\cr Class(es) of the response `Y`.
#' @param priors \[`data.frame`]\cr Prior definitions, or `NULL`, in which case
#'   the default priors are used.
#'
#' @srrstats {RE1.2} Checks for expected types and classes along with other
#'   `prepare_channel_*` functions.
#' @noRd
prepare_channel_default <- function(y, Y, channel, mean_gamma, sd_gamma,
                                    mean_y, sd_y, resp_class, priors) {
  if (is.null(priors)) {
    out <- default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)
    channel <- out$channel
    priors <- out$priors
  } else {
    priors <- priors[priors$response == y, ]
    types <- priors$type
    for (ptype in intersect(types,
      c("alpha", "tau_alpha", "sigma_nu", "sigma_lambda", "tau_psi"))) {
      pdef <- priors[priors$type == ptype, ]
      channel[[paste0(ptype, "_prior_distr")]] <- pdef$prior
    }
    for (ptype in intersect(types, c("beta", "delta", "tau"))) {
      channel <- prepare_prior(ptype, priors, channel)
    }
  }
  channel$write_beta <- channel$has_fixed &&
    identical(length(channel$beta_prior_distr), 1L)
  channel$write_delta <- channel$has_varying &&
    identical(length(channel$delta_prior_distr), 1L)
  channel$write_tau <- channel$has_varying &&
    identical(length(channel$tau_prior_distr), 1L)
  list(channel = channel, priors = priors)
}

#' @describeIn prepare_channel_default Prepare a Categorical Channel
#' @noRd
prepare_channel_categorical <- function(y, Y, channel, sd_x, resp_class,
                                        priors) {
  stopifnot_(
    "factor" %in% resp_class,
    c(
      "Response variable {.var {y}} is invalid:",
      `x` = "Categorical family supports only {.cls factor} variables."
    )
  )
  S_y <- length(attr(resp_class, "levels"))
  channel$S <- S_y
  if (is.null(priors)) {
    out <- default_priors_categorical(y, channel, sd_x, resp_class)
    channel <- out$channel
    priors <- out$priors
  } else {
    priors <- priors[priors$response == y, ]
    types <- priors$type
    for (ptype in intersect(types, c("alpha", "beta", "delta", "tau"))) {
      channel <- prepare_prior(ptype, priors, channel)
    }
    if ("tau_alpha" %in% types) {
      pdef <- priors[priors$type == "tau_alpha", ]
      channel$tau_alpha_prior_distr <- pdef$prior
    }
    priors <- check_priors(
      priors, default_priors_categorical(y, channel, sd_x, resp_class)$priors
    )
  }
  channel$write_alpha <-
    (channel$has_fixed_intercept || channel$has_varying_intercept) &&
      identical(length(channel$alpha_prior_distr), 1L)
  channel$write_beta <- channel$has_fixed &&
    identical(length(channel$beta_prior_distr), 1L)
  channel$write_delta <- channel$has_varying &&
    identical(length(channel$delta_prior_distr), 1L)
  channel$write_tau <- channel$has_varying &&
    identical(length(channel$tau_prior_distr), 1L)
  list(channel = channel, priors = priors)
}

#' @describeIn prepare_channel_default Prepare a Gaussian Channel
#' @noRd
prepare_channel_gaussian <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Gaussian", call = rlang::caller_env())
  }
  if (ncol(Y) > 1L) {
    sd_y <- mean(apply(Y, 1L, sd, na.rm = TRUE))
    mean_y <- mean(Y[1L, ], na.rm = TRUE)
  } else {
    sd_y <- sd(Y, na.rm = TRUE)
    mean_y <- Y[1L]
  }
  if (is.na(sd_y) || identical(sd_y, 0.0)) {
    sd_y <- 1.0
  }
  if (is.na(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 * sd_y / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    resp_class,
    priors
  )
  sigma_prior <- data.frame(
    parameter = paste0("sigma_", y),
    response = y,
    prior = paste0("exponential(", signif(1.0 / sd_y, 2L), ")"),
    type = "sigma",
    category = ""
  )
  if (is.null(priors)) {
    out$channel$sigma_prior_distr <- sigma_prior$prior
    out$priors <- rbind(out$priors, sigma_prior)
  } else {
    priors <- priors[priors$response == y, ]
    pdef <- priors[priors$type == "sigma", ]
    if (identical(nrow(pdef), 1L)) {
      out$channel$sigma_prior_distr <- pdef$prior
    }
    defaults <- rbind(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      sigma_prior
    )
    out$priors <- check_priors(priors, defaults)
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Binomial Channel
#' @noRd
prepare_channel_binomial <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Binomial", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0.0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(y, "Binomial",
      type = "integers",
      call = rlang::caller_env()
    )
  }
  sd_y <- 0.5
  mean_y <- 0.0
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    resp_class,
    priors
  )
  if (is.null(priors)) {
    out$priors <- check_priors(
      out$priors,
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors
    )
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Bernoulli Channel
#' @noRd
prepare_channel_bernoulli <- function(y, Y, channel, sd_x, resp_class,
                                      priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Bernoulli", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (!all(Y_obs %in% c(0L, 1L))) {
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
  if (any(Y_obs < 0.0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(y, "Poisson", type = "integers", call = rlang::caller_env())
  }
  sd_y <- 1.0
  if (ncol(Y) > 1L) {
    mean_y <- log(mean(Y[1L, ], na.rm = TRUE))
  } else {
    mean_y <- log(Y[1L])
  }
  if (is.na(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    resp_class,
    priors
  )
  if (is.null(priors)) {
    out$priors <- check_priors(
      out$priors,
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors
    )
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Negative Binomial Channel
#' @noRd
prepare_channel_negbin <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Negative binomial", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0.0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(y, "Negative binomial",
      type = "integers",
      call = rlang::caller_env()
    )
  }
  sd_y <- 1.0
  if (ncol(Y) > 1L) {
    mean_y <- log(mean(Y[1L, ], na.rm = TRUE))
  } else {
    mean_y <- log(Y[1L])
  }
  if (is.na(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    resp_class,
    priors
  )
  phi_prior <- data.frame(
    parameter = paste0("phi_", y),
    response = y,
    prior = "exponential(1)",
    type = "phi",
    category = ""
  )

  if (is.null(priors)) {
    out$channel$phi_prior_distr <- phi_prior$prior
    out$priors <- rbind(out$priors, phi_prior)
  } else {
    priors <- priors[priors$response == y, ]
    pdef <- priors[priors$type == "phi", ]
    if (identical(nrow(pdef), 1L)) {
      out$channel$phi_prior_distr <- pdef$prior
    }
    defaults <- rbind(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      phi_prior
    )
    out$priors <- check_priors(priors, defaults)
  }
  out
}

#' @describeIn prepare_channel_default Prepare an Exponential Channel
#' @noRd
prepare_channel_exponential <- function(y, Y, channel, sd_x, resp_class,
                                        priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Exponential", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs <= 0.0)) {
    abort_negative(y, "Exponential",
      type = "values",
      call = rlang::caller_env()
    )
  }
  sd_y <- 1.0
  if (ncol(Y) > 1L) {
    mean_y <- log(mean(Y[1L, ], na.rm = TRUE))
  } else {
    mean_y <- log(Y[1L])
  }
  if (is.na(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    resp_class,
    priors
  )
  if (is.null(priors)) {
    out$priors <- check_priors(
      out$priors,
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors
    )
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Gamma channel
#' @noRd
prepare_channel_gamma <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Gamma", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs <= 0.0)) {
    abort_negative(y, "Gamma", type = "values", call = rlang::caller_env())
  }
  sd_y <- 1.0
  if (ncol(Y) > 1L) {
    mean_y <- log(mean(Y[1L, ], na.rm = TRUE))
  } else {
    mean_y <- log(Y[1L])
  }
  if (is.na(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    resp_class,
    priors
  )
  phi_prior <- data.frame(
    parameter = paste0("phi_", y),
    response = y,
    prior = "exponential(1)",
    type = "phi",
    category = ""
  )

  if (is.null(priors)) {
    out$channel$phi_prior_distr <- phi_prior$prior
    out$priors <- rbind(out$priors, phi_prior)
  } else {
    priors <- priors[priors$response == y, ]
    pdef <- priors[priors$type == "phi", ]
    if (identical(nrow(pdef), 1L)) {
      out$channel$phi_prior_distr <- pdef$prior
    }
    defaults <- rbind(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      phi_prior
    )
    out$priors <- check_priors(priors, defaults)
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Beta Channel
#' @noRd
prepare_channel_beta <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Beta", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs <= 0.0) || any(Y_obs >= 1.0)) {
    abort_nonunit(y, "Beta", type = "values", call = rlang::caller_env())
  }
  sd_y <- 1.0
  if (ncol(Y) > 1L) {
    mean_y <- log(mean(Y[1L, ], na.rm = TRUE))
  } else {
    mean_y <- log(Y[1L])
  }
  if (is.na(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    resp_class,
    priors
  )
  phi_prior <- data.frame(
    parameter = paste0("phi_", y),
    response = y,
    prior = "exponential(1)",
    type = "phi",
    category = ""
  )

  if (is.null(priors)) {
    out$channel$phi_prior_distr <- phi_prior$prior
    out$priors <- rbind(out$priors, phi_prior)
  } else {
    priors <- priors[priors$response == y, ]
    pdef <- priors[priors$type == "phi", ]
    if (identical(nrow(pdef), 1L)) {
      out$channel$phi_prior_distr <- pdef$prior
    }
    defaults <- rbind(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      phi_prior
    )
    out$priors <- check_priors(priors, defaults)
  }
  out
}
#' Raise an error if factor type is not supported by a family
#'
#' @param y \[`character(1)`]\cr Response variable the error is related to.
#' @param family \[`character(1)`]\cr Family as a character string.
#' @param call \[`call`]\cr Call to be passed to [stop_()].
#' @noRd
abort_factor <- function(y, family, call) {
  stop_(
    c(
      "Response variable {.var {y}} is invalid:",
      `x` = "{family} family is not supported for {.cls factor} variables."
    ),
    call = call
  )
}

#' Raise an Error If Negative Values Are Not Supported by a Family
#'
#' @inheritParams abort_factor
#' @param type \[`character(1)`]\cr Value type of the family.
#' @noRd
abort_negative <- function(y, family, type, call) {
  stop_(
    c(
      "Response variable {.var {y}} is invalid:",
      `x` = "{family} family supports only non-negative {type}."
    ),
    call = call
  )
}

#' Raise an Error If Values Are Outside of the Unit Interval
#'
#' @inheritParams abort_negative
#' @noRd
abort_nonunit <- function(y, family, type, call) {
  stop_(
    c(
      "Response variable {.var {y}} is invalid:",
      `x` = "{family} family supports only {type} on the open interval (0, 1)."
    ),
    call = call
  )
}

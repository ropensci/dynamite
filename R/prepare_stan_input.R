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
  resp <- get_responses(dformula)
  resp_names <- get_names(dformula)
  resp_missing <- !(resp %in% names(data))
  stopifnot_(
    all(!resp_missing),
    "Can't find variable{?s} {.var {resp[resp_missing]}} in {.arg data}."
  )
  specials <- lapply(dformula, evaluate_specials, data = data)
  model_matrix <- full_model.matrix(dformula, data, group_var, fixed, verbose)
  cg <- attr(dformula, "channel_groups")
  n_cg <- n_unique(cg)
  n_channels <- length(resp_names)
  responses <- list()
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    y_cg <- paste(resp_names[cg_idx], collapse = "_")
    responses[[y_cg]] <- resp_names[cg_idx]
  }
  attr(responses, "resp_class") <- lapply(
    data[, .SD, .SDcols = resp],
    function(x) {
      cl <- class(x)
      attr(cl, "levels") <- levels(x)
      cl
    }
  )
  # A list of variables for Stan sampling without grouping by channel
  sampling_vars <- list()
  # A list of variables for Stan model construction
  model_vars <- list()
  empty_list <- stats::setNames(
    vector(mode = "list", length = n_channels),
    resp_names
  )
  # A list containing a list for each channel consisting of
  # variables used to construct the Stan model code
  channel_vars <- empty_list
  # A list of channel group specific variables
  channel_group_vars <- list()
  # A list for getting current prior definitions
  prior_list <- empty_list
  time <- sort(unique(data[[time_var]]))
  T_full <- length(time)
  T_idx <- seq.int(fixed + 1L, T_full)
  has_groups <- !is.null(group_var)
  group <- data[[group_var]]
  spline_def <- attr(dformula, "splines")
  random_def <- attr(dformula, "random_spec")
  lfactor_def <- attr(dformula, "lfactor")
  has_splines <- spline_def$has_splines
  sampling_vars$D <- spline_def$D
  sampling_vars$Bs <- spline_def$Bs
  has_lfactor <- lfactor_def$has_lfactor
  stopifnot_(
    !has_lfactor || has_splines,
    "Model contains time-varying latent factor(s) but splines have not been
    defined."
  )
  N <- n_unique(group)
  K <- ncol(model_matrix)
  X <- model_matrix[, ]
  dim(X) <- c(T_full - fixed, N, K)
  x_tmp <- X[1L, , , drop = FALSE]
  sd_x <- pmax(
    stats::setNames(apply(X, 3L, sd, na.rm = TRUE), colnames(model_matrix)),
    1.0
  )
  x_means <- colMeans(x_tmp, dims = 2L, na.rm = TRUE)
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
  for (i in seq_len(n_channels)) {
    y <- resp[i]
    y_name <- resp_names[i]
    Y <- as.numeric(data[[y]])
    dim(Y) <- c(T_full, N)
    Y <- Y[T_idx, , drop = FALSE]
    tmp <- initialize_univariate_channel(
      dformula = dformula[[i]],
      specials = specials[[i]],
      fixed_pars = fixed_pars[[i]],
      varying_pars = varying_pars[[i]],
      random_pars = random_pars[[i]],
      Y = Y,
      y = y,
      y_name = y_name,
      group = group,
      fixed = fixed,
      T_full = T_full,
      N = N,
      X_na = X_na,
      lb = spline_def$lb[i],
      shrinkage = spline_def$shrinkage,
      noncentered = spline_def$noncentered[i],
      has_splines = has_splines,
      has_lfactor = y %in% lfactor_def$responses,
      noncentered_psi = lfactor_def$noncentered_psi,
      flip_sign = lfactor_def$flip_sign,
      nonzero_lambda = lfactor_def$nonzero_lambda[i]
    )
    if (has_univariate(dformula[[i]]$family)) {
      prep <- do.call(
        paste0("prepare_channel_", get_univariate(dformula[[i]]$family)),
        list(
          y = y_name,
          Y = Y,
          channel = tmp$channel,
          sampling = tmp$sampling,
          sd_x = sd_x,
          resp_class = resp_classes[[y]],
          priors = priors
        )
      )
      prior_list[[y_name]] <- prep$priors
      channel_vars[[y_name]] <- prep$channel
      vectorizable_priors <- ifelse_(
        is_categorical(dformula[[i]]$family),
        unlist(
          lapply(
            prep$channel$categories[-1L],
            function(s) {
              extract_vectorizable_priors(
                prep$channel$prior_distr[[s]],
                paste0(y_name, "_", s)
              )
            }
          ),
          recursive = FALSE
        ),
        extract_vectorizable_priors(
          prep$channel$prior_distr,
          y_name
        )
      )
      model_vars$Ks <- c(model_vars$Ks, prep$channel$Ks)
      sampling_vars <- c(
        sampling_vars,
        prep$sampling,
        vectorizable_priors
      )
    } else {
      channel_vars[[y_name]] <- tmp$channel
      sampling_vars <- c(sampling_vars, tmp$sampling)
    }
  }
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    j <- cg_idx[1L]
    family <- dformula[[j]]$family
    channel <- list(family = family)
    y <- resp[cg_idx]
    y_name <- resp_names[cg_idx]
    y_cg <- paste(y_name, collapse = "_")
    if (is_multivariate(family)) {
      tmp <- initialize_multivariate_channel(
        y = y,
        y_cg = y_cg,
        y_name = y_name,
        cg_idx = cg_idx,
        channel = channel,
        family = family,
        dims = c(T_full - fixed, N),
        channel_vars = channel_vars,
        sampling_vars = sampling_vars
      )
      prep <- do.call(
        paste0("prepare_channel_", family$name),
        list(
          y = y_name,
          y_cg = y_cg,
          Y = tmp$sampling[[paste0("y_", y_cg)]],
          channel = tmp$channel,
          sampling = tmp$sampling,
          sd_x = sd_x,
          resp_class = resp_classes[y],
          priors = priors
        )
      )
      vectorizable_priors <- ifelse_(
        is_multinomial(dformula[[i]]$family),
        unlist(
          lapply(
            y_name[-1L],
            function(s) {
              extract_vectorizable_priors(
                prep$channel$prior_distr[[s]],
                s
              )
            }
          ),
          recursive = FALSE
        ),
        extract_vectorizable_priors(
          prep$channel$prior_distr,
          y_name
        )
      )
      model_vars$Ks <- c(model_vars$Ks, prep$channel$Ks)
      sampling_vars[paste0("y_", y_name)] <- NULL
      sampling_vars <- c(
        sampling_vars,
        prep$sampling,
        vectorizable_priors
      )
      channel <- prep$channel
      prior_list[[y_cg]] <- prep$priors
      prior_list$multivariate <- rbind(
        prior_list$multivariate,
        prep$mvpriors
      )
    }
    channel_group_vars[[y_cg]] <- channel
  }
  sampling_vars$N <- N
  sampling_vars$K <- K
  sampling_vars$X <- X
  sampling_vars$M <- random_def$M
  sampling_vars$P <- lfactor_def$P
  # avoid goodpractice warning, T is a Stan variable, not an R variable
  sampling_vars[["T"]] <- T_full - fixed
  sampling_vars$X_m <- as.array(x_means)
  prior_list$common_priors <- prepare_common_priors(
    priors = priors,
    M = sampling_vars$M,
    shrinkage = spline_def$shrinkage,
    correlated_nu = random_def$correlated,
    P = sampling_vars$P,
    correlated_lf = lfactor_def$correlated
  )
  # for stanblocks
  model_vars[["T"]] <- sampling_vars[["T"]]
  model_vars$M <- sampling_vars$M
  model_vars$P <- sampling_vars$P
  model_vars$D <- sampling_vars$D
  model_vars$K <- K
  model_vars$common_priors = prior_list$common_priors
  model_vars$spline_def = spline_def
  model_vars$random_def = random_def
  model_vars$lfactor_def = lfactor_def
  list(
    channel_vars = channel_vars,
    channel_group_vars = channel_group_vars,
    model_vars = model_vars,
    sampling_vars = sampling_vars,
    priors = prior_list,
    responses = responses,
    u_names = colnames(model_matrix),
    fixed = fixed
  )
}

initialize_univariate_channel <- function(dformula, specials, fixed_pars,
                                          varying_pars, random_pars,
                                          Y, y, y_name, group, fixed,
                                          T_full, N, X_na, lb,
                                          shrinkage, noncentered,
                                          has_lfactor, has_splines,
                                          noncentered_psi, flip_sign,
                                          nonzero_lambda) {
  channel <- list()
  Y_na <- is.na(Y)
  # Separate copy of Y for Stan, so that added zeros do not influence channel
  # preparation nor influence other checks related to response variables.
  Y_out <- Y
  # Placeholder for NAs in Stan
  Y_out[Y_na] <- 0.0
  channel$y <- y_name
  channel$family <- dformula$family
  indices <- list(
    K_fixed = length(fixed_pars),
    K_varying = length(varying_pars),
    # Note! Random intercept is counted to K_random but not to J_random...
    K_random = length(random_pars) + dformula$has_random_intercept,
    K = length(fixed_pars) + length(varying_pars),
    J_fixed = as.array(fixed_pars),
    J_varying = as.array(varying_pars),
    J = as.array(c(fixed_pars, varying_pars)),
    J_random = as.array(random_pars),
    L_fixed = as.array(seq_along(fixed_pars)),
    L_varying = as.array(length(fixed_pars) + seq_along(varying_pars))
  )
  channel <- c(channel, indices)
  sampling <- stats::setNames(indices, paste0(names(indices), "_", y_name))
  # evaluate specials such as offset and trials
  spec_na <- matrix(FALSE, nrow = T_full - fixed, ncol = N)
  for (spec in formula_special_funs) {
    if (!is.null(specials[[spec]])) {
      spec_idx <- seq.int(fixed + 1L, T_full)
      spec_array <- as.numeric(specials[[spec]])
      dim(spec_array) <- c(T_full, N)
      spec_na <- spec_na | is.na(spec_array[spec_idx, , drop = FALSE])
      spec_name <- paste0(spec, "_", y_name)
      sampling[[spec_name]] <- ifelse_(
        identical(spec, "offset"),
        t(spec_array[spec_idx, , drop = FALSE]),
        spec_array[spec_idx, , drop = FALSE]
      )
      channel[[paste0("has_", spec)]] <- TRUE
    } else {
      channel[[paste0("has_", spec)]] <- FALSE
    }
  }
  obs_idx <- array(0L, dim = c(N, T_full - fixed))
  obs_len <- integer(T_full - fixed)
  for (j in seq_len(T_full - fixed)) {
    x_na <- X_na[j, , channel$J, drop = FALSE]
    dim(x_na) <- c(N, channel$K)
    y_na <- Y_na[j, ]
    obs_XY <- which(
      apply(x_na, 1L, function(z) all(!z)) & !y_na & all(!spec_na[j, ])
    )
    obs_XY_len <- length(obs_XY)
    obs_idx[, j] <- c(obs_XY, rep(0L, N - obs_XY_len))
    obs_len[j] <- obs_XY_len
  }
  channel$has_missing <- any(obs_len < N)
  channel$has_fully_missing <- any(obs_len == 0L)
  sampling[[paste0("obs_", y_name)]] <- obs_idx
  sampling[[paste0("n_obs_", y_name)]] <- obs_len
  sampling[[paste0("t_obs_", y_name)]] <- as.array(which(obs_len > 0L))
  sampling[[paste0("T_obs_", y_name)]] <- length(which(obs_len > 0L))
  # obs selects complete cases if there are missing observations
  channel$obs <- ifelse_(
    channel$has_missing,
    glue::glue("obs_{y_name}[1:n_obs_{y_name}[t], t]"),
    ""
  )
  channel$has_fixed_intercept <- dformula$has_fixed_intercept
  channel$has_varying_intercept <- dformula$has_varying_intercept
  channel$has_random_intercept <- dformula$has_random_intercept
  channel$has_fixed <- channel$K_fixed > 0L
  channel$has_varying <- channel$K_varying > 0L
  # note! Random intercept is counted to K_random above, while has_random is
  # for checking non-intercept terms....
  channel$has_random <- channel$K_random > channel$has_random_intercept
  channel$lb <- lb
  channel$shrinkage <- shrinkage
  channel$noncentered <- noncentered
  channel$has_lfactor <- has_lfactor
  channel$noncentered_psi <- noncentered_psi
  channel$flip_sign <- flip_sign
  channel$nonzero_lambda <- nonzero_lambda
  stopifnot_(
    has_splines || !(channel$has_varying || channel$has_varying_intercept),
    "Model for response variable {.var {y}} contains time-varying
     definitions but splines have not been defined."
  )
  sampling[[paste0("y_", y_name)]] <- ifelse_(
    dformula$family$name %in%
      c("gaussian", "gamma", "exponential", "beta", "student"),
    t(Y_out),
    Y_out
  )
  list(channel = channel, sampling = sampling)
}

initialize_multivariate_channel <- function(y, y_cg, y_name, cg_idx,
                                            channel, family, dims,
                                            channel_vars, sampling_vars) {
  j <- cg_idx[1L]
  sampling <- list()
  merge_has <- c(
    "has_fixed_intercept",
    "has_varying_intercept",
    "has_random_intercept",
    "has_fixed",
    "has_varying",
    "has_random",
    "has_lfactor"
  )
  for (has in merge_has) {
    channel[[has]] <- unname(vapply(
      channel_vars[cg_idx], "[[", logical(1L), has
    ))
  }
  channel$has_missing <- any(
    vapply(channel_vars[cg_idx], "[[", logical(1L), "has_missing")
  )
  channel$has_fully_missing <- any(
    vapply(channel_vars[cg_idx], "[[", logical(1L), "has_fully_missing")
  )
  channel$y <- unname(vapply(channel_vars[cg_idx], "[[", character(1L), "y"))
  channel$y_cg <- y_cg
  channel$obs <- ifelse_(
    channel$has_missing,
    glue::glue("obs_{y_cg}[1:n_obs_{y_cg}[t], t]"),
    ""
  )
  z <- y_name[j]
  for (spec in formula_special_funs) {
    has_spec <- paste0("has_", spec)
    if (channel_vars[[z]][[has_spec]]) {
      channel[[has_spec]] <- TRUE
      sampling[[paste0(spec, "_", y_cg)]] <-
        sampling_vars[[paste0(spec, "_", z)]]
    } else {
      channel[[has_spec]] <- FALSE
    }
  }
  sampling[[paste0("obs_", y_cg)]] <- matrix_intersect(
    sampling_vars[paste0("obs_", y_name)]
  )
  sampling[[paste0("n_obs_", y_cg)]] <- apply(
    sampling[[paste0("obs_", y_cg)]],
    2L,
    function(x) { sum(x > 0L) }
  )
  sampling[[paste0("t_obs_", y_cg)]] <- which(
    sampling[[paste0("n_obs_", y_cg)]] > 0L
  )
  sampling[[paste0("T_obs_", y_cg)]] <- length(
    sampling[[paste0("t_obs_", y_cg)]]
  )
  O <- length(cg_idx)
  sampling[[paste0("O_", y_cg)]] <- O
  sampling[[paste0("y_", y_cg)]] <- array(
    unlist(sampling_vars[paste0("y_", y_name)]),
    c(dims, O)
  )
  if (is_multinomial(family)) {
    copy_indices <- c(
      "K_fixed",
      "K_varying",
      "K_random",
      "K",
      "J_fixed",
      "J_varying",
      "J",
      "J_random",
      "L_fixed",
      "L_varying"
    )
    copy_channel <- setdiff(
      names(channel_vars[[z]]),
      names(channel)
    )
    for (var in copy_channel) {
      channel[[var]] <- channel_vars[[z]][[var]]
    }
    for (idx in copy_indices) {
      sampling[[paste0(idx, "_", y_cg)]] <-
        sampling_vars[[paste0(idx, "_", z)]]
    }
    for (has in merge_has) {
      channel[[has]] <- channel[[has]][1L]
    }
  }
  list(channel = channel, sampling = sampling)
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
#' @param priors \[`data.frame`]\cr Prior definitions, or `NULL`, in which case
#'   the default priors are used.
#'
#' @srrstats {RE1.2} Checks for expected types and classes along with other
#'   `prepare_channel_*` functions.
#' @noRd
prepare_channel_default <- function(y, Y, channel, sampling,
                                    mean_gamma, sd_gamma, mean_y, sd_y,
                                    priors, category = "") {
  if (is.null(priors)) {
    out <- default_priors(
      y, channel, mean_gamma, sd_gamma, mean_y, sd_y, category
    )
    priors <- out$priors
    channel$prior_distr <- out$prior_distributions
  } else {
    priors <- ifelse_(
      nzchar(category),
      priors[priors$response == y & priors$category == category, ],
      priors[priors$response == y, ]
      )
    channel$prior_distr <- list()
    types <- priors$type
    loop_types <- intersect(
      types,
      c("alpha", "tau_alpha", "sigma_lambda", "psi", "kappa", "zeta")
    )
    for (ptype in loop_types) {
      pdef <- priors[priors$type == ptype, ]
      channel$prior_distr[[paste0(ptype, "_prior_distr")]] <- pdef$prior
    }
    loop_types <- c("beta", "delta", "tau", "sigma_nu")
    for (ptype in loop_types) {
      ifelse_(
        ptype %in% types,
        channel$prior_distr <- c(
          channel$prior_distr,
          create_vectorized_prior(ptype, priors, channel)
        ),
        channel$prior_distr[[paste0("vectorized_", ptype)]] <- FALSE
      )
    }
  }
  channel$Ks <- stats::setNames(
    channel$K_random,
    ifelse_(
      nzchar(category),
      paste0(y, "_", category),
      y
    )
  )
  list(channel = channel, sampling = sampling, priors = priors)
}

#' @describeIn prepare_channel_default Prepare a Bernoulli Channel
#' @noRd
prepare_channel_bernoulli <- function(y, Y, channel, sampling,
                                      sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Bernoulli", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  stopifnot_(
    all(Y_obs %in% c(0L, 1L)),
    c(
      "Response variable {.var {y}} is invalid:",
      `x` = "Bernoulli family supports only 0/1 integers."
    )
  )
  prepare_channel_binomial(y, Y, channel, sampling, sd_x, resp_class, priors)
}

#' @describeIn prepare_channel_default Prepare a Beta Channel
#' @noRd
prepare_channel_beta <- function(y, Y, channel, sampling,
                                 sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Beta", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs <= 0.0) || any(Y_obs >= 1.0)) {
    abort_nonunit(y, "Beta", type = "values", call = rlang::caller_env())
  }
  sd_y <- 1.0
  mean_y <- ifelse_(
    ncol(Y) > 1L,
    mean(Y[1L, ], na.rm = TRUE),
    Y[1L]
  )
  mean_y <- stats::qlogis(pmin(0.99, pmax(0.01, mean_y)))
  if (!is.finite(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    sampling,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
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
    out$channel$prior_distr$phi_prior_distr <- phi_prior$prior
    out$priors <- rbind(out$priors, phi_prior)
  } else {
    priors <- priors[priors$response == y, ]
    pdef <- priors[priors$type == "phi", ]
    out$channel$prior_distr$phi_prior_distr <- pdef$prior
    defaults <- rbind(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      phi_prior
    )
    check_priors(priors, defaults)
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Binomial Channel
#' @noRd
prepare_channel_binomial <- function(y, Y, channel, sampling,
                                     sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Binomial", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0.0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(
      y,
      "Binomial",
      type = "integers",
      call = rlang::caller_env()
    )
  }
  mean_y <- ifelse_(
    ncol(Y) > 1L,
    mean(Y[1L, ], na.rm = TRUE),
    Y[1L]
  )
  mean_y <- stats::qlogis(pmin(0.99, pmax(0.01, mean_y)))
  if (!is.finite(mean_y)) {
    mean_y <- 0.0
  }
  trials <- sampling[[paste0("trials_", y)]]
  if (!is.null(trials)) {
    trial_idx <- which(Y_obs < trials, arr.ind = TRUE)
    stopifnot_(
      nrow(trial_idx) > 0L,
      "Invalid number of trials at time index {trial_idx[1, 1]} for group
       {trial_idx[1, 2]}."
    )
  }
  sd_y <- 1
  mean_y <- 0.0
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    sampling,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    priors
  )
  if (!is.null(priors)) {
    check_priors(
      out$priors,
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors
    )
  }
  out
}


#' @describeIn prepare_channel_default Prepare a Categorical Channel
#' @noRd
prepare_channel_categorical <- function(y, Y, channel, sampling,
                                        sd_x, resp_class, priors) {
  stopifnot_(
    "factor" %in% resp_class,
    c(
      "Response variable {.var {y}} is invalid:",
      `x` = "Categorical family supports only {.cls factor} variables."
    )
  )
  resp_levels <- attr(resp_class, "levels")
  S_y <- length(resp_levels)
  channel$S <- S_y
  channel$categories <- resp_levels
  sampling[[paste0("S_", y)]] <- S_y
  sd_y <- 1.0
  mean_y <- 0.0
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  outcat <- lapply(
    resp_levels[-1L],
    function(s) {
      prepare_channel_default(
        y,
        Y,
        channel,
        sampling,
        mean_gamma,
        sd_gamma,
        mean_y,
        sd_y,
        priors,
        category = s
      )
    }
  )
  out <- outcat[[1L]]
  out$priors <- rbindlist_(lapply(outcat, "[[", "priors"))
  out$channel$prior_distr <- lapply(outcat, function(x) x$channel$prior_distr)
  names(out$channel$prior_distr) <- resp_levels[-1L]
  if (!is.null(priors)) {
    defaults <- rbindlist_(
      lapply(
        resp_levels[-1L],
        function(s) {
          default_priors(
            y,
            channel,
            mean_gamma,
            sd_gamma,
            mean_y,
            sd_y,
            category = s
          )$priors
        }
      )
    )
    check_priors(out$priors, defaults)
  }
  out$channel$Ks <- ulapply(outcat, function(x) x$channel$Ks)
  list(channel = out$channel, sampling = out$sampling, priors = out$priors)
}

#' @describeIn prepare_channel_default Prepare a Cumulative Channel
#' @noRd
prepare_channel_cumulative <- function(y, Y, channel, sampling,
                                       sd_x, resp_class, priors) {
  stopifnot_(
    "factor" %in% resp_class,
    c(
      "Response variable {.var {y}} is invalid:",
      `x` = "Cumulative family supports only {.cls factor} variables."
    )
  )
  resp_levels <- attr(resp_class, "levels")
  S_y <- length(resp_levels)
  channel$S <- S_y
  channel$categories <- seq_len(S_y - 1)
  sampling[[paste0("S_", y)]] <- S_y
  sd_y <- 1
  mean_y <- 0.0
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  fixed_cutpoints <- channel$has_fixed_intercept
  if (fixed_cutpoints) {
    cutpoint_priors <- data.frame(
      parameter = paste0("cutpoint_", y, "_", seq_len(S_y - 1)),
      response = y,
      prior = "std_normal()",
      type = "cutpoint",
      category = seq_len(S_y - 1)
    )
  }
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    sampling,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    priors
  )

  if (fixed_cutpoints) {
    out$priors <- out$priors[out$priors$type != "alpha", ]
    out$channel$prior_distr$alpha_prior_distr <- NULL
    if (is.null(priors)) {
      out$channel$prior_distr$cutpoint_prior_distr <- cutpoint_priors$prior
      names(out$channel$prior_distr$cutpoint_prior_distr) <-
        cutpoint_priors$category
      out$priors <- rbind(cutpoint_priors, out$priors)
    } else {
      priors <- priors[priors$response == y, ]
      pdef <- priors[priors$type == "cutpoint", ]
      out$channel$prior_distr$cutpoint_prior_distr <- pdef$prior
      default_priors <-
        default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors
      default_priors <- default_priors[default_priors$type != "alpha", ]
      defaults <- rbind(
        cutpoint_priors,
        default_priors
      )
      check_priors(priors, defaults)
    }
  } else {
    if (!is.null(priors)) {
      check_priors(
        priors[priors$response == y, ],
        default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors
      )
    }
  }
  out
}

#' @describeIn prepare_channel_default Prepare an Exponential Channel
#' @noRd
prepare_channel_exponential <- function(y, Y, channel, sampling,
                                        sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Exponential", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs <= 0.0)) {
    abort_negative(
      y,
      "Exponential",
      type = "values",
      call = rlang::caller_env()
    )
  }
  sd_y <- 1.0
  mean_y <- ifelse_(
    ncol(Y) > 1L,
    mean(Y[1L, ], na.rm = TRUE),
    Y[1L]
  )
  mean_y <- log(pmax(0.1, mean_y))
  if (!is.finite(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    sampling,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    priors
  )
  if (!is.null(priors)) {
    check_priors(
      out$priors,
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors
    )
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Gamma channel
#' @noRd
prepare_channel_gamma <- function(y, Y, channel, sampling,
                                  sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Gamma", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0.0)) {
    abort_negative(y, "Gamma", type = "values", call = rlang::caller_env())
  }
  sd_y <- 1.0
  mean_y <- ifelse_(
    ncol(Y) > 1L,
    mean(Y[1L, ], na.rm = TRUE),
    Y[1L]
  )
  mean_y <- log(pmax(0.1, mean_y))
  if (!is.finite(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    sampling,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
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
    out$channel$prior_distr$phi_prior_distr <- phi_prior$prior
    out$priors <- rbind(out$priors, phi_prior)
  } else {
    priors <- priors[priors$response == y, ]
    pdef <- priors[priors$type == "phi", ]
    out$channel$prior_distr$phi_prior_distr <- pdef$prior
    defaults <- rbind(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      phi_prior
    )
    check_priors(priors, defaults)
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Gaussian Channel
#' @noRd
prepare_channel_gaussian <- function(y, Y, channel, sampling,
                                     sd_x, resp_class, priors, ...) {
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
  if (!is.finite(sd_y) || identical(sd_y, 0.0)) {
    sd_y <- 1.0
  }
  sd_y <- max(1.0, sd_y)
  if (!is.finite(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 * sd_y / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    sampling,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
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
    out$channel$prior_distr$sigma_prior_distr <- sigma_prior$prior
    out$priors <- rbind(out$priors, sigma_prior)
  } else {
    priors <- priors[priors$response == y, ]
    pdef <- priors[priors$type == "sigma", ]
    out$channel$prior_distr$sigma_prior_distr <- pdef$prior
    defaults <- rbind(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      sigma_prior
    )
    check_priors(priors, defaults)
  }
  out
}

prepare_channel_multinomial <- function(y, y_cg, Y, channel, sampling,
                                        sd_x, resp_class, priors) {

  if (any("factor" %in% unlist(resp_class))) {
    abort_factor(y_cg, "Multinomial", call = rlang::caller_env())
  }
  obs <- sampling[[paste0("n_obs_", y_cg)]] > 0L
  Y_obs <- Y[obs, , ,drop = FALSE]
  if (any(Y_obs < 0.0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(
      y_cg,
      "Multinomial",
      type = "integers",
      call = rlang::caller_env()
    )
  }
  trials <- sampling[[paste0("trials_", y_cg)]][obs, , drop = FALSE]
  if (any(obs)) {
    trial_idx <- which(
      apply(Y, c(1L, 2L), sum, na.rm = TRUE) < trials,
      arr.ind = TRUE
    )
    stopifnot_(
      nrow(trial_idx) == 0L,
      "Invalid number of trials at time index {trial_idx[1, 1]} for group
       {trial_idx[1, 2]}."
    )
  }
  S_y <- dim(Y)[3L]
  channel$S <- S_y
  sampling[[paste0("S_", y_cg)]] <- S_y
  sd_y <- 1.0
  mean_y <- 0.0
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  outcat <- lapply(
    y[-1L],
    function(s) {
      prepare_channel_default(
        s,
        Y,
        channel,
        sampling,
        mean_gamma,
        sd_gamma,
        mean_y,
        sd_y,
        priors
      )
    }
  )
  out <- outcat[[1L]]
  out$priors <- rbindlist_(lapply(outcat, "[[", "priors"))
  out$channel$prior_distr <- lapply(outcat, function(x) x$channel$prior_distr)
  names(out$channel$prior_distr) <- y[-1L]
  if (!is.null(priors)) {
    defaults <- rbindlist_(
      lapply(
        y[-1L],
        function(s) {
          default_priors(s, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors
        }
      )
    )
    check_priors(out$priors, defaults)
  }
  out$channel$Ks <- ulapply(outcat, function(x) x$channel$Ks)
  list(channel = out$channel, sampling = out$sampling, priors = out$priors)
}

#' @describeIn prepare_channel_default Prepare a Multivariate Gaussian Channel
#' @noRd
prepare_channel_mvgaussian <- function(y_cg, channel, sampling, priors, ...) {
  L_prior <- data.frame(
    parameter = paste0("L_", y_cg),
    response = y_cg,
    prior = "lkj_corr_cholesky(1)",
    type = "L",
    category = ""
  )
  if (is.null(priors)) {
    mvpriors <- L_prior
    channel$prior_distr$L_prior_distr <- L_prior$prior
  } else {
    mvpriors <- priors[priors$response == y_cg, ]
    pdef <- priors[priors$type == "L", ]
    stopifnot_(
      identical(nrow(pdef), 1L),
      c(
        "Argument {.var priors} must contain all relevant parameters:",
        `x` = "Prior for parameter {.var L_{y_cg}} is not defined."
      )
    )
    # TODO some checks that prior distr makes sense
    channel$prior_distr$L_prior_distr <- pdef$prior
  }
  list(channel = channel, sampling = sampling, mvpriors = mvpriors)
}

#' @describeIn prepare_channel_default Prepare a Negative Binomial Channel
#' @noRd
prepare_channel_negbin <- function(y, Y, channel, sampling,
                                   sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Negative binomial", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0.0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(
      y,
      "Negative binomial",
      type = "integers",
      call = rlang::caller_env()
    )
  }
  sd_y <- 1.0
  mean_y <- ifelse_(
    ncol(Y) > 1L,
    mean(Y[1L, ], na.rm = TRUE),
    Y[1L]
  )
  mean_y <- log(pmax(0.1, mean_y))
  if (!is.finite(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    sampling,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
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
    out$channel$prior_distr$phi_prior_distr <- phi_prior$prior
    out$priors <- rbind(out$priors, phi_prior)
  } else {
    priors <- priors[priors$response == y, ]
    pdef <- priors[priors$type == "phi", ]
    out$channel$prior_distr$phi_prior_distr <- pdef$prior
    defaults <- rbind(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      phi_prior
    )
    check_priors(priors, defaults)
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Poisson channel
#' @noRd
prepare_channel_poisson <- function(y, Y, channel, sampling,
                                    sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Poisson", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0.0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(y, "Poisson", type = "integers", call = rlang::caller_env())
  }
  sd_y <- 1.0
  mean_y <- ifelse_(
    ncol(Y) > 1L,
    mean(Y[1L, ], na.rm = TRUE),
    Y[1L]
  )
  mean_y <- log(pmax(0.1, mean_y))
  if (!is.finite(mean_y)) {
    mean_y <- 0.0
  }

  sd_gamma <- 2.0 / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    sampling,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    priors
  )
  if (!is.null(priors)) {
    check_priors(
      out$priors,
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors
    )
  }
  out
}

#' @describeIn prepare_channel_default Prepare a Student-t Channel
#' @noRd
prepare_channel_student <- function(y, Y, channel, sampling,
                                    sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Student-t", call = rlang::caller_env())
  }
  if (ncol(Y) > 1L) {
    sd_y <- mean(apply(Y, 1L, sd, na.rm = TRUE))
    mean_y <- mean(Y[1L, ], na.rm = TRUE)
  } else {
    sd_y <- sd(Y, na.rm = TRUE)
    mean_y <- Y[1L]
  }
  if (!is.finite(sd_y) || identical(sd_y, 0.0)) {
    sd_y <- 1.0
  }
  sd_y <- max(1.0, sd_y)
  if (!is.finite(mean_y)) {
    mean_y <- 0.0
  }
  sd_gamma <- 2.0 * sd_y / sd_x
  mean_gamma <- rep(0.0, length(sd_gamma))
  out <- prepare_channel_default(
    y,
    Y,
    channel,
    sampling,
    mean_gamma,
    sd_gamma,
    mean_y,
    sd_y,
    priors
  )
  sigma_prior <- data.frame(
    parameter = paste0("sigma_", y),
    response = y,
    prior = paste0("exponential(", signif(1.0 / sd_y, 2L), ")"),
    type = "sigma",
    category = ""
  )
  phi_prior <- data.frame(
    parameter = paste0("phi_", y),
    response = y,
    prior = "gamma(2, 0.1)",
    type = "phi",
    category = ""
  )
  if (is.null(priors)) {
    out$channel$prior_distr$sigma_prior_distr <- sigma_prior$prior
    out$channel$prior_distr$phi_prior_distr <- phi_prior$prior
    out$priors <- rbind(out$priors, sigma_prior, phi_prior)
  } else {
    priors <- priors[priors$response == y, ]
    pdef <- priors[priors$type == "sigma", ]
    out$channel$prior_distr$sigma_prior_distr <- pdef$prior
    pdef <- priors[priors$type == "phi", ]
    out$channel$prior_distr$phi_prior_distr <- pdef$prior
    defaults <- rbind(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      sigma_prior,
      phi_prior
    )
    check_priors(priors, defaults)
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

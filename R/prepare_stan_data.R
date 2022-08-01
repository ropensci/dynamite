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
#' @param priors \[`data.frame`]\cr A data frame containing the prior
#'   definitions, or `NULL`, in which case default priors are used.
#' @param fixed \[`integer(1)`]\cr Number of fixed time points
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
prepare_stan_data <- function(data, dformula, group_var, time_var,
                              priors = NULL, fixed, verbose) {

  resp_names <- get_responses(dformula)
  missing_resp <- !(resp_names %in% names(data))
  stopifnot_(
    all(!missing_resp),
    "Can't find variable{?s} {.var {resp_names[missing_resp]}} in {.arg data}."
  )
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
  T_idx <- seq.int(fixed + 1L, T_full)
  groups <- !is.null(group_var)
  group <- onlyif(groups, data[[group_var]])
  has_splines <- !is.null(attr(dformula, "splines"))
  spline_defs <- prepare_splines(
    attr(dformula, "splines"),
    n_channels,
    T_idx
  )
  sampling_vars$D <- spline_defs$D
  sampling_vars$Bs <- spline_defs$Bs
  has_random <- attr(dformula, "random")$channels
  N <- ifelse_(groups, length(unique(group)), 1L)
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
  resp_classes <- attr(responses, "resp_class")
  for (i in seq_len(n_channels)) {
    channel <- list()
    resp <- resp_names[i]
    resp_split <- ifelse_(
      groups,
      split(responses[, resp], group),
      responses[, resp]
    )
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
    channel$K <- length(assigned[[i]])
    channel$K_fixed <- length(fixed_pars[[i]])
    channel$K_varying <- length(varying_pars[[i]])
    obs_idx <- array(0, dim = c(N, T_full - fixed))
    obs_len <- integer(T_full - fixed)
    for (j in seq_len(T_full - fixed)) {
      x_na <- X_na[j, , channel$J, drop = FALSE]
      dim(x_na) <- c(N, channel$K)
      y_na <- Y_na[j, ]
      obs_XY <- which(apply(x_na, 1, function(z) all(!z)) & !y_na)
      obs_XY_len <- length(obs_XY)
      obs_idx[, j] <- c(obs_XY, rep(0L, N - obs_XY_len))
      obs_len[j] <- obs_XY_len
    }
    channel$has_missing <- any(obs_len < N)
    sampling_vars[[paste0("obs_", resp)]] <- obs_idx
    sampling_vars[[paste0("n_obs_", resp)]] <- obs_len
    channel$obs <- ifelse_(
      channel$has_missing,
      glue::glue("obs_{resp}[1:n_obs_{resp}[t], t]"),
      ""
    )
    channel$has_fixed_intercept <- dformula[[i]]$has_fixed_intercept
    channel$has_varying_intercept <- dformula[[i]]$has_varying_intercept
    channel$has_random_intercept <- resp %in% has_random
    channel$has_fixed <- channel$K_fixed > 0L
    channel$has_varying <- channel$K_varying > 0L
    channel$lb <- spline_defs$lb[i]
    channel$shrinkage <- spline_defs$shrinkage
    channel$noncentered <- spline_defs$noncentered[i]
    if (!has_splines &&
        (channel$has_varying || channel$has_varying_intercept)) {
      stop_("Model for response variable {.var {resp}}
             contains time-varying definitions
             but splines have not been defined.")
    }
    for (spec in formula_special_funs) {
      if (!is.null(form_specials[[spec]])) {
        spec_split <- ifelse_(
          groups,
          split(form_specials[[spec]], group),
          form_specials[[spec]]
        )
        spec_array <- array(as.numeric(unlist(spec_split)), dim = c(T_full, N))
        sampling_vars[[paste0(spec, "_", resp)]] <-
          spec_array[seq.int(fixed + 1, T_full), , drop = FALSE]
        channel[[paste0("has_", spec)]] <- TRUE
      } else {
        channel[[paste0("has_", spec)]] <- FALSE
      }
    }
    family <- dformula[[i]]$family
    sampling_vars[[resp]] <- ifelse_(
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
  sampling_vars$M <- sum(unlist(lapply(model_vars, "[[",
    "has_random_intercept")))
  # avoid goodpractice warning, T is a Stan variable, not an R variable
  sampling_vars[["T"]] <- T_full - fixed
  sampling_vars$X_m <- as.array(x_means)

  if (spline_defs$shrinkage) {
    prior_list[["common_priors"]] <- ifelse_(
      is.null(priors),
      data.frame(
        parameter = "lambda",
        response = "",
        prior = "normal(0, 1)",
        type = "lambda",
        category = ""
      ),
      priors |> dplyr::filter(.data$type == "lambda")
    )
  }
  if (sampling_vars$M > 1 && attr(dformula, "random")$correlated) {
    prior_list[["common_priors"]] <- ifelse_(
      is.null(priors),
      dplyr::bind_rows(
        prior_list[["common_priors"]],
        data.frame(
          parameter = "L",
          response = "",
          prior = "lkj_corr_cholesky(1)",
          type = "L",
          category = ""
        )
      ),
      dplyr::bind_rows(
        prior_list[["common_priors"]],
        priors |> dplyr::filter(.data$type == "L")
      )
    )
  }
  # for stanblocks
  attr(model_vars, "common_priors") <- prior_list[["common_priors"]]
  list(
    model_vars = model_vars,
    sampling_vars = sampling_vars,
    priors = prior_list,
    responses = responses,
    u_names = colnames(model_matrix),
    fixed = fixed
  )
}

#' Prepare B-spline Parameters for Stan
#'
#' @param spline_defs A `splines` object.
#' @param n_channels Number of channels
#' @param T_idx An `integer` vector of time indices
#' @noRd
prepare_splines <- function(spline_defs, n_channels, T_idx) {
  out <- list()
  if (!is.null(spline_defs)) {
    out$shrinkage <- spline_defs$shrinkage
    out$bs_opts <- spline_defs$bs_opts
    out$bs_opts$x <- T_idx
    if (is.null(out$bs_opts$Boundary.knots)) {
      out$bs_opts$Boundary.knots <- range(out$bs_opts$x)
    }
    out$Bs <- t(do.call(splines::bs, args = out$bs_opts))
    out$D <- nrow(out$Bs)
    out$lb <- spline_defs$lb_tau
    if (length(out$lb) %in% c(1L, n_channels)) {
      out$lb <- rep(out$lb, length = n_channels)
    } else {
      stop_(
        "Length of the {.arg lb_tau} argument of {.fun splines} function
          is not equal to 1 or {n_channels}, the number of the channels."
      )
    }
    out$noncentered <- spline_defs$noncentered
    if (length(out$noncentered) %in% c(1L, n_channels)) {
      out$noncentered <- rep(out$noncentered, length = n_channels)
    } else {
      stop_(
        "Length of the {.arg noncentered} argument of {.fun splines} function
          is not equal to 1 or {n_channels}, the number of the channels."
      )
    }
  } else {
    out <- list(
      lb = numeric(n_channels),
      noncentered = logical(n_channels),
      shrinkage = logical(1)
    )
  }
  out
}

#' Construct a Prior Definition for a Time-varying Parameter
#'
#' @param ptype Type of the parameter
#' @param priors A `data.frame` defining the priors
#' @param channel A `list` of channel-specific variables for Stan sampling
#' @noRd
prepare_prior <- function(ptype, priors, channel) {
  pdef <- priors |> dplyr::filter(.data$type == ptype)
  channel[[paste0(ptype, "_prior_distr")]] <- pdef$prior
  dists <- sub("\\(.*", "", pdef$prior)
  if (nrow(pdef) > 0L && length(unique(dists)) == 1L) {
    pars <- strsplit(sub(".*\\((.*)\\).*", "\\1", pdef$prior), ",")
    pars <- do.call("rbind", lapply(pars, as.numeric))
    channel[[paste0(ptype, "_prior_npars")]] <- ncol(pars)
    channel[[paste0(ptype, "_prior_pars")]] <- pars
    channel[[paste0(ptype, "_prior_distr")]] <- dists[1L]
  }
  channel
}

#' Default channel preparation
#'
#' Computes default channel-specific variables for Stan sampling,
#' Stan model code construction, and prior definitions.
#'
#' @param y \[`character(1)`]\cr Name of the response variable of the channel.
#' @param Y \[`matrix()`]\cr A matrix of values of the response variable.
#' @param channel \[`list()`]\cr Channel-specific helper variables.
#' @param mean_gamma Prior mean betas and deltas (at time `fixed + 1`).
#' @param sd_gamma Prior SD betas and deltas (at time `fixed + 1`).
#' @param mean_y Mean of the response variable at time `fixed + 1`.
#' @param sd_y SD of the response variable at time `fixed + 1`.
#' @param resp_class \[`character()`]\cr Class(es) of the response `Y`.
#' @param priors \[`data.frame`]\cr A data frame containing the prior
#'   definitions, or `NULL`, in which case default priors are used.
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
    priors <- priors |> dplyr::filter(.data$response == y)
    for (ptype in c("alpha", "tau_alpha", "sigma_nu")) {
      pdef <- priors |> dplyr::filter(.data$type == ptype)
      channel[[paste0(ptype, "_prior_distr")]] <- pdef$prior
    }
    for (ptype in c("beta", "delta", "tau")) {
      channel <- prepare_prior(ptype, priors, channel)
    }
  }
  channel$write_beta <- channel$has_fixed &&
    length(channel$beta_prior_distr) == 1L
  channel$write_delta <- channel$has_varying &&
    length(channel$delta_prior_distr) == 1L
  channel$write_tau <- channel$has_varying &&
    length(channel$tau_prior_distr) == 1L
  list(channel = channel, priors = priors)
}

#' @describeIn prepare_channel_default Prepare a categorical channel
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
    priors <- priors |> dplyr::filter(.data$response == y)
    for (ptype in c("alpha", "beta", "delta", "tau")) {
      channel <- prepare_prior(ptype, priors, channel)
    }
    pdef <- priors |> dplyr::filter(.data$type == "tau_alpha")
    channel[["tau_alpha_prior_distr"]] <- pdef$prior
    priors <- check_priors(
      priors, default_priors_categorical(y, channel, sd_x, resp_class)$priors
    )
  }
  channel$write_alpha <-
    (channel$has_fixed_intercept || channel$has_varying_intercept) &&
    length(channel$alpha_prior_distr) == 1L
  channel$write_beta <- channel$has_fixed &&
    length(channel$beta_prior_distr) == 1L
  channel$write_delta <- channel$has_varying &&
    length(channel$delta_prior_distr) == 1L
  channel$write_tau <- channel$has_varying &&
    length(channel$tau_prior_distr) == 1L
  list(channel = channel, priors = priors)
}

#' @describeIn prepare_channel_default Prepare a gaussian channel
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
  if (is.na(sd_y) || sd_y == 0) {
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
    prior = paste0("exponential(", signif(1.0 / sd_y, 2), ")"),
    type = "sigma",
    category = ""
  )
  if (is.null(priors)) {
    out$channel$sigma_prior_distr <- sigma_prior$prior
    out$priors <- dplyr::bind_rows(out$priors, sigma_prior)
  } else {
    priors <- priors |> dplyr::filter(.data$response == y)
    pdef <- priors |> dplyr::filter(.data$type == "sigma")
    if (nrow(pdef) == 1L) {
      out$channel$sigma_prior_distr <- pdef$prior
    }
    defaults <- dplyr::bind_rows(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      sigma_prior)
    out$priors <- check_priors(priors, defaults)
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
  if (any(Y_obs < 0.0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(y, "Binomial", type = "integers",
                   call = rlang::caller_env())
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
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors)
  }
  out
}

#' @describeIn prepare_channel_default Prepare a bernoulli channel
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
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors)
  }
  out
}

#' @describeIn prepare_channel_default Prepare a negative binomial channel
#' @noRd
prepare_channel_negbin <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Negative binomial", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0.0) || any(Y_obs != as.integer(Y_obs))) {
    abort_negative(y, "Negative binomial", type = "integers",
                   call = rlang::caller_env())
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
    out$priors <- dplyr::bind_rows(out$priors, phi_prior)
  } else {
    priors <- priors |> dplyr::filter(.data$response == y)
    pdef <- priors |> dplyr::filter(.data$type == "phi")
    if (nrow(pdef) == 1L) {
      out$channel$phi_prior_distr <- pdef$prior
    }
    defaults <- dplyr::bind_rows(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      phi_prior)
    out$priors <- check_priors(priors, defaults)
  }
  out
}

#' @describeIn prepare_channel_default Prepare an exponential channel
#' @noRd
prepare_channel_exponential <- function(y, Y, channel, sd_x, resp_class,
                                        priors) {
  if ("factor" %in% resp_class) {
    abort_factor(y, "Exponential", call = rlang::caller_env())
  }
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs <= 0.0)) {
    abort_negative(y, "Exponential", type = "values",
                   call = rlang::caller_env())
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
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors)
  }
  out

}

#' @describeIn prepare_channel_default Prepare a gamma channel
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
    out$priors <- dplyr::bind_rows(out$priors, phi_prior)
  } else {
    priors <- priors |> dplyr::filter(.data$response == y)
    pdef <- priors |> dplyr::filter(.data$type == "phi")
    if (nrow(pdef) == 1L) {
      out$channel$phi_prior_distr <- pdef$prior
    }
    defaults <- dplyr::bind_rows(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      phi_prior)
    out$priors <- check_priors(priors, defaults)
  }
  out
}

#' @describeIn prepare_channel_default Prepare a beta channel
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
    out$priors <- dplyr::bind_rows(out$priors, phi_prior)
  } else {
    priors <- priors |> dplyr::filter(.data$response == y)
    pdef <- priors |> dplyr::filter(.data$type == "phi")
    if (nrow(pdef) == 1L) {
      out$channel$phi_prior_distr <- pdef$prior
    }
    defaults <- dplyr::bind_rows(
      default_priors(y, channel, mean_gamma, sd_gamma, mean_y, sd_y)$priors,
      phi_prior)
    out$priors <- check_priors(priors, defaults)
  }
  out
}
#' Raise an error if factor type is not supported by a family
#'
#' @param y Response variable the error is related to.
#' @param family Family as character vector.
#' @param call call to be passed to [stop_()].
#' @noRd
abort_factor <- function(y, family, call) {
  stop_(c(
    "Response variable {.var {y}} is invalid:",
    `x` = "{family} family is not supported for {.cls factor} variables."
  ), call = call)
}

#' Raise error if negative values are not supported by a family
#'
#' @param y Response variable the error is related to.
#' @param family Family as character vector.
#' @param type Value type of the family.
#' @param call call to be passed to [stop_()].
#' @noRd
abort_negative <- function(y, family, type, call) {
  stop_(c(
    "Response variable {.var {y}} is invalid:",
    `x` = "{family} family supports only non-negative {type}."
  ), call = call)
}

#' Raise error if values outside unit interval
#'
#' @param y Response variable the error is related to.
#' @param family Family as character vector.
#' @param type Value type of the family.
#' @param call call to be passed to [stop_()].
#' @noRd
abort_nonunit <- function(y, family, type, call) {
  stop_(c(
    "Response variable {.var {y}} is invalid:",
    `x` = "{family} family supports only {type} on open interval (0, 1)."
  ), call = call)
}

#' Prepare data for Stan
#'
#' Prepares data for Stan sampling, Stan model code construction and
#' default/user-modifiable prior definitions.
#'
#' @param e \[`environment`]\cr An environment containing the data
#' @param dformula \[`dynamiteformula`]\cr The model formula of stochastic
#'   channels
#' @param group_var \[`character(1)`]\cr The grouping variable name
#' @param time_var \[`character(1)`]\cr The time index variable name
#' @param priors TODO
#'
#' @noRd
prepare_stan_data <- function(e, dformula, group_var, time_var, priors = NULL) {

  responses <- e$data[, get_responses(dformula), drop = FALSE]
  # Needs sapply/lapply instead of apply to keep factors as factors
  attr(responses, "resp_class") <- lapply(responses, function(x) {
    cl <- class(x)
    attr(cl, "levels") <- levels(x)
    cl
  })
  model_matrix <- full_model.matrix(dformula, e$data)
  specials <- evaluate_specials(dformula, e$data)
  resp_names <- colnames(responses)
  n_channels <- length(resp_names)
  # A list of variables for stan sampling without grouping by channel
  sampling_vars <- list()
  empty_list <- setNames(vector(mode = "list", length = n_channels), resp_names)
  # A list containing a list for each channel consisting of
  # variables used to construct the stan model code
  model_vars <- empty_list
  # A list for getting current prior definitions
  prior_list <- empty_list

  group <- e$data[[group_var]]
  time <- sort(unique(e$data[[time_var]]))
  T_full <- length(time)
  groups <- !is.null(group)

  spline_defs <- attr(dformula, "splines")
  has_splines <- !is.null(spline_defs)
  lb <- ""
  shrinkage <- FALSE
  if (has_splines) {
    lb <- spline_defs$lb_tau
    shrinkage <- spline_defs$shrinkage
    bs_opts <- spline_defs$bs_opts
    bs_opts$x <- time
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
        "Length of the 'noncentered' argument of 'splines' function ",
        "is not equal to 1 or the number of the channels. Recycling."
      )
    }
    sampling_vars$D <- D
    sampling_vars$Bs <- Bs
  }
  N <- T_full
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
  sd_x <- setNames(pmax(0.5, apply(X, 3, sd, na.rm = TRUE)),
                   colnames(model_matrix))
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
    channel$L_fixed <- as.array(fixed_pars[[i]])
    channel$L_varying <- as.array(varying_pars[[i]])
    channel$J <- as.array(assigned[[i]])
    channel$J_fixed <- as.array(assigned[[i]][channel$L_fixed])
    channel$J_varying <- as.array(assigned[[i]][channel$L_varying])
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
    channel$has_fixed <- channel$K_fixed > 0
    channel$has_varying <- channel$K_varying > 0
    channel$lb <- lb
    channel$shrinkage <- shrinkage
    if (channel$has_varying) {
      if (!has_splines) {
        stop_("Model for response variable ", resp, " ",
              "contains time-varying definitions ",
              "but splines have not been defined.")
        # TODO switch back to warning after defining default splines?
        warn_nosplines <- TRUE
      }
      if (warn_nosplines) {
        channel$has_varying <- FALSE
        channel$noncentered <- FALSE
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
    if (is_gaussian(family)) {
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
  # TODO move this before prep
  if (warn_nosplines) {
    warning_("All channels will now default to time-constant ",
             "coefficients for all predictors.")
  }
  sampling_vars$N <- N
  sampling_vars$K <- K
  sampling_vars$X <- X
  sampling_vars$T <- T_full
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
#' @param Y \[`vector()`]\cr A vector of values of the response variable.
#' @param channel \[`list()`]\cr Channel-specific helper variables.
#' @param sd_gamma TODO
#' @param resp_class \[`character()`]\cr Class(es) of the response `Y`.
#' @param priors TODO
#'
#' @noRd
prepare_channel_default <- function(y, Y, channel,
                                    mean_gamma, sd_gamma, resp_class, priors) {
  if (is.null(priors)) {
    priors <- list()

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

  S_y <- length(unique(na.exclude(as.vector(Y))))
  channel$S <- S_y
  if (!("factor" %in% resp_class)) {
    stop_("Response variable ", y, " is invalid: ",
          "categorical family supports only factors.")
  }
  if (is.null(priors)) {
    # remove the first level which acts as reference
    resp_levels <- attr(resp_class, "levels")[-1]
    priors <- list()
    sd_gamma <- 2 / sd_x
    # TODO arbitrary, perhaps should depend on S
    sd_gamma[which(names(sd_gamma) == "(Intercept)")] <- 5

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
    for (ptype in c("beta", "delta", "tau")) {
      pdef <- priors |> dplyr::filter(.data$type == ptype)
      if (nrow(pdef) > 0) {
        dists <- sub("\\(.*", "", pdef$prior)
        vectorized <- length(unique(dists)) == 1
        if (vectorized) {
          pars <- strsplit(sub(".*\\((.*)\\).*", "\\1", pdef$prior), ",")
          pars <- do.call("rbind", lapply(pars, as.numeric))
          channel[[paste0(ptype, "_prior_npars")]] <- ncol(pars)
          channel[[paste0(ptype, "_prior_parsy")]] <- pars
          channel[[paste0(ptype, "_prior_distr")]] <- dists[1]
        } else {
          channel[[paste0(ptype, "_prior_distr")]] <- pdef$prior # write separate priors
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

#' @describeIn prepare_channel_default Prepare a gaussian channel
#' @noRd
prepare_channel_gaussian <- function(y, Y, channel, sd_x, resp_class, priors) {
  if ("factor" %in% resp_class) {
    stop_("Response variable ", y, " is invalid: ",
          "gaussian family is not supported for factors.")
  }
  if (ncol(Y) > 1) {
    s_y <- 1 / mean(apply(Y, 1, sd, na.rm = TRUE))
  } else {
    s_y <- sd(Y, na.rm = TRUE)
  }
  sd_gamma <- 2 * s_y / sd_x
  mean_gamma <- rep(0, length(sd_gamma))
  idx <- which(names(sd_gamma) == "(Intercept)")
  sd_gamma[idx] <- 2 * s_y
  mean_gamma[idx] <- mean(Y)

  out <- prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma,
                                 resp_class, priors)
  if (is.null(priors)) {

    out$channel$sigma_prior_distr <- paste0("exponential(", s_y, ")")
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
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0) ||
      any(is.logical(Y_obs)) ||
      any(Y_obs != as.integer(Y_obs))) {
    stop_("Response variable ", y, " is invalid: ",
          "binomial family supports only non-negative integers.")
  }
  if ("factor" %in% resp_class) {
    stop_("Response variable ", y, " is invalid: ",
          "binomial family is not supported for factors.")
  }
  sd_gamma <- 2 / sd_x
  sd_gamma[which(names(sd_gamma) == "(Intercept)")] <- 2.5
  mean_gamma <- rep(0, length(sd_gamma))
  prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma, resp_class,
                          priors)
}

#' @describeIn prepare_channel_default Prepare a bernoulli channel
#' @noRd
prepare_channel_bernoulli <- function(y, Y, channel, sd_x, resp_class,
                                      priors) {
  # TODO should bernoulli autoconvert logical to integer?
  # Maybe, if so, then autoconversions in categorical should be done as well?
  Y_obs <- Y[!is.na(Y)]
  if (!all(Y_obs %in% 0:1) || any(is.logical(Y_obs))) {
    stop_("Response variable ", y, " is invalid: ",
          "bernoulli family supports only 0/1 integers.")
  }
  if ("factor" %in% resp_class) {
    stop_("Response variable ", y, " is invalid: ",
          "bernoulli family is not supported for factors.")
  }
  prepare_channel_binomial(y, Y, channel, sd_x, resp_class, priors)
}

#' @describeIn prepare_channel_default Prepare a Poisson channel
#' @noRd
prepare_channel_poisson <- function(y, Y, channel, sd_x, resp_class, priors) {
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0) || any(Y_obs != as.integer(Y_obs))) {
    stop_("Response variable ", y, " is invalid: ",
          "Poisson family supports only non-negative integers.")
  }
  if ("factor" %in% resp_class) {
    stop_("Response variable ", y, " is invalid: ",
          "Poisson family is not supported for factors.")
  }
  sd_gamma <- 2 / sd_x
  sd_gamma[which(names(sd_gamma) == "(Intercept)")] <- 10
  mean_gamma <- rep(0, length(sd_gamma))
  prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma, resp_class,
                          priors)
}

#' @describeIn prepare_channel_default Prepare a negative binomial channel
#' @noRd
prepare_channel_negbin <- function(y, Y, channel, sd_x, resp_class, priors) {
  Y_obs <- Y[!is.na(Y)]
  if (any(Y_obs < 0) || any(Y_obs != as.integer(Y_obs))) {
    stop_("Response variable ", y, " is invalid: ",
          "negative binomial family supports only non-negative integers.")
  }
  if ("factor" %in% resp_class) {
    stop_("Response variable ", y, " is invalid: ",
          "negative binomial family is not supported for factors.")
  }
  sd_gamma <- 2 / sd_x
  sd_gamma[which(names(sd_gamma) == "(Intercept)")] <- 10
  mean_gamma <- rep(0, length(sd_gamma))
  out <- prepare_channel_default(y, Y, channel, mean_gamma, sd_gamma,
                                 resp_class, priors)
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
  warning_(
    "Found nonfinite prior standard deviation when using default priors ",
    "for regression coeffients for response ", y, " ",
    "indicating constant covariate: Switching to N(0, 0.01) prior."
  )
}

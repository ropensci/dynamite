# Convert data for Stan
convert_data <- function(formula, responses, specials, group, time,
    model_matrix, coef_names, priors = NULL) {

    # A list of variables for stan sampling without grouping by channel
    sampling_vars <- list()
    # A list containing a list for each channel consisting of variables used to construct the stan model code
    model_vars <- list()
    # A list for getting current prior definitions
    prior_list <- list()

    T_full <- length(time)
    groups <- !is.null(group)
    #free_obs <- (fixed + 1):T_full
    n_channels <- length(get_resp(formula))
    spline_defs <- attr(formula, "splines")
    has_splines <- !is.null(spline_defs)
    lb <- ""
    shrinkage <- FALSE
    if (has_splines) {
        lb <- attr(formula, "splines")$lb_tau
        shrinkage <- attr(formula, "splines")$shrinkage
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
            warning_(paste0(
                "Length of the 'noncentered' argument of 'splines' function ",
                "is not equal to 1 or the number of the channels. Recycling. "))
        }
        sampling_vars$D <- D
        sampling_vars$Bs <- Bs
    }
    N <- T_full
    if (groups) {
        N <- length(unique(group))
    }
    K <- ncol(model_matrix)
    X <- aperm(array(as.numeric(unlist(split(model_matrix, gl(T_full, 1, N * T_full)))),
                     dim = c(N, K, T_full)),
               c(3, 1, 2))
    X[is.na(X)] <- 0
    assigned <- attr(model_matrix, "assign")
    fixed_pars <- attr(model_matrix, "fixed")
    varying_pars <- attr(model_matrix, "varying")
    # if (N > 1) {
    #     sd_x <- apply(X[1, , ], 2, sd)
    # } else {
    #     sd_x <- apply(X[, 1, ], 2, sd)
    # }
    # sd_x[sd_x < 0.1] <- 0.1 # Intercept and other constants at time 1
    # use sd over all time points so we get reasonable prior scale for covariates which are constant at t=1
    sd_x <- pmax(0.5, apply(X, 3, sd, na.rm = TRUE))
    resp_classes <- attr(responses, "resp_class")
    warn_nosplines <- FALSE
    for (i in seq_len(n_channels)) {
        channel <- list()
        resp <- formula[[i]]$response
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
        obs_idx <- apply(X[,,channel$J], 1, function(x) {
            nc <- nrow(x)
            obs <- which(apply(x, 1, function(y) {
                all(!is.na(y))
            }))
            c(obs, rep(0, nc - length(obs)))
        })
        obs_len <- apply(obs_idx, 2, function(x) {
            sum(x > 0)
        })
        channel$has_missing <- any(obs_len < N)
        if (channel$has_missing) {
            sampling_vars[[paste0("obs_", resp)]] <- t(obs_idx)
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
                stop_("Model for response variable ", resp, " contains time-varying definitions, but splines have not been defined.")
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
        if (groups) {
            resp_split <- split(responses[, resp], group)
        } else {
            resp_split <- responses[, resp]
        }
        if (length(form_specials)) {
            for (spec in formula_special_funs) {
                if (!is.null(form_specials[[spec]])) {
                    if (groups) {
                        spec_split <- split(form_specials[[spec]], group)
                    } else {
                        spec_split <- form_specials[[spec]]
                    }
                    sampling_vars[[paste0(spec, "_", resp)]] <- array(as.numeric(unlist(spec_split)), dim = c(T_full, N))
                    channel[[paste0("has_", spec)]] <- TRUE
                } else {
                    channel[[paste0("has_", spec)]] <- FALSE
                }
            }
        }
        Y <- array(as.numeric(unlist(resp_split)), dim = c(T_full, N))
        family <- formula[[i]]$family
        if (is_gaussian(family)) {
            sampling_vars[[resp]] <- t(Y) # NxT matrix
        } else {
            sampling_vars[[resp]] <- Y # T*N array (needs to be integers)
        }
        prep <- do.call(paste0("prepare_channel_", family),
            list(
                y = resp,
                Y = Y,
                channel,
                sd_x = sd_x,
                resp_class = resp_classes[resp],
                coef_names = coef_names[[i]],
                priors = priors)
            )
        prior_list[[resp]] <- prep$priors
        model_vars[[i]] <- prep$channel
        sampling_vars <- c(sampling_vars, prep$sampling_vars)
    }
    # TODO move this before prep
    if (warn_nosplines) {
        warning_("All channels will now default to time-constant coefficients for all predictors.")
    }
    sampling_vars$N <- N
    sampling_vars$K <- K
    sampling_vars$X <- X
    sampling_vars$T <- T_full
    list(model_vars = model_vars, sampling_vars = sampling_vars, priors = prior_list)
}

prepare_channel_default <- function(y, Y, channel,
    sd_beta, resp_class, coef_names, priors) {

    if (is.null(priors)) {
        priors <- list()
        #default priors
        bnames <- gsub(paste0("^", y), "beta", coef_names)
        if (channel$has_fixed) {
            m <- rep(0, channel$K_fixed)
            s <- sd_beta[channel$J_fixed]
            channel$beta_fixed_prior_npars <- 2
            channel$beta_fixed_prior_pars <- cbind(m, s, deparse.level = 0)
            channel$beta_fixed_prior_distr <- "normal"

            priors$beta_fixed <-
                data.frame(parameter = bnames[channel$L_fixed],
                    response = y,
                    prior = paste0("normal(", m, ", ", s, ")"),
                    type = "beta_fixed",
                    category = "")
        }
        if (channel$has_varying) {
            m <- rep(0, channel$K_varying)
            s <- sd_beta[channel$J_varying]
            channel$beta_varying_prior_npars <- 2
            channel$beta_varying_prior_pars <- cbind(m, s, deparse.level = 0)
            channel$beta_varying_prior_distr <- "normal"

            priors$beta_varying <-
                data.frame(parameter = bnames[channel$L_varying],
                    response = y,
                    prior = paste0("normal(", m, ", ", s, ")"),
                    type = "beta_varying",
                    category = "")

            channel$tau_prior_npars <- 2
            channel$tau_prior_pars <- cbind(0, rep(1, channel$K_varying))
            channel$tau_prior_distr <- "normal"
            bnames <- gsub(paste0("^", y), "", coef_names)
            priors$tau <-
                data.frame(parameter = paste0("tau", bnames[channel$L_varying]),
                    response = y,
                    prior = "normal(0, 1)",
                    type = "tau",
                    category = "")
        }
        priors <- dplyr::bind_rows(priors)
    } else {
        # TODO add a warning to documentation that the only the 'prior' column
        # of the priors data.frame should be altered (i.e. there's no checks for names or reordering of rows)
        # Or arrange...
        priors <- priors |> dplyr::filter(response == y)
        for (ptype in c("beta_fixed", "beta_varying", "tau")) {
            pdef <- priors |> dplyr::filter(type == ptype)
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
    channel$write_beta_fixed <- channel$has_fixed && length(channel$beta_fixed_prior_distr) == 1
    channel$write_beta_varying <- channel$has_varying && length(channel$beta_varying_prior_distr) == 1
    channel$write_tau <- channel$has_varying && length(channel$tau_prior_distr) == 1
    list(channel = channel, priors = priors)
}

prepare_channel_categorical <- function(y, Y, channel,
    sd_x, resp_class, coef_names, priors) {

    S_y <- length(unique(na.exclude(as.vector(Y))))
    channel$S <- S_y
    if (!("factor" %in% resp_class)) {
        stop_("Response variable ", y, " is invalid: categorical family supports only factors.")
    }
    if (is.null(priors)) {
        priors <- list()
        #default priors
        bnames <- gsub(paste0("^", y), "beta", attr(coef_names, "simplified")$names)
        levels_ <- attr(coef_names, "simplified")$levels

        sd_beta <- 2 / sd_x
        k <- grep("(Intercept)", bnames)
        if (!is.null(k)) sd_beta[k] <- 5 # TODO arbitrary, perhaps should depend on S
        if (any(!is.finite(sd_beta))) { # never happens due to pmax(0.5,sd_x)
            msg <- paste0("Found nonfinite prior standard deviation when using default priors for regression coeffients for response ",
                y, ", indicating constant covariate: Switching to N(0, 0.01) prior.")
            sd_beta[!is.finite(sd_beta)] <- 0.1
            warning(msg)
        }

        if (channel$has_fixed > 0) {
            m <- rep(0, channel$K_fixed * (S_y - 1))
            s <- rep(sd_beta[channel$J_fixed], S_y - 1) # match with binomial in case S_y=2

            channel$beta_fixed_prior_npars <- 2
            channel$beta_fixed_prior_distr <- "normal"
            channel$beta_fixed_prior_pars <- cbind(m, s, deparse.level = 0)
            priors$beta_fixed <-
                data.frame(parameter = bnames[channel$L_fixed],
                    response = y,
                    prior = paste0("normal(", m, ", ", s, ")"),
                    type = "beta_fixed",
                    category = levels_)
        }

        if (channel$has_varying) {
            m <- rep(0, channel$K_varying * (S_y - 1))
            s <- rep(sd_beta[channel$J_varying], S_y - 1)
            channel$beta_varying_prior_npars <- 2
            channel$beta_varying_prior_pars <- cbind(m, s, deparse.level = 0)
            channel$beta_varying_prior_distr <- "normal"

            priors$beta_varying <-
                data.frame(parameter = bnames[channel$L_varying],
                    response = y,
                    prior = paste0("normal(", m, ", ", s, ")"),
                    type = "beta_varying",
                    category = levels_)

            channel$tau_prior_npars <- 2
            channel$tau_prior_pars <- cbind(0, rep(1, channel$K_varying))
            channel$tau_prior_distr <- "normal"
            bnames <- gsub(paste0("^", y), "", attr(coef_names, "simplified")$names)
            priors$tau <-
                data.frame(parameter = paste0("tau", bnames[channel$L_varying]),
                    response = y,
                    prior = "normal(0, 1)",
                    type = "tau",
                    category = "")
        }
        priors <- dplyr::bind_rows(priors)
    } else {
        # TODO add a warning to documentation that the only the 'prior' column
        # of the priors data.frame should be altered (i.e. there's no checks for names or reordering of rows)
        # Or arrange...
        priors <- priors |> dplyr::filter(response == y)
        for (ptype in c("beta_fixed", "beta_varying", "tau")) {
            pdef <- priors |> dplyr::filter(type == ptype)
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
    channel$write_beta_fixed <- channel$has_fixed && length(channel$beta_fixed_prior_distr) == 1
    channel$write_beta_varying <- channel$has_varying && length(channel$beta_varying_prior_distr) == 1
    channel$write_tau <- channel$has_varying && length(channel$tau_prior_distr) == 1
    list(channel = channel, priors = priors)
}

prepare_channel_gaussian <- function(y, Y, channel,
    sd_x, resp_class, coef_names, priors) {

    if ("factor" %in% resp_class) {
        stop_("Response variable ", y, " is invalid: gaussian family is not supported for factors.")
    }
    sd_beta <- 2 / sd_x
    k <- grep("(Intercept)", coef_names)
    if (!is.null(k)) sd_beta[k] <- 10 # Wider prior for intercept as we are not centering X
    if (any(!is.finite(sd_beta))) {
        msg <- paste0("Found nonfinite prior standard deviation when using default priors for regression coeffients for response ",
            y, ", indicating constant covariate: Switching to N(0, 0.01) prior.")
        sd_beta[!is.finite(sd_beta)] <- 0.1
        warning(msg)
    }
    out <- prepare_channel_default(y, Y, channel,
        sd_beta, resp_class, coef_names, priors)
    if (is.null(priors)) {
        s <- 1 / mean(apply(Y, 1, sd, na.rm = TRUE))
        out$channel$sigma_prior_distr <- paste0("exponential(", s, ")")
        out$priors <- dplyr::bind_rows(out$priors,
            data.frame(parameter = "sigma",
                response = y,
                prior = out$channel$sigma_prior_distr,
                type = "sigma",
                category = ""))
    } else {
        pdef <- priors |> dplyr::filter(response == y && type == "sigma")
        if (nrow(pdef) == 1) {
            out$channel$sigma_prior_distr <- pdef$prior
        }
    }
    out
}

prepare_channel_binomial <- function(y, Y, channel,
    sd_x, resp_class, coef_names, priors) {

    if (any(Y < 0) || any(is.logical(Y)) || any(Y != as.integer(Y))) {
        stop_("Response variable ", y, " is invalid: binomial family supports only non-negative integers.")
    }
    if ("factor" %in% resp_class) {
        stop_("Response variable ", y, " is invalid: binomial family is not supported for factors.")
    }
    sd_beta <- 2 / sd_x
    k <- grep("(Intercept)", coef_names)
    if (!is.null(k)) sd_beta[k] <- 2.5
    if (any(!is.finite(sd_beta))) {
        msg <- paste0("Found nonfinite prior standard deviation when using default priors for regression coeffients for response ",
            y, ", indicating constant covariate: Switching to N(0, 0.01) prior.")
        sd_beta[!is.finite(sd_beta)] <- 0.1
        warning_(msg)
    }
    prepare_channel_default(y, Y, channel,
        sd_beta, resp_class, coef_names, priors)
}

prepare_channel_bernoulli <- function(y, Y, channel,
    sd_x, resp_class, coef_names, priors) {
    if (!all(Y %in% 0:1) || any(is.logical(Y))) {
        stop_("Response variable ", y, " is invalid: bernoulli family supports only 0/1 integers.")
    }
    if ("factor" %in% resp_class) {
        stop_("Response variable ", y, " is invalid: bernoulli family is not supported for factors.")
    }
    prepare_channel_binomial(y, Y, channel,
        sd_x, resp_class, coef_names, priors)
}

prepare_channel_poisson <- function(y, Y, channel,
    sd_x, resp_class, coef_names, priors) {
    if (any(Y < 0) || any(Y != as.integer(Y))) {
        stop_("Response variable ", y, " is invalid: Poisson family supports only non-negative integers.")
    }
    if ("factor" %in% resp_class) {
        stop_("Response variable ", y, " is invalid: Poisson family is not supported for factors.")
    }
    sd_beta <- 2 / sd_x
    k <- grep("(Intercept)", coef_names)
    if (!is.null(k)) sd_beta[k] <- 10 # wider prior for intercept as we are not centering X
    if (any(!is.finite(sd_beta))) {
        msg <- paste0("Found nonfinite prior standard deviation when using default priors for regression coeffients for response ",
            y, ", indicating constant covariate: Switching to N(0, 0.01) prior.")
        sd_beta[!is.finite(sd_beta)] <- 0.1
        warning_(msg)
    }
    prepare_channel_default(y, Y, channel,
        sd_beta, resp_class, coef_names, priors)
}

prepare_channel_negbin <- function(y, Y, channel,
    sd_x, resp_class, coef_names, priors) {
    if (any(Y < 0) || any(Y != as.integer(Y))) {
        stop_("Response variable ", y, " is invalid: negative binomial family supports only non-negative integers.")
    }
    if ("factor" %in% resp_class) {
        stop_("Response variable ", y, " is invalid: negative binomial family is not supported for factors.")
    }
    sd_beta <- 2 / sd_x
    k <- grep("(Intercept)", coef_names)
    if (!is.null(k)) sd_beta[k] <- 10 # wider prior for intercept as we are not centering X
    if (any(!is.finite(sd_beta))) {
        msg <- paste0("Found nonfinite prior standard deviation when using default priors for regression coeffients for response ",
            y, ", indicating constant covariate: Switching to N(0, 0.01) prior.")
        sd_beta[!is.finite(sd_beta)] <- 0.1
        warning(msg)
    }
    out <- prepare_channel_default(y, Y, channel,
        sd_beta, resp_class, coef_names, priors)
    if (is.null(priors)) {
        out$channel$phi_prior_distr <- "exponential(1)"
        out$priors <- dplyr::bind_rows(out$priors,
            data.frame(parameter = "phi",
                response = y,
                prior = out$channel$phi_prior_distr,
                type = "phi",
                category = ""))
    } else {
        pdef <- priors |> dplyr::filter(response == y && type == "phi")
        if (nrow(pdef) == 1) {
            out$channel$phi_prior_distr <- pdef$prior
        }
    }
    out
}

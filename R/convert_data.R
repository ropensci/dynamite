# Convert data for Stan
convert_data <- function(formula, responses, specials, group, time, fixed,
                         model_matrix, coef_names, priors = NULL) {
    T_full <- length(time)
    groups <- !is.null(group)
    free_obs <- (fixed + 1):T_full
    spline_defs <- attr(formula, "splines")
    has_splines <- !is.null(spline_defs)
    if (has_splines) {
        bs_opts <- spline_defs$bs_opts
        bs_opts$x <- time[free_obs]
        if (is.null(bs_opts$Boundary.knots)) {
            bs_opts$Boundary.knots <- range(bs_opts$x)
        }
        Bs <- t(do.call(splines::bs, args = bs_opts))
        D <- nrow(Bs)
    }
    N <- T_full
    if (groups) {
        N <- length(unique(group))
    }
    K <- ncol(model_matrix)
    C <- length(get_resp(formula))
    X <- aperm(array(as.numeric(unlist(split(model_matrix, gl(T_full, 1, N * T_full)))),
                     dim = c(N, K, T_full))[,,free_obs, drop = FALSE],
               c(3, 1, 2))
    assigned <- attr(model_matrix, "assign")
    fixed_pars <- attr(model_matrix, "fixed")
    varying_pars <- attr(model_matrix, "varying")
    channel_vars <- list()
    if (N > 1) {
        sd_x <- apply(X[1, , ], 2, sd)
    } else {
        sd_x <- apply(X[, 1, ], 2, sd)
    }
    sd_x[sd_x < 0.1] <- 0.1 # Intercept and other constants at time 1
    resp_classes <- attr(responses, "resp_class")
    helpers <- list()
    warn_nosplines <- FALSE
    prior_list <- list()
    for (i in seq_along(formula)) {
        resp <- formula[[i]]$response
        form_specials <- specials[[i]]
        Ls <- c(paste0(c("L_fixed_", "L_varying_"), resp))
        Js <- c(paste0(c("J_", "J_fixed_", "J_varying_"), resp))
        Ks <- c(paste0(c("K_", "K_fixed_", "K_varying_"), resp))
        channel_vars[[Ls[1]]] <- as.array(fixed_pars[[i]])
        channel_vars[[Ls[2]]] <- as.array(varying_pars[[i]])
        channel_vars[[Js[1]]] <- as.array(assigned[[i]])
        channel_vars[[Js[2]]] <- as.array(assigned[[i]][channel_vars[[Ls[1]]]])
        channel_vars[[Js[3]]] <- as.array(assigned[[i]][channel_vars[[Ls[2]]]])
        channel_vars[[Ks[1]]] <- length(assigned[[i]])
        channel_vars[[Ks[2]]] <- length(fixed_pars[[i]])
        channel_vars[[Ks[3]]] <- length(varying_pars[[i]])
        helpers[[i]] <- list(has_fixed = channel_vars[[Ks[2]]] > 0,
                             has_varying = channel_vars[[Ks[3]]] > 0)
        if (helpers[[i]]$has_varying && !has_splines) {
            warning_("Model for response variable ", resp, " contains time-varying definitions, but splines have not been defined.")
            warn_nosplines <- TRUE
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
                    channel_vars[[paste0(spec, "_", resp)]] <- array(as.numeric(unlist(spec_split)), dim = c(T_full, N))[free_obs, , drop = FALSE]
                    helpers[[i]][[paste0("has_", spec)]] <- TRUE
                } else {
                    helpers[[i]][[paste0("has_", spec)]] <- FALSE
                }
            }
        }
        Y <- array(as.numeric(unlist(resp_split)), dim = c(T_full, N))[free_obs, , drop = FALSE]
        channel_vars[[resp]] <- Y
        prep <- do.call(paste0("prepare_channel_vars_", formula[[i]]$family),
                        list(i = resp,
                             Y = Y,
                             J_fixed = channel_vars[[Js[2]]],
                             J_varying = channel_vars[[Js[3]]],
                             L_fixed = channel_vars[[Ls[1]]],
                             L_varying = channel_vars[[Ls[2]]],
                             K_fixed = channel_vars[[Ks[2]]],
                             K_varying = channel_vars[[Ks[3]]],
                             sd_x = sd_x,
                             resp_class = resp_classes[resp],
                             coef_names = coef_names[[i]],
                             priors = priors))
        prior_list[[resp]] <- prep$priors
        channel_vars <- c(channel_vars, prep$channel_vars)
    }


    if (warn_nosplines) {
        warning_("All channels will now default to time-constant coefficients for all predictors.")
        for (i in seq_along(formula)) {
            helpers[[i]]$has_varying <- FALSE
        }
    }
    T <- T_full - fixed

    out_data <- c(named_list(T, N, C, K, X), channel_vars)
    if (has_splines) {
        out_data$D <- D
        out_data$Bs <- Bs
    }
    list(data = out_data, helpers = helpers, priors = prior_list)
}

prepare_channel_vars_default <- function(i, Y, J_fixed, J_varying, L_fixed,
                                         L_varying, K_fixed, K_varying,
                                         sd_x, resp_class, coef_names, priors) {

    channel_vars <- list()
    if (is.null(priors)) {
        priors <- list()
        #default priors
        bnames <- gsub(paste0("^", i), "beta", coef_names)
        if (K_fixed > 0) {
            m <- rep(0, K_fixed)
            s <- 2 / sd_x[J_fixed]
            channel_vars[[paste0("beta_fixed_prior_npars_", i)]] <- 2
            channel_vars[[paste0("beta_fixed_prior_pars_", i)]] <- cbind(m, s, deparse.level = 0)
            channel_vars[[paste0("beta_fixed_prior_distr_", i)]] <- "normal"

            priors$beta_fixed <-
                data.frame(parameter = bnames[L_fixed],
                           response = i,
                           prior = paste0("normal(", m, ", ", s, ")"),
                           type = "beta_fixed")
        }

        if (K_varying > 0) {
            m <- rep(0, K_varying)
            s <- 2 / sd_x[J_varying]
            channel_vars[[paste0("beta_varying_prior_npars_", i)]] <- 2
            channel_vars[[paste0("beta_varying_prior_pars_", i)]] <- cbind(m, s, deparse.level = 0)
            channel_vars[[paste0("beta_varying_prior_distr_", i)]] <- "normal"

            priors$beta_varying <-
                data.frame(parameter = bnames[L_varying],
                           response = i,
                           prior = paste0("normal(", m, ", ", s, ")"),
                           type = "beta_varying")

            channel_vars[[paste0("tau_prior_npars_", i)]] <- 2
            channel_vars[[paste0("tau_prior_pars_", i)]] <- cbind(0, rep(1, K_varying))
            channel_vars[[paste0("tau_prior_distr_", i)]] <- "normal"
            priors$tau <-
                data.frame(parameter = paste0("tau", bnames[L_varying]),
                           response = i,
                           prior = "normal(0, 1)",
                           type = "tau")
        }
        priors <- dplyr::bind_rows(priors)
    } else {
        # TODO add a warning to documentation that the only the 'prior' column
        # of the priors data.frame should be altered (i.e. there's no checks for names or reordering of rows)
        # Or arrange...
        priors <- priors |> dplyr::filter(response == i)
        for (ptype in c("beta_fixed", "beta_varying", "tau")) {
            pdef <- priors |> dplyr::filter(type == ptype)
            if (nrow(pdef) > 0) {
                dists <- sub("\\(.*", "", pdef$prior)
                vectorized <- length(unique(dists)) == 1
                if (vectorized) {
                    pars <- strsplit(sub(".*\\((.*)\\).*", "\\1", pdef$prior), ",")
                    pars <- do.call("rbind", lapply(pars, as.numeric))
                    channel_vars[[paste0(ptype, "_prior_npars_", i)]] <- ncol(pars)
                    channel_vars[[paste0(ptype, "_prior_pars_", i)]] <- pars
                    channel_vars[[paste0(ptype, "_prior_distr_", i)]] <- dists[1]
                } else {
                    channel_vars[[paste0(ptype, "_prior_distr_", i)]] <- pdef$prior # write separate priors
                }
            }
        }
    }

    list(channel_vars = channel_vars, priors = priors)
}

prepare_channel_vars_categorical <- function(i, Y, J_fixed, J_varying, L_fixed,
                                             L_varying, K_fixed, K_varying,
                                             sd_x, resp_class, coef_names, priors) {
    S_i <- length(unique(as.vector(Y)))
    channel_vars <- list()
    channel_vars[[paste0("S_", i)]] <- S_i
    if (is.null(priors)) {
        priors <- list()
        #default priors
        bnames <- gsub(paste0("^", i), "beta", attr(coef_names, "simplified")$names)
        levels_ <- attr(coef_names, "simplified")$levels
        if (K_fixed > 0) {
            m <- rep(0, K_fixed * (S_i - 1))
            s <- rep(2 / sd_x[J_fixed], S_i - 1)
            channel_vars[[paste0("beta_fixed_prior_npars_", i)]] <- 2
            channel_vars[[paste0("beta_fixed_prior_pars_", i)]] <- cbind(m, s, deparse.level = 0)
            channel_vars[[paste0("beta_fixed_prior_distr_", i)]] <- "normal"

            priors$beta_fixed <-
                data.frame(parameter = bnames[L_fixed],
                           response = paste0(i, "_", rep(levels_, each = K_fixed)),
                           prior = paste0("normal(", m, ", ", s, ")"),
                           type = "beta_fixed")
        }

        if (K_varying > 0) {
            m <- rep(0, K_varying * (S_i - 1))
            s <- rep(2 / sd_x[J_varying], S_i - 1)
            channel_vars[[paste0("beta_varying_prior_npars_", i)]] <- 2
            channel_vars[[paste0("beta_varying_prior_pars_", i)]] <- cbind(m, s, deparse.level = 0)
            channel_vars[[paste0("beta_varying_prior_distr_", i)]] <- "normal"

            priors$beta_varying <-
                data.frame(parameter = bnames[L_varying],
                           response = paste0(i, "_", rep(levels_, each = K_varying)),
                           prior = paste0("normal(", m, ", ", s, ")"),
                           type = "beta_varying")

            channel_vars[[paste0("tau_prior_npars_", i)]] <- 2
            channel_vars[[paste0("tau_prior_pars_", i)]] <- cbind(0, rep(1, K_varying))
            channel_vars[[paste0("tau_prior_distr_", i)]] <- "normal"
            priors$tau <-
                data.frame(parameter = paste0("tau", bnames[L_varying]),
                           response = i,
                           prior = "normal(0, 1)",
                           type = "tau")
        }
        priors <- dplyr::bind_rows(priors)
    } else {
        # TODO add a warning to documentation that the only the 'prior' column
        # of the priors data.frame should be altered (i.e. there's no checks for names or reordering of rows)
        # Or arrange...
        priors <- priors |> dplyr::filter(response == i)
        for (ptype in c("beta_fixed", "beta_varying", "tau")) {
            pdef <- priors |> dplyr::filter(type == ptype)
            if (nrow(pdef) > 0) {
                dists <- sub("\\(.*", "", pdef$prior)
                vectorized <- length(unique(dists)) == 1
                if (vectorized) {
                    pars <- strsplit(sub(".*\\((.*)\\).*", "\\1", pdef$prior), ",")
                    pars <- do.call("rbind", lapply(pars, as.numeric))
                    channel_vars[[paste0(ptype, "_prior_npars_", i)]] <- ncol(pars)
                    channel_vars[[paste0(ptype, "_prior_pars_", i)]] <- pars
                    channel_vars[[paste0(ptype, "_prior_distr_", i)]] <- dists[1]
                } else {
                    channel_vars[[paste0(ptype, "_prior_distr_", i)]] <- pdef$prior # write separate priors
                }
            }
        }
    }

    list(channel_vars = channel_vars, priors = priors)
}

prepare_channel_vars_gaussian <- function(i, Y, J_fixed, J_varying, L_fixed,
                                          L_varying, K_fixed, K_varying,
                                          sd_x, resp_class, coef_names, priors) {
    if ("factor" %in% resp_class) {
        stop_("Response variable ", i, " is invalid: gaussian family is not supported for factors.")
    }
    out <- prepare_channel_vars_default(i, Y, J_fixed, J_varying, L_fixed,
                                        L_varying, K_fixed, K_varying,
                                        sd_x, resp_class, coef_names, priors)
    if (is.null(priors)) {
        s <- 1 / mean(apply(Y, 1, sd))
        out$channel_vars[[paste0("sigma_prior_distr_", i)]] <- paste0("exponential(", s, ")")
        out$priors <- dplyr::bind_rows(out$priors,
                                       data.frame(parameter = "sigma",
                                                  response = i,
                                                  prior = paste0("exponential(", s, ")"),
                                                  type = "sigma"))
    } else {
        pdef <- priors |> dplyr::filter(response == i && type == "sigma")
        if (nrow(pdef) == 1) {
            out$channel_vars[[paste0("sigma_prior_distr_", i)]] <- pdef$prior
        }
    }
    out
}

prepare_channel_vars_binomial <- function(i, Y, J_fixed, J_varying, L_fixed,
                                          L_varying, K_fixed, K_varying,
                                          sd_x, resp_class, coef_names, priors) {
    if (any(Y < 0) || any(Y != as.integer(Y))) {
        stop_("Response variable ", i, " is invalid: binomial family supports only non-negative integers.")
    }
    if ("factor" %in% resp_class) {
        stop_("Response variable ", i, " is invalid: binomial family is not supported for factors.")
    }
    prepare_channel_vars_default(i, Y, J_fixed, J_varying, L_fixed,
                                 L_varying, K_fixed, K_varying,
                                 sd_x, resp_class, coef_names, priors)
}

prepare_channel_vars_bernoulli <- function(i, Y, J_fixed, J_varying, L_fixed,
                                           L_varying, K_fixed, K_varying,
                                           sd_x, resp_class, coef_names, priors) {
    if (!all(Y %in% 0:1)) {
        stop_("Response variable ", i, " is invalid: bernoulli family supports only 0/1 integers.")
    }
    if ("factor" %in% resp_class) {
        stop_("Response variable ", i, " is invalid: bernoulli family is not supported for factors.")
    }
    prepare_channel_vars_binomial(i, Y, J_fixed, J_varying, L_fixed,
                                  L_varying, K_fixed, K_varying,
                                  sd_x, resp_class, coef_names, priors)
}

prepare_channel_vars_poisson <- function(i, Y, J_fixed, J_varying, L_fixed,
                                         L_varying, K_fixed, K_varying,
                                         sd_x, resp_class, coef_names, priors) {
    if (any(Y < 0) || any(Y != as.integer(Y))) {
        stop_("Response variable ", i, " is invalid: Poisson family supports only non-negative integers.")
    }
    if ("factor" %in% resp_class) {
        stop_("Response variable ", i, " is invalid: Poisson family is not supported for factors.")
    }
    prepare_channel_vars_binomial(i, Y, J_fixed, J_varying, L_fixed,
                                  L_varying, K_fixed, K_varying,
                                  sd_x, resp_class, coef_names, priors)
}

prepare_channel_vars_negbin <- function(i, Y, J_fixed, J_varying, L_fixed,
                                        L_varying, K_fixed, K_varying,
                                        sd_x, resp_class, coef_names, priors) {
    if (any(Y < 0) || any(Y != as.integer(Y))) {
        stop_("Response variable ", i, " is invalid: negative binomial family supports only non-negative integers.")
    }
    if ("factor" %in% resp_class) {
        stop_("Response variable ", i, " is invalid: negative binomial family is not supported for factors.")
    }
    out <- prepare_channel_vars_default(i, Y, J_fixed, J_varying, L_fixed,
                                        L_varying, K_fixed, K_varying,
                                        sd_x, resp_class, coef_names, priors)
    if (is.null(priors)) {
        out$channel_vars[[paste0("phi_prior_distr_", i)]] <- "exponential(1)"
        out$priors <- dplyr::bind_rows(out$priors,
                                       data.frame(parameter = "phi",
                                                  response = i,
                                                  prior = paste0("exponential(1)"),
                                                  type = "phi"))
    } else {
        pdef <- priors |> dplyr::filter(response == i && type == "phi")
        if (nrow(pdef) == 1) {
            out$channel_vars[[paste0("phi_prior_distr_", i)]] <- pdef$prior
        }
    }
    out
}

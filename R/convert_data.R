# Convert data for Stan
convert_data <- function(formula, responses, specials, group, time, fixed, model_matrix) {
    T_full <- 0
    groups <- !is.null(group)
    if (groups) {
        id_tab <- table(group)
        if (!all(id_tab == id_tab[1])) {
            # TODO In theory we could allow unequal number of observations
            # Just need to define time from smallest t to the largest,
            # and use proper indexing.
            stop_("Unequal number of time points")
        }
        T_full <- as.integer(id_tab[1])
    } else {
        T_full <- nrow(responses)
    }
    if (is.null(time)) {
        time <- 1:T_full
    }
    free_obs <- (fixed + 1):T_full
    bs_opts <- attr(formula, "splines")$bs_opts
    bs_opts$x <- time[free_obs]
    if (is.null(bs_opts$Boundary.knots)) {
        bs_opts$Boundary.knots <- range(bs_opts$x)
    }
    Bs <- t(do.call(splines::bs, args = bs_opts))
    D <- nrow(Bs)
    N <- T_full
    if (groups) {
        N <- length(id_tab)
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
    sd_x[sd_x < 1] <- 1 # Intercept and other constants at time 1
    resp_classes <- attr(responses, "resp_class")
    helpers <- list()
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
                             K_fixed = channel_vars[[Ks[2]]],
                             K_varying = channel_vars[[Ks[3]]],
                             sd_x = sd_x,
                             resp_class = resp_classes[resp]))
        channel_vars <- c(channel_vars, prep)
    }
    T <- T_full - fixed
    list(data = c(named_list(T, N, C, K, X, D, Bs), channel_vars),
         helpers = helpers)
}

prepare_channel_vars_categorical <- function(i, Y, J_fixed, J_varying, K_fixed, K_varying, sd_x, resp_class) {
    S_i <- length(unique(as.vector(Y)))
    channel_vars <- list()
    channel_vars[[paste0("S_", i)]] <- S_i
    channel_vars[[paste0("beta_prior_mean_", i)]] <- matrix(0, K_fixed, S_i - 1)
    channel_vars[[paste0("beta_prior_sd_", i)]] <- matrix(5, K_fixed, S_i - 1) # TODO better initial values
    channel_vars[[paste0("a_prior_mean_", i)]] <- matrix(0, K_varying, S_i - 1)
    channel_vars[[paste0("a_prior_sd_", i)]] <- matrix(2 / sd_x[J_varying], K_varying, S_i - 1)
    channel_vars
}

prepare_channel_vars_gaussian <- function(i, Y, J_fixed, J_varying, K_fixed, K_varying, sd_x, resp_class) {
    if ("factor" %in% resp_class) {
        stop_("Response variable ", resp, " is invalid: gaussian family is not supported for factors.")
    }
    channel_vars <- list()
    channel_vars[[paste0("beta_prior_mean_", i)]] <- array(0, K_fixed)
    channel_vars[[paste0("beta_prior_sd_", i)]] <- array(5, K_fixed) # TODO better initial values
    channel_vars[[paste0("a_prior_mean_", i)]] <- array(0, K_varying)
    # TODO adjust prior mean for the intercept term under the assumption that other betas/x are 0
    channel_vars[[paste0("a_prior_sd_", i)]] <- array(2 * sd(Y[1, ]) / sd_x[J_varying], K_varying)
    channel_vars[[paste0("sigma_scale_", i)]] <- 1 / mean(apply(Y, 1, sd))
    channel_vars
}

prepare_channel_vars_binomial <- function(i, Y, J_fixed, J_varying, K_fixed, K_varying, sd_x, resp_class) {
    channel_vars <- list()
    channel_vars[[paste0("beta_prior_mean_", i)]] <- array(0, K_fixed)
    channel_vars[[paste0("beta_prior_sd_", i)]] <- array(5, K_fixed) # TODO better initial values
    channel_vars[[paste0("a_prior_mean_", i)]] <- array(0, K_varying)
    channel_vars[[paste0("a_prior_sd_", i)]] <- array(2 / sd_x[J_varying], K_varying)
    channel_vars
}

prepare_channel_vars_bernoulli <- function(...) {
    prepare_channel_vars_binomial(...)
}

prepare_channel_vars_poisson <- function(i, resp_class, ...) {
    if ("factor" %in% resp_class) {
        stop_("Response variable ", i, " is invalid: Poisson family is not supported for factors.")
    }
    prepare_channel_vars_binomial(i = i, resp_class = resp_class, ...)
}

prepare_channel_vars_negbin <- function(i, Y, J_fixed, J_varying, K_fixed, K_varying, sd_x, resp_class) {
    if ("factor" %in% resp_class) {
        stop_("Response variable ", resp, " is invalid: negative binomial family is not supported for factors.")
    }
    channel_vars <- list()
    channel_vars[[paste0("beta_prior_mean_", i)]] <- array(0, K_fixed)
    channel_vars[[paste0("beta_prior_sd_", i)]] <- array(5, K_fixed) # TODO better initial values
    channel_vars[[paste0("a_prior_mean_", i)]] <- array(0, K_varying)
    channel_vars[[paste0("a_prior_sd_", i)]] <- array(2 / sd_x[J_varying], K_varying)
    channel_vars[[paste0("phi_scale_", i)]] <- 1
    channel_vars
}

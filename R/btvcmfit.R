#' @include utilities.R

# btvcmfit class
# TODO add options to return just the created model code or data without compiling & sampling for debugging etc
btvcmfit <- function(formula, data, group, time, ...) {
    dots <- list(...)
    if (missing(group)) {
        group <- NULL
    } else {
        # TODO allow group = c(ID1, ID2, ...) etc.
        group_var <- deparse(substitute(group))
        group <- data[[group_var]] # had , drop = FALSE, is it needed somewhere?
        if (is.factor(group)) group <- droplevels(group)
    }
    # TODO is there a better way?
    if (missing(time))
        stop("Argument 'time' is missing. ")
    time_var <- deparse(substitute(time))
    data <- dplyr::arrange(data, dplyr::across(dplyr::all_of(c(group_var, time_var))))
    time <- unique(data[[time_var]])
    resp_all <- get_resp(formula)
    n_rows <- nrow(data)
    n_resp <- length(resp_all)
    lag_map <- extract_lags(unlist(get_pred(formula)))
    data_names <- names(data)
    fixed <- 0L
    # Process lag terms defined via lags()
    if (!is.null(lag_all <- attr(formula, "lags"))) {
        fixed <- lag_all$k
        pred_lag <- character(fixed * n_resp)
        for (i in seq_len(fixed)) {
            for (j in seq_len(n_resp)) {
                ix <- (i-1)*n_resp + j
                pred_lag[ix] <- paste0("I(lag_(", resp_all[j], ", ", i, "))")
            }
        }
        for (j in seq_len(n_resp)) {
            c(formula[[j]]$predictors) <- pred_lag
            icpt <- attr(terms(formula[[j]]$formula), "intercept")
            formula[[j]]$formula <- reformulate(formula[[j]]$predictors, intercept = icpt)
        }
        lag_map <- lag_map[lag_map$def %in% resp_all & lag_map$k <= fixed,]
    }
    # Process lag terms defined via lag() in formulas
    if (nrow(lag_map)) {
        fixed <- max(max(lag_map$k), fixed)
        for (i in seq_along(lag_map)) {
            if (lag_map$k[i] <= 0) {
                stop_("Only positive shift values are allowed in lag()")
            }
            if (is_as_is(lag_map$src[i])) {
                # Swap to internal lag function
                lag_map$scr[i] <- gsub("lag", "lag_", lag_map$src[i])
            }
            if (is_as_is(lag_map$def[i])) {
                # Remove 'as is' from inside of lag terms
                lag_map$def[i] <- gsub("I\\((.*)\\)", "\\1", lag_map$def[i])
            }
            pred_lag <- paste0("I(lag_(", lag_map$def[i], ", ", lag_map$k[i], "))")
            for (j in seq_len(n_resp)) {
                formula[[j]]$predictors <- gsub(lag_map$src[i], pred_lag, formula[[j]]$predictors, fixed = TRUE)
                icpt <- attr(terms(formula[[j]]$formula), "intercept")
            }
        }
    }
    #    model_matrices <- list()
    #    for (i in seq_len(n_resp)) {
    #        model_matrices[[i]] <- model.matrix(formula[[i]]$formula, data)
    #        if (any(ind <- formula[[i]]$predictors %in% names(lag_resp_map[[i]]))) {
    #            ind <- which(ind)
    #            assign_i <- attr(model_matrices[[i]], "assign")
    #            cols_i <- colnames(model_matrices[[i]])
    #            for (j in seq_along(ind)) {
    #                pred <- formula[[i]]$predictors[ind[j]]
    #                lag_cols <- which(assign_i == ind[j])
    #                lag_resp_map[[i]][[pred]]$cols <- cols_i[lag_cols]
    #            }
    #        }
    #    }
    #    model_matrix <- do.call(cbind, model_matrices)
    #    u_names <- unique(colnames(model_matrix))
    #    # TODO: Is intercept always named as (Intercept)?
    #    # checking assign attributes for 0 is an another option
    #    if (any(ind <- u_names == "(Intercept)")) {
    #        u_names <- u_names[-which(ind)]
    #    }
    #    model_matrix <- model_matrix[, u_names]
    #    assigned <- lapply(model_matrices, function(x) {
    #        which(u_names %in% colnames(x))
    #    })
    #    for (i in seq_len(n_resp)) {
    #        for (j in seq_along(lag_resp_map[[i]])) {
    #            lag_resp_map[[i]][[j]]$cols <- which(u_names %in% lag_resp_map[[i]][[j]]$cols)
    #        }
    # }
    responses <- data[, resp_all, drop = FALSE]
    model_matrix <- full_model.matrix(formula, data, resp_all)
    model_data <- convert_data(formula, responses, group, time, fixed, model_matrix)
    model_code <- create_blocks(formula, indent = 2L, model_data)
    debug <- dots$debug
    model <- if (isTRUE(debug$no_compile)) {
        NULL
    } else {
        message("Compiling Stan model")
        rstan::stan_model(model_code = model_code)
    }
    stanfit <- if (isTRUE(debug$no_compile) || isTRUE(debug$no_sampling)) NULL else rstan::sampling(model, data = model_data, ...)
    # TODO return the function call for potential update method?
    out <- structure(
        list(
            stanfit = stanfit,
            # TODO could just extract these from prediction_basis$past
            coef_names = attr(model_matrix, "coef_names"),
            # TODO what else do we need to return?
            time = time,
            time_var = time_var,
            group_var = group_var,
            levels = lapply(droplevels(responses), levels),
            # TODO: extract only D for as.data.frame and J for predict
            #model_data = model_data,
            data = data,
            # spline = list(B = model_data$Bs,
            #     D = model_data$D),
            prediction_basis = list(
                formula = formula,
                fixed = fixed,
                past = model_matrix[(n_rows - fixed):n_rows,],
                start = model_matrix[1:fixed,], # Needed for some posterior predictive checks?
                ord = data_names[!data_names %in% c(group_var, time_var)],
                J = attr(model_matrix, "assign")
            )
        ),
        class = "btvcmfit"
    )
    for (opt in names(debug)) {
        if (debug[[opt]]) {
            got <- try(get(x = opt), silent = TRUE)
            if (!"try-error" %in% class(got)) {
                out[[opt]] <- got
            }
        }
    }
    out
}

# Combine model.matrix objects of all formulas of a btvcmformula into one
full_model.matrix <- function(formula, data, resp_all) {
    model_matrices <- lapply(lapply(formula, "[[", "formula"), model.matrix, data)
    model_matrix <- do.call(cbind, model_matrices)
    u_names <- unique(colnames(model_matrix))
    model_matrix <- model_matrix[, u_names, drop = FALSE]
    attr(model_matrix, "assign") <- lapply(model_matrices, function(x) {
        which(u_names %in% colnames(x))
    })
    # TODO: simplify I(lag(variable, 1)) to something shorter, e.g. lag_1(variable)?
    # TODO: shorten variable and or channel name if they are very long?
    attr(model_matrix, "coef_names") <- lapply(seq_along(resp_all), function(i) {
        paste0(resp_all[i], "_", u_names[attr(model_matrix, "assign")[[i]]])
    })
    model_matrix
}
# For prediction
full_model.matrix_fast <- function(formula, data, u_names) {
    model_matrices <- lapply(lapply(formula, "[[", "formula"), model.matrix, data)
    model_matrix <- do.call(cbind, model_matrices)
    model_matrix[, u_names, drop = FALSE]
}

# Convert data for Stan
convert_data <- function(formula, responses, group, time, fixed, model_matrix) {
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
    channel_vars <- list()
    if (N > 1) {
        sd_x <- apply(X[1, , ], 2, sd)
    } else {
        sd_x <- apply(X[, 1, ], 2, sd)
    }
    sd_x[sd_x < 1] <- 1 # Intercept and other constants at time 1
    for (i in seq_along(formula)) {
        channel_vars[[paste0("J_", i)]] <- assigned[[i]]
        channel_vars[[paste0("K_", i)]] <- length(assigned[[i]])
        if (groups) {
            resp_split <- split(responses[,formula[[i]]$response], group)
        } else {
            resp_split <- responses[,formula[[i]]$response]
        }
        Y <- array(as.numeric(unlist(resp_split)), dim = c(T_full, N))[free_obs, , drop = FALSE]
        channel_vars[[paste0("response_", i)]] <- Y
        if (is_categorical(formula[[i]]$family)) {
            S_i <- length(unique(as.vector(Y)))
            K_i <- length(assigned[[i]])
            prior_sds <- matrix(2 / sd_x[assigned[[i]]], K_i, S_i - 1)
            channel_vars[[paste0("S_", i)]] <- S_i
            channel_vars[[paste0("a_prior_mean_", i)]] <- matrix(0, K_i, S_i - 1)
            channel_vars[[paste0("a_prior_sd_", i)]] <- prior_sds
        }
        if (is_gaussian(formula[[i]]$family)) {
            channel_vars[[paste0("a_prior_mean_", i)]] <- rep(0, length(assigned[[i]]))
            # TODO adjust prior mean for the intercept term under the assumption that other betas/x are 0
            #if(model_matrix[])
            channel_vars[[paste0("a_prior_sd_", i)]] <- 2 * sd(Y[1, ]) / sd_x[assigned[[i]]]
            channel_vars[[paste0("sigma_scale_", i)]] <- 1 / mean(apply(Y, 1, sd))
        }
        # TODO switch case other families, with Gamma shape, negbin dispersion
        # TODO or something like do.call(paste0("prepare_channel_vars_", formula[[i]]$family, list(Y, S_i, K_i, i ...))) for branchless

    }
    T <- T_full - fixed
    c(named_list(T, N, C, K, X, D, Bs), channel_vars)
}

#' @include utilities.R

# btvcmfit class
# TODO add options to return just the created model code or data without compiling & sampling for debugging etc
btvcmfit <- function(formula, data, group, time, ...) {
    dots <- list(...)
    if (missing(group)) {
        group <- NULL
    } else {
        # TODO allow group = c(ID1, ID2, ...) etc.
        group <- data[,deparse(substitute(group)),drop = FALSE]
    }
    if (missing(time)) {
        time <- NULL
    } else {
        # TODO is there a better way?
        time <- sort(unique(data[,deparse(substitute(time))]))
    }
    resp_all <- get_resp(formula)
    pred_all <- unlist(get_pred(formula))
    lag_map <- extract_lags(pred_all)
    n_rows <- nrow(data)
    n_resp <- length(resp_all)
    data_vars <- names(data)
    fixed <- 0L
    # Process lag terms defined via lags()
    if (!is.null(lag_all <- attr(formula, "lags"))) {
        fixed <- lag_all$k
        pred_lag <- character(fixed * n_resp)
        for (i in 1:fixed) {
            for (j in 1:n_resp) {
                ix <- (i-1)*n_resp + j
                pred_lag[ix] <- paste0("I(lag_(", resp_all[j], ", ", i, "))")
            }
        }
        for (j in 1:n_resp) {
            c(formula[[j]]$predictors) <- pred_lag
            icpt <- attr(terms(formula[[j]]$formula), "intercept")
            formula[[j]]$formula <- reformulate(formula[[j]]$predictors, intercept = icpt)
        }
        #c(pred_all) <- pred_lag
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
            for (j in 1:n_resp) {
                formula[[j]]$predictors <- gsub(lag_map$src[i], pred_lag, formula[[j]]$predictors, fixed = TRUE)
                icpt <- attr(terms(formula[[j]]$formula), "intercept")
                formula[[j]]$formula <- reformulate(formula[[j]]$predictors, intercept = icpt)
            }
            #pred_all <- gsub(lag_map$src[i], pred_lag, pred_all, fixed = TRUE)
        }
    }
    #all_rhs_vars <- unique(pred_all)
    # resp_all <- get_resp(formula)
    model_matrix <- NULL
    assigned <- list()
    if (n_resp > 1) {
        model_matrices <- lapply(lapply(formula, "[[", "formula"), model.matrix, data)
        model_matrix <- do.call(cbind, model_matrices)
        u_names <- unique(colnames(model_matrix))
        model_matrix <- model_matrix[, u_names, drop = FALSE]
        assigned <- lapply(model_matrices, function(x) {
            which(u_names %in% colnames(x))
        })

    } else {
        model_matrix <- model.matrix(formula[[1]]$formula, data)
        u_names <- colnames(model_matrix)
        assigned <- list(1:ncol(model_matrix))
    }
    # TODO: simplify I(lag(variable, 1)) to something shorter eg lag_1(variable)?
    coef_names <- lapply(1:n_resp, function(i) {
        paste0(resp_all[i], "_", u_names[assigned[[i]]])
    })
    # Place lags last
    # all_lags <- find_lags(all_rhs_vars, processed = TRUE)
    # all_rhs_vars <- c(all_rhs_vars[all_lags], all_rhs_vars[!all_lags])
    # all_rhs_formula <- reformulate(all_rhs_vars, intercept = FALSE)
    # model_matrix <- model.matrix(all_rhs_formula, data)
    responses <- data[, resp_all, drop = FALSE]
    model_data <- convert_data(formula, responses, group, time, model_matrix, assigned, fixed)
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
            coef_names = coef_names,
            # TODO what else do we need to return?
            time = if(is.null(time)) 1:model_data$T else time[-(1:fixed)],
            model_data = model_data # TODO: remove X and responses?
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

# Convert data for Stan
convert_data <- function(formula, responses, group, time, model_matrix, assigned, fixed) {
    T_full <- 0
    groups <- !is.null(group)
    if (groups) {
        id_tab <- table(group)
        if (!all(id_tab == id_tab[1])) {
            stop_("Unequal number of time points")
        }
        T_full <- as.integer(id_tab[1])
    } else {
        T_full <- nrow(responses)
    }
    if (is.null(time)) {
        time <- 1:T_full
    }
    # TODO take spline definition into account here
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
    # assigned <- attr(model_matrix, "assign")
    channel_vars <- list()
    sd_x <- apply(X[1, , ], 2, sd)
    sd_x[sd_x == 0] <- 1 # Intercept and other constants at time 1
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
            prior_sds <- matrix(2 / sd_x[assigned[[i]]], K_i, S_i)
            channel_vars[[paste0("S_", i)]] <- S_i
            channel_vars[[paste0("a_prior_mean_", i)]] <- matrix(0, K_i, S_i)
            channel_vars[[paste0("a_prior_sd_", i)]] <- prior_sds
        }
        if (is_gaussian(formula[[i]]$family)) {
            sd_y <- sd(Y)
            channel_vars[[paste0("a_prior_mean_", i)]] <- rep(0, length(assigned[[i]]))
            # TODO adjust prior mean for the intercept term under the assumption that other betas/x are 0
            #if(model_matrix[])
            channel_vars[[paste0("a_prior_sd_", i)]] <- 2 * sd_y / sd_x[assigned[[i]]]
            channel_vars[[paste0("sigma_scale_", i)]] <- 1 / sd_y
        }
        # TODO switch case other families, with Gamma shape, negbin dispersion

    }
    T <- T_full - fixed
    c(named_list(T, N, C, K, X, D, Bs), channel_vars)
}

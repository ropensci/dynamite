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
        }
        c(pred_all) <- pred_lag
        lag_map <- lag_map[lag_map$var %in% resp_all & lag_map$k <= fixed,]
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
            }
            pred_all <- gsub(lag_map$src[i], pred_lag, pred_all, fixed = TRUE)
        }
    }
    all_rhs_vars <- unique(pred_all)
    # Place lags last
    all_lags <- find_lags(all_rhs_vars, processed = TRUE)
    all_rhs_vars <- c(all_rhs_vars[all_lags], all_rhs_vars[!all_lags])
    all_rhs_formula <- reformulate(all_rhs_vars, intercept = FALSE)
    model_matrix <- model.matrix(all_rhs_formula, data)
    responses <- data[,resp_all,drop = FALSE]
    model_data <- convert_data(formula, responses, group, time, model_matrix, all_rhs_vars, fixed)
    model_code <- create_blocks(formula, indent = 2L)
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
            stanfit = stanfit
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
convert_data <- function(formula, responses, group, time, model_matrix, all_rhs_vars, fixed) {
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
    # knots <- seq(time[1 + fixed], timer[T_full], length.out = min(10, T_full - fixed))
    # knots <- knots[2:(length(knots)-1)]
    Bs <- do.call(splines::bs, args = bs_opts)
    # Use sum-to-zero constraint so that we can separate mean beta and spline effect
    # based on Wood (2006) section 5.4.1 (QR version)
    Z <- qr.Q(qr(colSums(Bs)), complete = TRUE)[, 2:ncol(Bs)]
    Bs <- t(Bs %*% Z)
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
    for (i in seq_along(formula)) {
        pred_pos <- which(all_rhs_vars %in% formula[[i]]$predictors)
        pred_ix <- which(assigned %in% pred_pos)
        channel_vars[[paste0("J_", i)]] <- pred_ix
        channel_vars[[paste0("K_", i)]] <- length(pred_ix)
        if (groups) {
            resp_split <- split(responses[,formula[[i]]$response], group)
        } else {
            resp_split <- responses[,formula[[i]]$response]
        }
        Y <- array(as.numeric(unlist(resp_split)), dim = c(T_full, N))[free_obs, , drop = FALSE]
        channel_vars[[paste0("response_", i)]] <- Y
        if (is_categorical(formula[[i]]$family)) {
            channel_vars[[paste0("S_", i)]] <- length(unique(as.vector(Y)))
        }
    }
    T <- T_full - fixed
    c(named_list(T, N, C, K, X, D, Bs), channel_vars)
}

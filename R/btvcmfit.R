#' @include utilities.R

# btvcmfit class
# TODO add options to return just the created model code or data without compiling & sampling for debugging etc
btvcmfit <- function(formula, data, group, time, ...) {
    if (missing(group)) {
        # TODO do we need to allow this?
        # Does it make sense to have a single individual
        # Not a priority but maybe reasonable special case, and allows comparison with the walker package
        stop_("No groups specified")
    }
    # TODO allow group = c(ID1, ID2, ...) etc.
    group <- deparse(substitute(group))
    resp_all <- get_resp(formula)
    pred_all <- unlist(get_pred(formula))
    fixed <- 0L
    if (!is.null(lag_all <- attr(formula, "lags"))) {
        # TODO remove all manual lag terms from formula
        fixed <- lag_all$k
        n_rows <- nrow(data)
        n_resp <- length(resp_all)
        pred_lag <- character(fixed * n_resp)
        for (i in 1:fixed) {
            remove <- n_rows:(n_rows-i+1)
            for (j in 1:n_resp) {
                ix <- (i-1)*n_resp + j
                pred_lag[ix] <- paste0(resp_all[j], "_lag", i)
                data[,pred_lag[ix]] <- c(rep(0L, i), data[-remove,resp_all[j]])
            }
        }
        # TODO add check for name conflicts, and replace conflicting names
        for (j in 1:n_resp) {
            c(formula[[j]]$predictors) <- pred_lag
        }
        c(pred_all) <- pred_lag
    } else {
        # TODO process manual lag terms, if any
    }
    all_rhs_vars <- unique(pred_all)
    all_rhs_formula <- reformulate(all_rhs_vars, intercept = FALSE)
    moma <- model.matrix(all_rhs_formula, data)
    model_data <- convert_data(formula, data, group, time, moma, all_rhs_vars, fixed)
    model_code <- create_blocks(formula)
    message("Compiling Stan model")
    model <- rstan::stan_model(model_code = model_code)
    stanfit <- rstan::sampling(model, data = model_data, ...)
    # TODO return the function call for potential update method?
    structure(
        list(stanfit = stanfit),
        class = "btvcmfit"
    )
}

# Convert data for Stan
# TODO might not need to copy data here, investigate
convert_data <- function(formula, data, group, time, moma, all_rhs_vars, fixed) {
    id_tab <- table(data[,group])
    if (!all(id_tab == id_tab[1])) {
        stop_("Unequal number of time points")
    }
    T_full <- as.integer(id_tab[1])
    if (missing(time)) {
        time_var <- 1:T_full
    } else {
        # TODO is there a better way?
        time_var <- sort(unique(data[[deparse(substitute(time))]]))
    }
    # TODO take spline definition into account here
    # TODO take fixed time points into account in spline
    free_obs <- (fixed+1):T_full
    knots <- seq(time_var[2 + fixed], time_var[T_full-1], length.out = min(10, T_full - 2 - fixed))
    Bs <- t(splines::bs(time_var[free_obs], knots = knots, degree = 3, intercept = TRUE))
    D <- nrow(Bs)
    N <- length(id_tab)
    K <- ncol(moma)
    C <- length(get_resp(formula))
    moma_split <- split(moma, gl(T, 1, N * T))
    X <- aperm(array(as.numeric(unlist(moma_split)), dim = c(N, K, T_full))[,,free_obs], c(3, 1, 2))
    assigned <- attr(moma, "assign")
    channel_vars <- list()
    for (i in seq_along(formula)) {
        pred_pos <- which(all_rhs_vars %in% formula[[i]]$predictors)
        pred_ix <- which(assigned %in% pred_pos)
        channel_vars[[paste0("J_", i)]] <- pred_ix
        channel_vars[[paste0("K_", i)]] <- length(pred_ix)
        resp_split <- split(data[,formula[[i]]$response], group)
        Y <- array(as.numeric(unlist(resp_split)), dim = c(T_full, N))[free_obs, ]
        channel_vars[[paste0("response_", i)]] <- Y
        if (is_categorical(formula[[i]]$family)) {
            channel_vars[[paste0("S_", i)]] <- length(unique(as.vector(Y)))
        }
    }
    T <- T_full - fixed
    c(named_list(T, N, C, K, X, D, Bs), channel_vars)
}

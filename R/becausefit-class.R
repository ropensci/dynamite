#' @include utilities.R

# becausefit class
becausefit <- function(formula, data, group, time, ...) {
    if (missing(group)) {
        # TODO do we need to allow this?
        # Does it make sense to have a single individual
        stop_("No groups specified")
    }
    # TODO allow group = c(ID1, ID2, ...) etc.
    group <- deparse(substitute(group))
    resp_all <- get_resp(formula)
    pred_all <- get_pred(formula)
    all_rhs_vars <- unique(c(unlist(pred_all), resp_all))
    all_rhs_formula <- reformulate(all_rhs_vars, intercept = FALSE)
    moma <- model.matrix(all_rhs_formula, data)
    model_data <- convert_data(formula, data, group, time, moma, all_rhs_vars)
    model_code <- create_blocks(formula)
    message("Compiling Stan model")
    model <- rstan::stan_model(model_code = model_code)
    out <- rstan::sampling(model, data = model_data, ...)
    class(out) <- "becausefit"
    out
}

# Convert data for Stan
# TODO might not need to copy data here, investigate
convert_data <- function(formula, data, group, time, moma, all_rhs_vars) {
    id_tab <- table(data[[group]])
    if (!all(id_tab == id_tab[1])) {
        stop_("Unequal number of time points")
    }
    T <- as.integer(id_tab[1])
    if (missing(time)) {
        time_var <- 1:T
    } else {
        # TODO is there a better way?
        time_var <- sort(unique(data[[deparse(substitute(time))]]))
    }
    # TODO take spline definition into account here
    knots <- seq(time_var[1], time_var[T], length.out = min(10, T))
    Bs <- t(splines::bs(time_var, knots = knots, degree = 3, intercept = TRUE))
    D <- nrow(Bs)
    N <- length(id_tab)
    K <- ncol(moma)
    C <- length(get_resp(formula))
    moma_split <- split(moma, gl(T, 1, N * T))
    X <- aperm(array(as.numeric(unlist(moma_split)), dim = c(N, K, T)), c(3, 1, 2))
    assigned <- attr(moma, "assign")
    channel_vars <- list()
    for (i in seq_along(formula)) {
        pred_pos <- which(all_rhs_vars %in% formula[[i]]$predictors)
        pred_ix <- which(assigned %in% pred_pos)
        channel_vars[[paste0("J_", i)]] <- pred_ix
        channel_vars[[paste0("K_", i)]] <- length(pred_ix)
        resp_split <- split(data[[formula[[i]]$response]], group)
        Y <- array(as.numeric(unlist(resp_split)), dim = c(T, N))
        channel_vars[[paste0("response_", i)]] <- Y
        if (is_categorical(formula[[i]]$family)) {
            resp_pos <- which(all_rhs_vars %in% formula[[i]]$response)
            resp_ix <- which(assigned %in% resp_pos)
            channel_vars[[paste0("S_", i)]] <- length(resp_ix) + 1L
        }
    }
    c(named_list(T, N, C, K, X, D, Bs), channel_vars)
}

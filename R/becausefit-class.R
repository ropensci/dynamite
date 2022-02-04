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
    all_rhs <- reformulate(unique(c(unlist(formula$pred), unlist(formula$resp))),
                           intercept = FALSE)
    moma <- model.matrix(all_rhs, data)
    model_data <- convert_data(formula, data, group, time, moma)
    model_code <- create_blocks(formula)
    model <- rstan::stan_model(model_code = model_code)
    out <- rstan::sampling(model, data = model_data, ...)
    class(out) <- "becausefit"
    out
}

# Convert data for Stan
# TODO might not need to copy data here, investigate
convert_data <- function(formula, data, group, time,  moma) {
    id_tab <- table(data[[group]])
    if (!all(id_tab == id_tab[1])) {
        stop_("Unequal number of time points")
    }
    T <- as.integer(id_tab[1])
    if (missing(time)) {
        time <- 1:T
    }
    N <- length(id_tab)
    K <- ncol(moma)
    C <- length(formula$resp)
    moma_split <- split(moma, gl(T, 1, N * T))
    X <- aperm(array(unlist(moma_split), dim = c(N, K, T)), c(3, 1, 2))
    pred_names <- names(moma)
    channel_vars <- list()
    for (i in seq_along(formula$resp)) {
        pred_i <- which(pred_names %in% formula$pred[[i]])
        channel_vars[[paste0("L_", i)]] <- pred_i
        channel_vars[[paste0("K_", i)]] <- length(pred_i)
        if (formula$families[[i]] == "categorical") {
            channel_vars[[paste0("S_", i)]] <- length(unique(moma[,formula$resp[[i]]]))
        }
    }
    c(named_list(T, N, K, C, X), channel_vars)
}

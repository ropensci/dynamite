#' @include utilities.R

# becausefit class
becausefit <- function(formula, data, ...) {
    all_rhs <- reformulate(unique(c(unlist(formula$pred), unlist(formula$resp))),
                           intercept = FALSE)
    print(all_rhs)
    moma <- model.matrix(all_rhs, data)
    model_data <- convert_data(formula, data, moma)
    model_code <- create_blocks(formula)
    model <- rstan::stan_model(model_code = model_code)
    out <- rstan::sampling(model, data = model_data, ...)
    class(out) <- "becausefit"
    out
}

# Convert data for Stan
# TODO might not need to copy data here, investigate
convert_data <- function(formula, data, moma) {
    # TODO need some canonical way to define individuals
    # For now, assume that data is in long format with
    # column ID identifying the individuals
    id_tab <- table(data$ID)
    if (!all(id_tab == id_tab[1])) {
        stop_("Unequal number of time points")
    }
    T <- as.integer(id_tab[1])
    N <- length(id_tab)
    M <- formula$hidden$n
    K <- ncol(moma)
    C <- length(formula$resp)
    T_fixed_start <- formula$hidden$fixed[1]
    T_fixed_end <- formula$hidden$fixed[2]
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
    c(named_list(T, N, M, K, C, T_fixed_start, T_fixed_end, X),
      channel_vars)
}

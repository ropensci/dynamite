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
        stop_("Argument 'time' is missing.")
    time_var <- deparse(substitute(time))
    # Pipe for readability, not really needed if we need to support older R versions
    data <- data |>
        # Convert character types to factors
        dplyr::mutate(dplyr::across(tidyselect:::where(is.character), as.factor)) |>
        dplyr::arrange(data, dplyr::across(dplyr::all_of(c(group_var, time_var))))
    time <- unique(data[[time_var]])
    resp_all <- get_resp(formula)
    n_rows <- nrow(data)
    n_resp <- length(resp_all)
    lag_map <- extract_lags(unlist(get_pred(formula)))
    unprocessed_lags <- rep(TRUE, nrow(lag_map))
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
        }
        unprocessed_lags[lag_map$def %in% resp_all & lag_map$k <= fixed] <- FALSE
    }
    # Process lag terms defined via lag() in formulas
    if (any(unprocessed_lags)) {
        fixed <- max(max(lag_map$k), fixed)
        for (i in which(unprocessed_lags)) {
            if (lag_map$k[i] <= 0) {
                stop_("Only positive shift values are allowed in lag().")
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
            }
        }
    }
    for (j in seq_len(n_resp)) {
        formula[[j]]$formula <- reformulate(termlabels = formula[[j]]$predictors,
                                            response = resp_all[j],
                                            intercept = attr(terms(formula[[j]]$formula), "intercept"))
    }
    responses <- data[, resp_all, drop = FALSE]
    attr(responses, "resp_class") <- apply(responses, 2, class)
    model_matrix <- full_model.matrix(formula, data)
    resp_levels <- lapply(droplevels(responses), levels)
    # TODO: simplify I(lag(variable, 1)) to something shorter, e.g. lag_1(variable)?
    # TODO: shorten variable and or channel name if they are very long?
    # TODO: NOTE! you can use lag_map to get the variables within the complicated definition for formatting
    u_names <- colnames(model_matrix)
    coef_names <- lapply(seq_along(resp_all), function(i) {
        x <- paste0(resp_all[i], "_", u_names[attr(model_matrix, "assign")[[i]]])
        if (is_categorical(formula[[i]]$family)) {
            x <- paste0(x, "_", rep(resp_levels[[i]][-length(resp_levels[[i]])], each = length(x)))
        }
        x
    })
    specials <- evaluate_specials(formula, data)
    converted <- convert_data(formula, responses, specials, group, time, fixed, model_matrix)
    model_data <- converted$data
    model_helpers <- converted$helpers
    model_code <- create_blocks(formula, indent = 2L, resp = resp_all, helpers = model_helpers)
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
            time = time,
            time_var = time_var,
            group_var = group_var,
            levels = resp_levels,
            specials = specials,
            # TODO: extract only D for as.data.frame and J for predict
            #model_data = model_data,
            data = data,
            spline = list(B = model_data$Bs,
                D = model_data$D),
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

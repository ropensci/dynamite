#' @export
predict.btvcmfit <- function(object, newdata = NULL, mode = c("counterfactual", "forecast"),
    type = c("response", "mean", "link"), n_draws = NULL) {
    mode <- match.arg(mode)
    type <- match.arg(type)
    do.call(paste0("predict.btvcmfit_", mode), list(object = object, newdata = newdata, type = type, n_draws = n_draws))
}

predict.btvcmfit_counterfactual <- function(object, newdata, type, n_draws = NULL) {

    if (is.null(n_draws)) {
        n_draws <- ndraws(object)
    }
    # TODO needs to support different starting value? Although nothing is done for the time points without NA so perhaps not necessary?

    if (is.null(newdata)) {
        newdata <- object$data
        group <- newdata[[object$group_var]]
        time <- object$time
    } else {
        if (!(object$group_var %in% names(newdata)))
            stop(paste("Grouping variable", object$group_var, "not found in 'newdata'."))
        group <- newdata[[object$group_var]]
        if (is.factor(group)) group <- droplevels(group)
        # TODO doesn't really matter at least at the moment
        if (!all(group %in% object$data[[object$group_var]]))
            stop(paste("Grouping variable", object$group_var, "contains new levels not found in the original data."))

        if (!(object$time_var %in% names(newdata)))
            stop(paste("Time index variable", object$time_var, "not found in 'newdata'."))
        # sort just in case
        newdata <- dplyr::arrange(newdata, dplyr::across(dplyr::all_of(c(object$group_var, object$time_var))))
        # TODO just use the original time points starting from start_time
        time <- unique(newdata[[object$time_var]])
        if (!all(time %in% object$time))
            stop(paste("Timing variable", object$time_var, "contains time points not found in the original data."))
    }

    basis <- object$prediction_basis
    fixed <- basis$fixed
    n_time <- length(time)
    n_id <- length(unique(group))
    if (n_time <= fixed) {
        stop_("Model definition implies ", fixed, " fixed time points, but 'newdata' has only ", n_time, " time points.")
    }
    resp_all <- get_resp(basis$formula)
    # TODO check that fixed time points do not contain NAs
    for (resp in resp_all) {
        if (is.null(newdata[[resp]])) {
            newdata[[resp]] <- NA
        }
    }

    if (type != "response") {
        # create separate column for each level of categorical variables
        for (i in seq_along(resp_all)) {
            resp <- resp_all[i]
            if (is_categorical(basis$formula[[i]]$family)) {
                resp_levels <- object$levels[[resp]]
                newdata[, c(glue::glue("{resp}_{resp_levels}"))] <- NA # TODO: glued names to formula?
            } else {
                newdata[[c(glue::glue("{resp}_store"))]] <- newdata[[resp]]
            }
        }
    }
    newdata <- data.frame(newdata, draw = rep(1:n_draws, each = nrow(newdata)))
    samples <- rstan::extract(object$stanfit)
    u_names <- unique(names(basis$start))
    n <- nrow(newdata)
    for (i in (fixed + 1):n_time) {
        idx <- rep(seq(i - fixed, n, by = n_time), 1 + fixed) + rep(0:fixed, each = n_id * n_draws)
        idx_i <- seq(i, n, by = n_time)
        model_matrix <- full_model.matrix_fast(basis$formula, newdata[idx, ], u_names)[-1,] #remove extra due to NAs
        for (j in seq_along(resp_all)) {
            resp <- resp_all[j]
            if (any(is.na(newdata[idx_i, resp]))) { #TODO partial missingness?

                sim <- do.call(paste0("predict_", basis$formula[[j]]$family),
                    list(model_matrix = model_matrix[, basis$J[[j]], drop = FALSE],
                        samples = samples, resp, time =  i - fixed, type))


                if (is_categorical(basis$formula[[j]]$family)) {
                    resp_levels <- object$levels[[resp]]
                    if (type != "response") {
                        idx_resp <- which(names(newdata) %in% c(glue::glue("{resp}_{resp_levels}")))
                        newdata[idx_i, idx_resp] <- sim$mean_or_link
                    }
                    newdata[idx_i, resp] <- resp_levels[sim$response]
                } else {
                    if (type != "response") {
                        newdata[idx_i, c(glue::glue("{resp}_store"))] <- sim$mean_or_link
                    }
                    newdata[idx_i, resp] <- sim$response
                }
            }
        }
    }
    if (type != "response") {
        for (i in seq_along(resp_all)) {
            resp <- resp_all[i]
            if (is_categorical(basis$formula[[i]]$family)) {
                newdata[[resp]] <- NULL
            } else {
                newdata[[resp]] <- newdata[[c(glue::glue("{resp}_store"))]]
                newdata[[c(glue::glue("{resp}_store"))]] <- NULL
            }
        }
    }

    # TODO return either full newdata or just the responses
    #newdata[, resp_all, drop = FALSE]
    newdata
}

predict.btvcmfit_forecast <- function(object, newdata, type, n_draws) {
    stop_("Forecasting is not yet supported")
}

#' @export
predict.btvcmfit <- function(object, newdata = NULL, mode = c("counterfactual", "forecast"), n_draws = NULL) {
    mode <- match.arg(mode)
    do.call(paste0("predict.btvcmfit_", mode), list(object = object, newdata = newdata, n_draws = n_draws))
}
# TODO: What about the posterior distribution at the means/probability level?
predict.btvcmfit_counterfactual <- function(object, newdata, n_draws = NULL) {

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
    n_resp <- length(resp_all)
    # TODO check that first fixed time points do not contain NAs
    for (resp in resp_all) {
        if (is.null(newdata[[resp]])) {
            newdata[[resp]] <- NA
        }
    }
    newdata <- data.frame(newdata, draw = rep(1:n_draws, each = nrow(newdata)))
    samples <- rstan::extract(object$stanfit)
    u_names <- unique(names(basis$start))
    n <- nrow(newdata)
    for (i in (fixed + 1):n_time) {
        idx <- rep(seq(i - fixed, n, by = n_time), 1 + fixed) + rep(0:fixed, each = n_id * n_draws)
        idx_i <- seq(i, n, by = n_time)
        # TODO: revert back to switch when we support for unequal series
        # idx <- which(dplyr::between(newdata[[object$time_var]], i - fixed, i)) # TODO don't assume time_var is 1,2,...
        # idx_i <- which(newdata[[object$time_var]] == i)

        model_matrix <- full_model.matrix_fast(basis$formula,
            newdata[idx, ], u_names)[-1,] #remove extra due to NAs
        # drop previous time points
        # model_matrix <- model_matrix[-seq(i, nrow(model_matrix), length = n_id * n_draws), ]
        for (j in seq_along(resp_all)) {
            resp <- resp_all[j]
            if (any(is.na(newdata[idx_i, resp]))) { #TODO partial missingness?
                # TODO instead of if clauses call rng_xx(xbeta, samples, k)
                # repeated assignments to data.frame row is slow
                if (is_gaussian(basis$formula[[j]]$family)) {
                    sim <- numeric(length(idx_i))
                    for(k in seq_len(n_draws)) {
                        idx_k <- ((k - 1) * n_id + 1):(k * n_id)
                        xbeta <- model_matrix[idx_k, basis$J[[j]], drop = FALSE] %*% samples[[paste0("beta_", j)]][k, i - fixed, ]
                        sigma <- samples[[paste0("sigma_", j)]][k]
                        sim[idx_k] <- rnorm(n_id, xbeta, sigma)
                    }
                    newdata[idx_i, resp] <- sim
                }
                if (is_categorical(basis$formula[[j]]$family)) {
                    #TODO check that the order is correct
                    y_levels <- levels(droplevels(object$data[[resp]])) #TODO, could be already stored in basis$formula$levels
                    S <- length(y_levels)
                    sim <- integer(length(idx_i))
                    for(k in seq_len(n_draws)) {
                        idx_k <- ((k - 1) * n_id + 1):(k * n_id)
                        xbeta <- model_matrix[idx_k, basis$J[[j]], drop = FALSE] %*% samples[[paste0("beta_", j)]][k, i - fixed, , ]
                        sim[idx_k] <- max.col(xbeta - log(-log(runif(S * n_id))))
                    }
                    newdata[idx_i, resp] <- y_levels[sim]
                }
                # for(k in seq_len(n_draws)) {
                #     idx_k <- ((k - 1) * n_id + 1):(k * n_id)
                #     xbeta <- model_matrix[idx_k, basis$J[[j]], drop = FALSE] %*% samples[[paste0("beta_", j)]][k, i - fixed, ]
                #     if (is_gaussian(basis$formula[[j]]$family)) {
                #         sigma <- samples[[paste0("sigma_", j)]][k]
                #         sim <- rnorm(n_id, xbeta, sigma)
                #         newdata[idx_i[idx_k], resp] <- sim
                #     } else {
                #         if (is_categorical(basis$formula[[j]]$family)) {
                #             y_levels <- basis$formula$levels #TODO
                #             newdata[start:end, resp] <- do.call(basis$rngs[[resp]],
                #                 list(y_levels, model_matrix[i, basis$J], samples[[paste0("beta_", j)]][k, i, ,]))
                #         }
                #     }
                # }
            }
        }
    }


    # TODO return either full newdata or just the responses
    #newdata[, resp_all, drop = FALSE]
    newdata
}

predict.btvcmfit_forecast <- function(object, newdata, n_draws) {
    stop_("Forecasting is not yet supported")
}

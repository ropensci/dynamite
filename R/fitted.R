#' Extract fitted values of btvcmfit
#'
#' TODO Note that these are conditional on the observed data i.e., we don't
#' simulate new lagged values for covariates, so we underestimate the uncertainty.
#' Better to use predict with type = "mean" TODO do it
#' @export
#' @param object An object of class \code{btvcmfit}.
#' @param ... Ignored.
#' @importFrom stats fitted
fitted.btvcmfit <- function(object, newdata = NULL, n_draws = NULL, ...) {

    if (is.null(n_draws)) {
        n_draws <- ndraws(object)
    }

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
        newdata <- dplyr::arrange(newdata, dplyr::across(object$time_var))
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

    samples <- rstan::extract(object$stanfit)
    u_names <- unique(names(basis$start))

    model_matrix <- full_model.matrix_fast(basis$formula,
        newdata, u_names)
    # create separate column for each level of categorical variables
    for (i in seq_along(resp_all)) {
        resp <- resp_all[i]
        if (is_categorical(basis$formula[[i]]$family)) {
            resp_levels <- object$levels[[resp]]
            newdata[, c(glue::glue("{resp}_{resp_levels}"))] <- NA # TODO: glued names to formula?
            newdata[[resp]] <- NULL
        }
    }
    newdata <- data.frame(newdata, draw = rep(1:n_draws, each = nrow(newdata)), check.names = FALSE)
    n <- nrow(newdata)
    for (i in (fixed + 1):n_time) {

        idx_i <- seq(i, n, by = n_time) # TODO fix for a case of unequal number of time points per group
        idx_i2 <- seq(i, nrow(model_matrix), by = n_time)
        for (j in seq_along(resp_all)) {
            resp <- resp_all[j]
            if (is_gaussian(basis$formula[[j]]$family)) {
                xbeta <- model_matrix[idx_i2, basis$J[[j]], drop = FALSE] %*% t(samples[[paste0("beta_", j)]][1:n_draws, i - fixed, ])
                newdata[idx_i, resp] <- c(xbeta)
            }
            if (is_categorical(basis$formula[[j]]$family)) {
                sim <- matrix(NA, length(idx_i), length(object$levels[[resp]]))
                for(k in seq_len(n_draws)) {
                    idx_k <- ((k - 1) * n_id + 1):(k * n_id)
                    xbeta <- model_matrix[idx_i2, basis$J[[j]], drop = FALSE] %*% samples[[paste0("beta_", j)]][k, i - fixed, , ]
                    # sim[idx_k, ] <- exp(xbeta)/rowSums(exp(beta))
                    maxs <- apply(xbeta, 1, max)
                    sim[idx_k, ] <- exp(xbeta - (maxs + log(rowSums(exp(xbeta - maxs)))))
                }
                idx_resp <- which(names(newdata) %in% c(glue::glue("{resp}_{resp_levels}")))
                newdata[idx_i, idx_resp] <- sim
            }

        }
    }


    newdata
}

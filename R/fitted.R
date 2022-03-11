#' Extract fitted values of btvcmfit
#'
#' TODO Note that these are conditional on the observed data i.e., we don't simulate new lagged values for covariates
#'
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
        group <- newdata[, object$group_var]
        time <- object$time
    } else {
        if (!(object$group_var %in% names(newdata)))
            stop(paste("Grouping variable", object$group_var, "not found in 'newdata'."))
        group <- newdata[, object$group_var]
        if (is.factor(group)) group <- droplevels(group)
        # TODO doesn't really matter at least at the moment
        if (!all(group %in% object$data[, object$group_var]))
            stop(paste("Grouping variable", object$group_var, "contains new levels not found in the original data."))

        if (!(object$time_var %in% names(newdata)))
            stop(paste("Time index variable", object$time_var, "not found in 'newdata'."))
        # sort just in case
        newdata <- dplyr::arrange(newdata, dplyr::across(object$time_var))
        # TODO just use the original time points starting from start_time
        time <- unique(newdata[, object$time_var])
        if (!all(time %in% object$time))
            stop(paste("Timing variable", object$time_var, "contains time points not found in the original data."))
    }

    # For each ID independently
    if (length(unique(group)) > 1) {
        newdata <- dplyr::bind_rows(lapply(split(newdata, group), function(x)
            fitted.btvcmfit(object, x, n_draws)))
    } else {
        basis <- object$prediction_basis
        fixed <- basis$fixed
        resp_all <- get_resp(basis$formula)
        responses <- newdata[, resp_all, drop = FALSE]
        model_matrix <- full_model.matrix(basis$formula, newdata, resp_all)
        model_data <- convert_data(basis$formula, responses, group, time, fixed, model_matrix)
        samples <- rstan::extract(object$stanfit)
        # extract returns permuted samples so we can just pick first n_draws
        #idx <- 1:n_draws #sample(ndraws(bf), size = n_draws)

        # remove fixed time points (TODO: Fix, now assumes time is of form 1,2,...)
        newdata <- dplyr::filter(newdata, time > fixed)
        # create separate column for each level of categorical variables
        for (i in seq_along(resp_all)) {
            resp <- resp_all[i]
            if (is_categorical(basis$formula[[i]]$family)) {
                resp_levels <- basis$formula[[i]]$levels
                newdata[, c(glue::glue("{resp}[{resp_levels}]"))] <- NA # TODO: glued names to formula?
                newdata[[resp]] <- NULL
            }
        }
        n <- nrow(newdata)
        newdata <- data.frame(newdata, draw = rep(1:n_draws, each = nrow(newdata)))
        for (i in seq_len(n)) {
            for (j in seq_along(resp_all)) {
                resp <- resp_all[j]
                if (is_gaussian(basis$formula[[j]]$family)) {
                    for(k in seq_len(n_draws)) {
                        newdata[i + (k - 1) * n, resp] <- do.call(basis$mean[[resp]],
                            list(model_matrix[i, basis$J], samples[[paste0("beta_", j)]][k, i, ]))
                    }
                } else {
                    if (is_categorical(basis$formula[[j]]$family)) {
                        resp_levels <- basis$formula[[i]]$levels
                        idx <- which(names(newdata) %in% c(glue::glue("{resp}[{resp_levels}]")))
                        for(k in seq_len(n_draws)) {
                            newdata[i + (k - 1) * n, idx] <- do.call(basis$mean[[resp]],
                                list(model_matrix[i, basis$J], samples[[paste0("beta_", j)]][k, i, ,]))
                        }
                    }
                }

            }

        }
    }
    newdata
}

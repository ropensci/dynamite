#' @export
predict.btvcmfit <- function(object, newdata, mode = c("counterfactual", "forecast"), n_draws = NULL) {
    mode <- match.arg(mode)
    do.call(paste0("predict.btvcmfit_", mode), list(bf = object, newdata = newdata, n_draws = n_draws))
}
# TODO: What about the posterior distribution at the means/probability level?
predict.btvcmfit_counterfactual <- function(bf, newdata, n_draws = NULL) {

    if (is.null(n_draws)) {
        n_draws <- ndraws(bf)
    }
    # TODO needs some logic and testing for time points, does time points in newdata match with original?
    # could be a (continuous? or can be made continuous) subset, e.g. original time was 1:100
    # we might want to predict for times 2:100, 5:50, 2:10 etc

    # TODO do we assume that newdata corresponds to one ID only? Simplifies things,
    # but one might be interested in marginal effects so that we could use original individuals (multiple ids)
    # and their exogenous covariates
    # We could also just call this function multiple times after grouping by ID
    basis <- bf$prediction_basis
    fixed <- basis$fixed
    n_new <- nrow(newdata)
    if (n_new <= fixed) {
        stop_("Model definition implies ", fixed, " fixed time points, but 'newdata' has only ", n_new, " rows")
    }
    resp_all <- get_resp(basis$formula)
    n_resp <- length(resp_all)
    for (resp in resp_all) {
        if (is.null(newdata[[resp]])) {
            newdata[[resp]] <- NA
        }
    }
    newdata <- newdata[, basis$ord]
    samples <- rstan::extract(bf$stanfit)
    # Prediction loop
    idx <- sample(ndraws(bf), size = n_draws)
    # create list-columns for responses, there's probably better way
    newdata <- tibble::as_tibble(newdata)
    for (resp in resp_all) {
        newdata[[resp]] <- lapply(1:nrow(newdata), function(i) rep(newdata[[resp]][i], n_draws))
    }
    for(k in 1:n_draws) {
        for (i in seq_len(n_new)) {
            newdata_k <- slice_tibble(newdata[(i-fixed):i,], k)
            model_matrix <- full_model.matrix(basis$formula, newdata_k, resp_all)
            for (j in seq_along(resp_all)) {
                resp <- resp_all[j]
                if (is.na(newdata_k[i,resp])) {
                    # TODO instead of if clauses call rng(model_matrix[, J_i], samples, idx[k])
                    # TODO except we need additional arguments depending on the family
                    if (is_gaussian(basis$formula[[j]]$family)) {
                        newdata[[resp]][[i]][k] <- do.call(basis$rng[[resp]], # TODO match time!
                            list(model_matrix[i, basis$J], samples[[paste0("beta_", j)]][idx[k], i - 1, ], samples[[paste0("sigma_", j)]][idx[k]]))
                    } else {
                        if (is_categorical(basis$formula[[j]]$family)) {
                            y_levels <- basis$formula$levels #TODO
                                newdata[[resp]][[i]][k] <- do.call(basis$rngs[[resp]],
                                list(y_levels, model_matrix[i, basis$J], samples[[paste0("beta_", j)]][idx[k], i - 1, ,]))
                        }
                    }
                }
            }
        }
    }
    # TODO return either full newdata or just the responses
    #newdata[, resp_all, drop = FALSE]
    newdata
}

predict.btvcmfit_forecast <- function(bf, newdata, n_draws) {
    stop_("Forecasting is not yet supported")
}

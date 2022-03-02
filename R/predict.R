#' @export
predict.btvcmfit <- function(object, newdata, mode = c("counterfactual", "forecast")) {
    mode <- match.arg(mode)
    do.call(paste0("predict.btvcmfit_", mode), list(bf = object, newdata = newdata))
}

predict.btvcmfit_counterfactual <- function(bf, newdata) {
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
    # TODO: This is now for a single sample for posterior, we would like to more
    for (i in seq_len(n_new)) {
        model_matrix <- full_model.matrix(basis$formula, newdata[(i-fixed):i,], resp_all)
        for (j in seq_along(resp_all)) {
            resp <- resp_all[j]
            if (is.na(newdata[i,resp])) {
                # TODO actual model fit based prediction logic, just something for now
                if (is_gaussian(basis$formula$family)) {
                    #TODO loop over iters
                    #TODO select correct columns using J_1 etc
                    #TODO expose stan functions when building the model?
                  newdata[i,resp] <- do.call(paste0("basis$rngs$response_", j, "_rng"),
                      list(model_matrix, samples[[paste0("beta_",j)]][1,i,], samples[[paste0("sigma_",j)]][1]))
                }
            }
        }
    }
    newdata[, resp_all, drop = FALSE]
}

predict.btvcmfit_forecast <- function(bf, newdata) {
    stop_("Forecasting is not yet supported")
}

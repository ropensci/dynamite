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
    # Prediction loop
    for (i in seq_len(n_new)) {
        model_matrix <- full_model.matrix(basis$formula, newdata[(i-fixed):i,])
        for (resp in resp_all) {
            if (is.na(newdata[i,resp])) {
                # TODO actual model fit based prediction logic, just something for now
                newdata[i,resp] <- rnorm(1)
            }
        }
    }
    newdata[,resp_all]
}

predict.btvcmfit_forecast <- function(bf, newdata) {
    stop_("Forecasting is not yet supported")
}

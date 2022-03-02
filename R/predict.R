#' @export
predict.btvcmfit <- function(object, newdata, mode = c("counterfactual", "forecast"), n_samples = 100) {
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
    idx <- sample(rstantools::nsamples(bf$stanfit), size = n_samples)
    #TODO: Is it easier to just create tibble with list-columns afterwards by combining separate normal data frames?
    # create list-columns for responses, there's probably better way
    newdata <- tibble::as_tibble(newdata)
    for (resp in resp_all) {
        newdata[[resp]] <- list(rep(newdata[[resp]], n_samples))
    }
    for(k in 1:n_samples) {
        for (i in seq_len(n_new)) {
            # TODO: how do you extract rows of tibble with list-columns converted to normal columns based on ith element in said list?
            model_matrix <- full_model.matrix(basis$formula, newdata[(i-fixed):i,], resp_all)
            for (resp in resp_all) {
                if (is.na(newdata[i,resp])) {
                    if (is_gaussian(basis$formula$family)) {
                        #TODO select correct columns using J_1 etc
                        newdata[i,resp] <- do.call(basis$rngs[[resp]],
                            list(model_matrix, samples[[paste0("beta_",j)]][1,i,], samples[[paste0("sigma_",j)]][1]))
                    }
                }
            }
        }
    }
    newdata[, resp_all, drop = FALSE]
}

predict.btvcmfit_forecast <- function(bf, newdata) {
    stop_("Forecasting is not yet supported")
}


categorical_rng <- function(y, x, beta, J) {
    sample(y, 1, prob = x[J] %*% beta)
}
gaussian_rng <- function(x, beta, sigma, J) {
    rnorm(1, x[J] %*% beta, sigma)
}
create_predict_functions <- function(x) {
    resp_all <- get_resp(x)
    rng <- vector("list", length(resp_all))
    names(rng) <- resp_all
    families <- lapply(x, function(x) x$family$name)
    for (i in seq_along(resp_all)) {
        rng[[i]] <- get(paste0(families[i], "_rng"))
    }
    rng
}

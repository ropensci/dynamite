#' @export
predict.btvcmfit <- function(object, newdata, response) {
    basis <- object$prediction_basis
    n_resp <- length(basis$map)
    all_lag_vars <- unique(unlist(lapply(basis$map, names)))
    fixed <- nrow(basis$past)
    n_new <- nrow(newdata)
    resp_all <- get_resp(basis$formula)
    resp_data <- setNames(replicate(n_resp, rep(0, n_new)), resp_all)
    newdata <- cbind(newdata, resp_data)[object$ord,]
    pred_data <- rbind(pred_past, newdata)
    pred_resp <- replicate(n_resp, numeric(n_new), simplify = FALSE)
    formulas <- get_form(basis_formula)

    u_names <- unique(colnames(model_matrix))
    # TODO: Is intercept always named as (Intercept)?
    # checking assign attributes for 0 is an another option
    if (any(ind <- u_names == "(Intercept)")) {
        u_names <- u_names[-which(ind)]
    }
    model_matrix <- model_matrix[, u_names]
    # Prediction loop
    for (i in (fixed + seq_len(n_new))) {
        model_matrices <- lapply(formulas, function(x) {
            model.matrix(x, pred_data[(i - fixed):i,])[fixed + 1, drop = FALSE]
        }
        model_matrix <- do.call(cbind, model_matrices)
        for (j in seq_len(n_resp)) {
            # TODO actual prediction logic
            pred_data[]
        }
    }
}

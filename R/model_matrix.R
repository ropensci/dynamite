#' Combine model.matrix objects of all formulas of a dynamiteformula into one
#'
#' @param dformula A `dynamiteformula` object
#' @param data A `data.table` containing the variables in the model
#'
#' @srrstats {RE1.3, RE1.3a} `full_model.matrix` preserves relevant attributes.
#' @noRd
full_model.matrix <- function(dformula, data) {
  formulas <- get_formulas(dformula)
  model_matrices <- vector(mode = "list", length = length(formulas))
  for (i in seq_along(formulas)) {
    mm <- model.matrix.lm(formulas[[i]], data = data, na.action = na.pass)
    test_collinearity(dformula[[i]]$resp, mm, data)
    # Intercept is not part of X
    model_matrices[[i]] <- remove_intercept(mm)
  }
  model_matrix <- do.call(cbind, model_matrices)
  u_names <- unique(colnames(model_matrix))
  model_matrix <- model_matrix[, u_names, drop = FALSE]
  n_models <- length(model_matrices)
  y_names <- get_responses(dformula)
  empty_list <- setNames(vector(mode = "list", length = n_models), y_names)
  attr(model_matrix, "assign") <- empty_list
  attr(model_matrix, "fixed") <- empty_list
  attr(model_matrix, "varying") <- empty_list
  for (i in seq_along(model_matrices)) {
    cols <- colnames(model_matrices[[i]])
    assign <- match(cols, u_names)
    assign_i <- attr(model_matrices[[i]], "assign")
    attr(model_matrix, "assign")[[i]] <- sort(assign)
    fixed <- assign[assign_i %in% dformula[[i]]$fixed]
    attr(model_matrix, "fixed")[[i]] <- setNames(fixed, u_names[fixed])
    varying <- assign[assign_i %in% dformula[[i]]$varying]
    attr(model_matrix, "varying")[[i]] <- setNames(varying, u_names[varying])
  }
  model_matrix
}

#' Test collinearity within a channel
#'
#' @param y The response variable of the channel
#' @param mm Model matrix based on the channel model formula
#' @param data A `data.table` containing the variables in the model
#' @srrstats {BS3.1, BS3.2, RE2.4, RE2.4a, RE2.4b} Collinearity is tested.
#' @noRd
test_collinearity <- function(y, mm, data) {
  nc <- ncol(mm)
  mm_obs <- stats::complete.cases(mm)
  # check for n < p
  n <- sum(mm_obs)
  if (n < nc) {
    warning_(
      "Number of non-missing observations {sum(mm_obs)} in channel {.var {y}}
      is less than {nc}, the number of predictors (including possible
      intercept)."
    )
  }
  mm_names <- colnames(mm)
  if (any(mm_obs) && !identical(qr(mm[mm_obs, ])$rank, min(n, nc))) {
    zero_col <- apply(mm[mm_obs, , drop = FALSE], 2L, function(x) all(x == 0))
    if (any(zero_col)) {
      k <- sum(zero_col)
      warning_(
        "{cli::qty(k)} Predictor{?s} {.var {cs(mm_names[zero_col])}}
         {cli::qty(k)} contain{?s/} only zeros in the complete case rows
         of the design matrix for the channel {.var {y}}."
      )
    } else {
      warning_(
        "Perfect collinearity found between predictor variables of
         channel {.var {y}}."
      )
    }
  }
  for (j in seq_len(nc)) {
    pred_resp <- cbind(as.numeric(mm[, j]), as.numeric(data[[y]]))
    pred_resp_obs <- stats::complete.cases(pred_resp)
    if (any(pred_resp_obs) &&
        !identical(qr(pred_resp[pred_resp_obs, ])$rank, 2L)) {
      warning_(c(
        "Perfect collinearity found between response and predictor variable:",
        `i` = "Response variable {.var {y}} is perfectly
               collinear with predictor variable {.var {mm_names[j]}}."
      ))
    }
  }
}

#' A streamlined version of full_model.matrix for prediction
#'
#' @param formula_list A `list` of `formula` objects
#' @param newdata A `data.frame` containing the variables in the model
#' @param idx An `integer` vector of row indices to subset by
#' @param u_names A `character` vector of unique column names of the resulting
#'   matrix
#' @noRd
full_model.matrix_predict <- function(formula_list, newdata, idx, u_names) {
  newdata_sub <- newdata[idx, ]
  model_matrices <- lapply(formula_list, function(x) {
    model.matrix.lm(x, newdata_sub, na.action = na.pass)
  })
  model_matrices <- lapply(model_matrices, remove_intercept)
  model_matrix <- do.call(cbind, model_matrices)
  model_matrix[, u_names, drop = FALSE]
}

#' Remove Intercept from the Model Matrix
#'
#' @param x A model matrix from `model.matrix.lm`
#' @noRd
remove_intercept <- function(x) {
  idx <- which(attr(x, "assign") == 0L)
  if (length(idx) > 0L) {
    a <- attr(x, "assign")[-idx]
    x <- x[, -idx,  drop = FALSE]
    attr(x, "assign") <- a
  }
  x
}

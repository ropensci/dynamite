#' Combine model.matrix objects of all formulas of a dynamiteformula into one
#'
#' @param dformula A `dynamiteformula` object
#' @param data A `data.frame` containing the variables in the model
#'
#' @noRd
full_model.matrix <- function(dformula, data) {
  model_matrices <- lapply(get_formulas(dformula), model.matrix.lm,
                           data = data, na.action = na.pass)
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
    attr(model_matrix, "assign")[[i]] <- sort(assign)
    fixed <- assign[attr(model_matrices[[i]], "assign") %in% dformula[[i]]$fixed]
    attr(model_matrix, "fixed")[[i]] <- setNames(fixed, u_names[fixed])
    varying <- assign[attr(model_matrices[[i]], "assign") %in% dformula[[i]]$varying]
    attr(model_matrix, "varying")[[i]] <- setNames(varying, u_names[varying])
  }
  model_matrix
}

#' A version of full_model.matrix for prediction
#'
#' @param formula_list A `list` of `formula` objects
#' @param newdata A `data.frame` containing the variables in the model
#' @param u_names TODO
#'
#' @noRd
full_model.matrix_predict <- function(formula_list, newdata, idx, u_names) {
  newdata_sub <- as.data.frame(newdata[idx, ])
  model_matrices <- lapply(formula_list, function(x) {
    model.matrix.lm(x, newdata_sub, na.action = na.pass)
  })
  model_matrix <- do.call(cbind, model_matrices)
  model_matrix[, u_names, drop = FALSE]
}

#' A fast version of full_model.matrix using formulas directly
#'
#' @param formula_list A `list` of `formula` objects
#' @param data A `data.frame` containing the variables in the model
#' @param u_names A character vector of unique predictor names
#'
#' @noRd
full_model.matrix_fast <- function(formula_list, data, u_names) {
  model_matrices <- lapply(formula_list, model.matrix, data)
  model_matrix <- do.call(cbind, model_matrices)
  model_matrix[, u_names, drop = FALSE]
}

#' A pseudo version of full_model.matrix, where the evaluation is assumed
#' 'as.is', i.e., the model tilde is assumed to represent a mathematical
#' equality
#'
#' @param formula_list A list of `formula` objects
#' @param data A `data.frame` containing the variables in the model
#' @param initial
#'
#' @noRd
full_model.matrix_pseudo <- function(formula_list, data) {
  data_env <- list2env(data)
  model_matrices <- lapply(formula_list, eval_formula, envir = data_env)
  out <- do.call(cbind, model_matrices)
  if (nrow(out) == 1) {
    # TODO warn about recycling?
    out[rep(1, nrow(data)),]
  } else {
    out
  }
}

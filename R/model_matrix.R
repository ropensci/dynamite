#' Combine model.matrix objects of all formulas of a dynamiteformula into one
#'
#' @param dformula A `dynamiteformula` object
#' @param data A `data.table` containing the variables in the model
#'
#' @srrstats {RE1.3} *Regression Software which passes or otherwise transforms aspects of input data onto output structures should ensure that those output structures retain all relevant aspects of input data, notably including row and column names, and potentially information from other `attributes()`.*
#' @srrstats {BS3.1} *Implement pre-processing routines to diagnose perfect collinearity, and provide appropriate diagnostic messages or warnings*
#' @srrstats {BS3.2} *Provide distinct routines for processing perfectly collinear data, potentially bypassing sampling algorithms*
#' @srrstats {RE2.4} *Regression Software should implement pre-processing routines to identify whether aspects of input data are perfectly collinear, notably including:*
#' @srrstats {RE2.4a} *Perfect collinearity among predictor variables*
#' @srrstats {RE2.4b} *Perfect collinearity between independent and dependent variables*
#' @noRd
full_model.matrix <- function(dformula, data) {
  formulas <- get_formulas(dformula)
  model_matrices <- vector(mode = "list", length = length(formulas))
  for (i in seq_along(formula)) {
    y <- dformula[[i]]$resp
    mm <- model.matrix.lm(formulas[[i]], data = data, na.action = na.pass) |>
      remove_intercept()
    nc <- ncol(mm)
    mm_obs <- stats::complete.cases(mm)
    if (any(mm_obs) && !identical(qr(mm[mm_obs, ])$rank, nc)) {
      warning_(
        "Perfect collinearity found between predictor variables of
         channel {.var {y}}."
      )
    }
    mm_names <- colnames(mm)
    for (j in seq_len(ncol(mm))) {
      pred_resp <- cbind(as.numeric(mm[,j]), as.numeric(data[[y]]))
      pred_resp_obs <- stats::complete.cases(pred_resp)
      if (any(pred_resp_obs) &&
          !identical(qr(pred_resp[pred_resp_obs, ])$rank, 2L)) {
        x <- mm_names[j]
        warning_(c(
          "Perfect collinearity found between response and predictor variable:",
          `i` = "Response variable {.var {y}} is perfectly
                 collinear with predictor variable {.var {x}}."
        ))
      }
    }
    model_matrices[[i]] <- mm
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

#' A version of full_model.matrix for prediction
#'
#' @param formula_list A `list` of `formula` objects
#' @param newdata A `data.frame` containing the variables in the model
#' @param idx An `integer` vector of row indices to subset by
#' @param u_names A `character` vector of unique column names of the resulting
#'   matrix
#'
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

#' A fast version of full_model.matrix using formulas directly
#'
#' @param formula_list A `list` of `formula` objects
#' @param newdata A `data.frame` containing the variables in the model
#' @param u_names A character vector of unique predictor names
#'
#' @noRd
full_model.matrix_fast <- function(formula_list, newdata, u_names) {
  model_matrices <- lapply(formula_list, function(x) {
    model.matrix.lm(x, newdata, na.action = na.pass)
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
  idx <- which(attr(x, "assign") == 0)
  if (length(idx) > 0) {
    a <- attr(x, "assign")[-idx]
    x <- x[, -idx,  drop = FALSE]
    attr(x, "assign") <- a
  }
  x
}

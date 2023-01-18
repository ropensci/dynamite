#' Combine `model.matrix` Objects of All Formulas of a `dynamiteformula`
#'
#' @inheritParams dynamite
#' @srrstats {RE1.3, RE1.3a} `full_model.matrix` preserves relevant attributes.
#' @noRd
full_model.matrix <- function(dformula, data, verbose) {
  formulas <- get_formulas(dformula)
  model_matrices <- vector(mode = "list", length = length(formulas))
  for (i in seq_along(formulas)) {
    mm <-
      stats::model.matrix.lm(formulas[[i]], data = data, na.action = na.pass)
    if (verbose) {
      test_collinearity(dformula[[i]]$resp, mm, data)
    }
    # Intercept is not part of X
    model_matrices[[i]] <- remove_intercept(mm)
  }
  model_matrix <- do.call(cbind, model_matrices)
  u_names <- unique(colnames(model_matrix))
  model_matrix <- model_matrix[, u_names, drop = FALSE]
  y_names <- get_responses(dformula)
  empty_list <- setNames(
    vector(mode = "list", length = length(model_matrices)),
    y_names
  )
  attr(model_matrix, "assign") <- empty_list
  attr(model_matrix, "fixed") <- empty_list
  attr(model_matrix, "varying") <- empty_list
  attr(model_matrix, "random") <- empty_list
  for (i in seq_along(model_matrices)) {
    cols <- colnames(model_matrices[[i]])
    assign <- match(cols, u_names)
    assign_i <- attr(model_matrices[[i]], "assign")
    attr(model_matrix, "assign")[[i]] <- sort(assign)
    fixed <- assign[assign_i %in% dformula[[i]]$fixed]
    attr(model_matrix, "fixed")[[i]] <- setNames(fixed, u_names[fixed])
    varying <- assign[assign_i %in% dformula[[i]]$varying]
    attr(model_matrix, "varying")[[i]] <- setNames(varying, u_names[varying])
    random <- assign[assign_i %in% dformula[[i]]$random]
    attr(model_matrix, "random")[[i]] <- setNames(random, u_names[random])
  }
  model_matrix
}

#' Test Collinearity Within a Channel
#'
#' @param y \[`character(1)`]\cr The response variable of the channel.
#' @param mm \[`matrix`]\cr Model matrix based on the channel model formula.
#' @param data \[`data.table`]\cr Data containing the variables in the model.
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

#' A streamlined Version of `full_model.matrix` for Prediction
#'
#' @param formula_list \[`list`]\cr A `list` of `formula` objects.
#' @param newdata \[`data.table`]\cr Data containing the variables in the model.
#' @param idx \[`integer()`]\cr A vector of row indices to subset with.
#' @param u_names \[`character()`]\cr A vector of unique column names of
#'   the resulting matrix.
#' @noRd
full_model.matrix_predict <- function(formula_list, newdata_resp, newdata_pred,
                                      idx_resp, idx_pred, n_draws, u_names) {
  sub_resp <- newdata_resp[idx_resp, ]
  sub_pred <- newdata_pred[idx_pred, ]
  sub <- cbind(sub_resp, sub_pred)
  model_matrices <- lapply(formula_list, function(x) {
    stats::model.matrix.lm(x, sub, na.action = na.pass)
  })
  model_matrices <- lapply(model_matrices, remove_intercept)
  model_matrix <- do.call(cbind, model_matrices)
  model_matrix <- model_matrix[,
    !duplicated(colnames(model_matrix)),
    drop = FALSE
  ]
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
    x <- x[, -idx, drop = FALSE]
    attr(x, "assign") <- a
  }
  x
}

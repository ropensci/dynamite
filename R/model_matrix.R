#' Combine `model.matrix` Objects of All Formulas of a `dynamiteformula`
#'
#' @inheritParams prepare_stan_input
#' @srrstats {RE1.3, RE1.3a} `full_model.matrix` preserves relevant attributes.
#' @noRd
full_model.matrix <- function(dformula, data, group_var, fixed, verbose) {
  # avoid NSE notes from R CMD check
  group <- NULL
  model_matrices <- vector(mode = "list", length = length(dformula))
  model_matrices_type <- vector(mode = "list", length = length(dformula))
  types <- c("fixed", "varying", "random")
  idx <- data[,
    .I[base::seq.int(fixed + 1L, .N)],
    by = group,
    env = list(fixed = fixed, group = group_var)
  ]$V1
  data_nonfixed <- droplevels(data[idx, , env = list(idx = idx)])
  for (i in seq_along(dformula)) {
    mm <- stats::model.matrix.lm(
      dformula[[i]]$formula,
      data = data_nonfixed,
      na.action = na.pass
    )
    if (verbose) {
      test_collinearity(dformula[[i]]$resp, mm, data_nonfixed)
    }
    model_matrices_type[[i]] <- list()
    for (type in c("fixed", "varying", "random")) {
      type_formula <- get_type_formula(dformula[[i]], type)
      if (!is.null(type_formula)) {
        model_matrices_type[[i]][[type]] <- stats::model.matrix.lm(
          type_formula,
          data = data_nonfixed,
          na.action = na.pass
        )
      }
    }
    tmp <- do.call(cbind, model_matrices_type[[i]])
    ifelse_(identical(length(tmp), 0L),
      model_matrices[[i]] <- matrix(nrow = nrow(mm), ncol = 0L),
      model_matrices[[i]] <- tmp
    )
  }
  model_matrix <- do.call(cbind, model_matrices)
  u_names <- setdiff(unique(colnames(model_matrix)), "(Intercept)")
  model_matrix <- model_matrix[, u_names, drop = FALSE]
  y_names <- get_responses(dformula)
  empty_list <- stats::setNames(
    vector(mode = "list", length = length(model_matrices)),
    y_names
  )
  #attr(model_matrix, "assign") <- empty_list
  attr(model_matrix, "fixed") <- empty_list
  attr(model_matrix, "varying") <- empty_list
  attr(model_matrix, "random") <- empty_list
  for (i in seq_along(model_matrices)) {
    for (type in types) {
      if (!is.null(model_matrices_type[[i]][[type]])) {
        cols <- which(u_names %in% colnames(model_matrices_type[[i]][[type]]))
        attr(model_matrix, type)[[i]] <- stats::setNames(cols, u_names[cols])
      } else {
        attr(model_matrix, type)[[i]] <- integer(0L)
      }
    }
    # attr(model_matrix, "assign")[[i]] <- sort(
    #   unique(
    #     c(
    #       attr(model_matrix, "fixed")[[i]],
    #       attr(model_matrix, "varying")[[i]],
    #       attr(model_matrix, "random")[[i]]
    #     )
    #   )
    # )
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

#' A Streamlined Version of `full_model.matrix` for Prediction
#'
#' @param dformula \[`dynamiteformula`]\cr Formulas for stochastic channels.
#' @param sub \[`data.table`]\cr Subset of data containing
#'   the variables in the model.
#' @param u_names \[`character()`]\cr A vector of unique column names of
#'   the resulting matrix.
#' @noRd
full_model.matrix_predict <- function(dformula, sub, u_names) {
  model_matrices <- lapply(dformula, function(x) {
    model_matrices_type <- list()
    for (type in c("fixed", "varying", "random")) {
      type_formula <- get_type_formula(x, type)
      if (!is.null(type_formula)) {
        model_matrices_type[[type]] <- stats::model.matrix.lm(
          type_formula,
          data = sub,
          na.action = na.pass
        )
      }
    }
    do.call(cbind, model_matrices_type)
  })
  model_matrix <- do.call(cbind, model_matrices)
  model_matrix[, u_names, drop = FALSE]
}

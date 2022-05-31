#' Create code blocks for the Stan model
#'
#' @param dformula A `dynamiteformula` defining the model
#' @param ... TODO
#'
#' @export
create_blocks <- function(dformula, ...) {
  UseMethod("create_blocks")
}

#' @export
create_blocks.default <- function(dformula, indent = 2L, vars, ...) {
  idt <- indenter_(indent)
  functions <- create_functions(dformula, idt, vars)
  data <- create_data(dformula, idt, vars)
  transformed_data <- create_transformed_data(dformula, idt, vars)
  parameters <- create_parameters(dformula, idt, vars)
  transformed_parameters <- create_transformed_parameters(dformula, idt, vars)
  model <- create_model(dformula, idt, vars)
  generated_quantities <- create_generated_quantities(dformula, idt, vars)
  # combine above text blocks
  paste_rows(functions, data, transformed_data, parameters,
    transformed_parameters, model, generated_quantities,
    .parse = FALSE
  )
}

#' Create the 'Functions' block of the Stan model code
#'
#' @param dformula A `dynamiteformula` defining the model
#' @param idt An indeter function created by [indenter_()]
#' @param vars `model_vars` component of [convert_data()] output
#'
#' @noRd
create_functions <- function(dformula, idt, vars) {
  NULL
}

#' @describeIn create_function Create the 'Data' block of the Stan model code
#' @noRd
create_data <- function(dformula, idt, vars) {
  has_splines <- any(unlist(lapply(vars, "[[", "has_varying")))
  mtext <- paste_rows(
    "int<lower=1> T; // number of time points",
    "int<lower=1> N; // number of individuals",
    "int<lower=1> K; // total number of covariates across all channels",
    "matrix[N, K] X[T]; // covariates as an array of N x K matrices",
    "row_vector[K] X_m; // Means of all covariates at first time point",
    onlyif(has_splines, "int<lower=0> D; // number of B-splines"),
    onlyif(has_splines, "matrix[D, T] Bs; // B-spline basis matrix"),
    .indent = idt(1),
    .parse = FALSE
  )
  datatext <- character(length(dformula))
  for (i in seq_along(dformula)) {
    channel <- vars[[i]]
    y <- channel$resp
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = y, idt = idt), channel)
    datatext[i] <- lines_wrap("data", family, line_args)
  }
  paste_rows("data {", mtext, datatext, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Transformed Data'
#'   block of the Stan model code
#' @noRd
create_transformed_data <- function(dformula, idt, vars) {
  tr_data <- character(length(dformula))
  for (i in seq_along(dformula)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    tr_data[i] <- lines_wrap("transformed_data", family, line_args)
  }
  paste_rows("transformed data {", tr_data, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Parameters'
#'   block of the Stan model code
#' @noRd
create_parameters <- function(dformula, idt, vars) {
  splinetext <- ""
  if (!is.null(spline_defs <- attr(dformula, "splines"))) {
    splinetext <- paste_rows(
      "// Spline parameters",
      onlyif(spline_defs$shrinkage,
             "vector<lower=0>[D - 1] lambda; // shrinkage parameter"),
      .indent = idt(c(1, 1))
    )
    # TODO handle centered case where spline is not defined but user inserts varying(.) terms
  }
  pars <- character(length(dformula))
  for (i in seq_along(dformula)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    pars[i] <- lines_wrap("parameters", family, line_args)
  }
  paste_rows("parameters {", splinetext, pars, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Transformed Parameters'
#'   block of the Stan model code
#' @noRd
create_transformed_parameters <- function(dformula, idt, vars) {
  spline_defs <- attr(dformula, "splines")
  tr_pars <- character(length(dformula))
  for (i in seq_along(dformula)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    tr_pars[i] <- lines_wrap("transformed_parameters", family, line_args)
  }
  paste_rows("transformed parameters {", tr_pars, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Model' block of the Stan model code
#' @noRd
create_model <- function(dformula, idt, vars) {
  mod <- character(length(dformula))
  for (i in seq_along(dformula)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    mod[i] <- lines_wrap("model", family, line_args)
  }
  paste_rows("model {", mod, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Generated Quantities'
#'   block of the Stan model code
#' @noRd
create_generated_quantities <- function(dformula, idt, vars) {
  gen <- character(length(dformula))
  for (i in seq_along(dformula)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    gen[i] <- lines_wrap("generated_quantities", family, line_args)
  }
  if (any(nzchar(gen))) {
    paste_rows("generated quantities {", gen, "}", .parse = FALSE)
  } else {
    NULL
  }
}

#' Create code blocks for the Stan model
#'
#' @param becauseformula defining model components
#' @export
create_blocks <- function(formula, ...) {
  UseMethod("create_blocks")
}

#' @export
create_blocks.default <- function(formula, indent = 2L, vars, ...) {
  idt <- indenter_(indent)
  functions <- create_functions(formula, idt, vars)
  data <- create_data(formula, idt, vars)
  transformed_data <- create_transformed_data(formula, idt, vars)
  parameters <- create_parameters(formula, idt, vars)
  transformed_parameters <- create_transformed_parameters(formula, idt, vars)
  model <- create_model(formula, idt, vars)
  generated_quantities <- create_generated_quantities(formula, idt, vars)
  # combine above text blocks
  paste_rows(functions, data, transformed_data, parameters,
    transformed_parameters, model, generated_quantities,
    .parse = FALSE
  )
}
#'
#' @export
create_functions <- function(formula, idt, vars) {
  NULL
}
#'
#' @export
create_data <- function(formula, idt, vars) {
  has_splines <- any(unlist(lapply(vars, "[[", "has_varying")))
  mtext <- paste_rows(
    "int<lower=1> T; // number of time points",
    "int<lower=1> N; // number of individuals",
    "int<lower=1> K; // total number of covariates across all channels",
    "matrix[N, K] X[T]; // all covariates as an array of N x K matrices",
    onlyif(has_splines, "int<lower=0> D; // number of B-splines"),
    onlyif(has_splines, "matrix[D, T] Bs; // B-spline basis matrix"),
    .indent = idt(1),
    .parse = FALSE
  )

  # loop over channels
  datatext <- character(length(formula))
  for (i in seq_along(formula)) {
    channel <- vars[[i]]
    y <- channel$resp
    line_args <- c(list(y = y, idt = idt), channel)
    datatext[i] <- lines_wrap("data", formula[[i]], line_args)
  }
  paste_rows("data {", mtext, datatext, "}", .parse = FALSE)
}

#'
#' @export
create_transformed_data <- function(formula, idt, vars) {
  tr_data <- character(length(formula))
  for (i in seq_along(formula)) {
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    tr_data[i] <- lines_wrap("transformed_data", formula[[i]], line_args)
  }
  paste_rows("transformed data {", tr_data, "}", .parse = FALSE)
}

#'
#' @export
create_parameters <- function(formula, idt, vars) {
  splinetext <- ""
  if (!is.null(spline_defs <- attr(formula, "splines"))) {
    splinetext <- paste_rows(
      "// Spline parameters",
      onlyif(spline_defs$shrinkage,
             "vector<lower=0>[D - 1] lambda; // shrinkage parameter"),
      .indent = idt(c(1, 1))
    )
    # TODO handle centered case where spline is not defined but user inserts varying(.) terms
  }
  pars <- character(length(formula))
  for (i in seq_along(formula)) {
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    pars[i] <- lines_wrap("parameters", formula[[i]], line_args)
  }
  paste_rows("parameters {", splinetext, pars, "}", .parse = FALSE)
}

#'
#' @export
create_transformed_parameters <- function(formula, idt, vars) {
  spline_defs <- attr(formula, "splines")
  tr_pars <- character(length(formula))
  for (i in seq_along(formula)) {
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    tr_pars[i] <- lines_wrap("transformed_parameters", formula[[i]], line_args)
  }
  paste_rows("transformed parameters {", tr_pars, "}", .parse = FALSE)
}

#'
#' @export
create_model <- function(formula, idt, vars) {
  # TODO: Without global shrinkage prior it probably makes sense to use user-defined prior for tau
  # With lambda&tau, need more testing if this is fine or do we need to support other forms
  # e.g. as in https://arxiv.org/abs/1611.01310 and https://www.mdpi.com/2225-1146/8/2/20
  # priors <- character(0)

  # if (!is.null(spline_defs <- attr(formula, "splines"))) {
  #     if (spline_defs$shrinkage) {
  #         priors <- paste0(idt(1), "lambda ~ std_normal();  // prior for shrinkage terms")
  #     }
  # }
  spline_defs <- attr(formula, "splines")
  mod <- character(length(formula))
  for (i in seq_along(formula)) {
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    mod[i] <- lines_wrap("model", formula[[i]], line_args)
  }
  paste_rows("model {", mod, "}", .parse = FALSE)
}

#'
#' @export
create_generated_quantities <- function(formula, idt, vars) {
  gen <- character(length(formula))
  for (i in seq_along(formula)) {
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    gen[i] <- lines_wrap("generated_quantities", formula[[i]], line_args)
  }
  if (any(nzchar(gen))) {
    paste_rows("generated quantities {", gen, "}", .parse = FALSE)
  } else {
    NULL
  }
}

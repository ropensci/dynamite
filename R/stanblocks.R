#' Create Code Blocks for the Stan Model
#'
#' @param dformula A `dynamiteformula` defining the model
#' @param ... Additional arguments for block generation methods.
#' @noRd
create_blocks <- function(dformula, ...) {
  UseMethod("create_blocks")
}

#' Default Stan Blocks
#'
#' @inheritParams create_blocks
#' @param indent \[`integer(1)`] How many units of indentation to use for the
#'   code generation. One unit is equal to one space.
#' @param vars The `model_vars` component of [prepare_stan_data()] output.
#' @param ... Not used.
#' @noRd
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
  paste_rows(
    functions,
    data,
    transformed_data,
    parameters,
    transformed_parameters,
    model,
    generated_quantities,
    .parse = FALSE
  )
}

#' Create the 'Functions' block of the Stan model code
#'
#' @param dformula A `dynamiteformula` defining the model
#' @param idt An indeter function created by [indenter_()]
#' @param vars The `model_vars` component of [prepare_stan_data()] output.
#' @noRd
create_functions <- function(dformula, idt, vars) {
  NULL
}

#' @describeIn create_function Create The 'Data' Block of the Stan Model Code
#' @noRd
create_data <- function(dformula, idt, vars) {
  has_splines <- any(unlist(lapply(vars, "[[", "has_varying"))) ||
    any(unlist(lapply(vars, "[[", "has_varying_intercept")))
  mtext <- paste_rows(
    "int<lower=1> T; // number of time points",
    "int<lower=1> N; // number of individuals",
    "int<lower=0> K; // total number of covariates across all channels",
    "matrix[N, K] X[T]; // covariates as an array of N x K matrices",
    "row_vector[K] X_m; // Means of all covariates at first time point",
    onlyif(has_splines, "int<lower=1> D; // number of B-splines"),
    onlyif(has_splines, "matrix[D, T] Bs; // B-spline basis matrix"),
    "int<lower=0> M; // number of channels with random intercept",
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
#'   Block of the Stan Model Code
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
#'   Block of the Stan Model Code
#' @noRd
create_parameters <- function(dformula, idt, vars) {
  splinetext <- ""
  if (!is.null(spline_defs <- attr(dformula, "splines"))) {
    splinetext <- paste_rows(
      "// Spline parameters",
      onlyif(
        spline_defs$shrinkage,
        "vector<lower=0>[D - 1] lambda; // shrinkage parameter"
      ),
      .indent = idt(c(1, 1))
    )
  }
  randomtext <- ""
  has_nu <- any(unlist(lapply(dformula, "[[", "has_random_intercept")))
  if (has_nu) {
    randomtext <- paste_rows(
      "// Random intercepts",
      onlyif(attr(dformula, "correlated_nu"),
        "cholesky_factor_corr[M] L; // Cholesky for correlated intercepts"),
      "matrix[M, N] nu_raw;",
      .indent = idt(c(1, 1, 1))
    )
  }
  pars <- character(length(dformula))
  for (i in seq_along(dformula)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    pars[i] <- lines_wrap("parameters", family, line_args)
  }
  paste_rows("parameters {", splinetext, randomtext, pars, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Transformed Parameters'
#'   Block of the Stan Model Code
#' @noRd
create_transformed_parameters <- function(dformula, idt, vars) {
  randomtext <- ""
  has_nu <- any(unlist(lapply(dformula, "[[", "has_random_intercept")))
  if (has_nu) {
    if (attr(dformula, "correlated_nu")) {
      randomtext <- paste_rows(
        "matrix[N, M] nu = (L * nu_raw)';",
        .indent = idt(1))
      for (i in seq_along(dformula)) {
        y <- dformula[[i]]$response
        randomtext <- paste_rows(
          randomtext,
          glue::glue("vector[N] nu_{y} = sigma_nu_{y} * nu[, {i}];"),
          .indent = idt(c(0, 1))
        )
      }
    } else {
      for (i in seq_along(dformula)) {
        y <- dformula[[i]]$response
        randomtext <- paste_rows(randomtext,
          glue::glue("vector[N] nu_{y} = sigma_nu_{y} * nu_raw[, {i}];"))
      }
    }
  }
  tr_pars <- character(length(dformula))
  for (i in seq_along(dformula)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    tr_pars[i] <- lines_wrap("transformed_parameters", family, line_args)
  }
  paste_rows("transformed parameters {", randomtext, tr_pars, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Model' Block of the Stan Model Code
#' @noRd
create_model <- function(dformula, idt, vars) {
  splinetext <- ""
  if (!is.null(spline_defs <- attr(dformula, "splines"))) {
    paste_rows(
      splinetext <- onlyif(
        spline_defs$shrinkage,
        "lambda ~ std_normal();" #TODO user-defined prior
      ), .indent = idt(1))
  }
  randomtext <- ""
  has_nu <- any(unlist(lapply(dformula, "[[", "has_random_intercept")))
  if (has_nu) {
    randomtext <- paste_rows(
      "to_vector(nu_raw) ~ std_normal();",
      onlyif(attr(dformula, "correlated_nu"), "L ~ lkj_corr_cholesky(1);"), #TODO user-defined prior
      .indent = idt(c(1, 1))
    )
  }
  mod <- character(length(dformula))
  for (i in seq_along(dformula)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    mod[i] <- lines_wrap("model", family, line_args)
  }
  paste_rows("model {", splinetext, randomtext, mod, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Generated Quantities'
#'   Block of the Stan Model Code
#' @noRd
create_generated_quantities <- function(dformula, idt, vars) {
  gen <- ""
  has_nu <- unlist(lapply(dformula, "[[", "has_random_intercept"))
  M <- sum(has_nu)
  if (M > 0 && attr(dformula, "correlated_nu")) {
    # evaluate number of corrs to avoid Stan warning about integer division
    gen <- paste_rows(
      "corr_matrix[M] Sigma_nu = multiply_lower_tri_self_transpose(L);",
      "vector<lower=-1,upper=1>[{M*(M-1)/2}] nu_corr;",
      "for (k in 1:M) {{",
      "  for (j in 1:(k - 1)) {{",
      "    nu_corr[choose(k - 1, 2) + j] = Sigma_nu[j, k];",
      "  }}",
      "}}", .indent = idt(c(1, 1, 1, 2, 3, 2, 1))
      )
  }
  if (any(nzchar(gen))) {
    paste_rows("generated quantities {", gen, "}", .parse = FALSE)
  } else {
   NULL
  }
  # uncomment if needed in the future
  #gen <- character(length(dformula))
  #for (i in seq_along(dformula)) {
  #  family <- dformula[[i]]$family$name
  #  line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
  #  gen[i] <- lines_wrap("generated_quantities", family, line_args)
  #}
}

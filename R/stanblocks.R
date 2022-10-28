#' Create Code Blocks for the Stan Model
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula.
#' @param ... Additional arguments for block generation methods.
#' @noRd
create_blocks <- function(dformula, ...) {
  UseMethod("create_blocks")
}

#' Default Stan Blocks
#'
#' @inheritParams create_blocks
#' @inheritParams dynamite
#' @param indent \[`integer(1)`] How many units of indentation to use for the
#'   code generation. One unit is equal to one space.
#' @param vars \[`list()`]\cr The `model_vars` component of
#'   [prepare_stan_input()] output.
#' @param ... Not used.
#' @noRd
create_blocks.default <- function(dformula, indent = 2L, vars, backend, ...) {
  idt <- indenter_(indent)
  paste_rows(
    create_functions(dformula, idt, vars),
    create_data(dformula, idt, vars),
    create_transformed_data(dformula, idt, vars),
    create_parameters(dformula, idt, vars),
    create_transformed_parameters(dformula, idt, vars),
    create_model(dformula, idt, vars, backend),
    create_generated_quantities(dformula, idt, vars),
    .parse = FALSE
  )
}

#' Create the 'Functions' Block of the Stan Model Code
#'
#' @inheritParams create_blocks.default
#' @param idt \[`function`] An indentation function created by [indenter_()]
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
  n <- length(dformula)
  declarations <- character(n)
  statements <- character(n)
  for (i in seq_len(n)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    tr_data <- lines_wrap("transformed_data", family, line_args)
    declarations[i] <- tr_data$declarations
    statements[i] <- tr_data$statements
  }
  paste_rows(
    "transformed data {",
    "  // Parameters for vectorized priors",
    declarations,
    statements,
    "}",
    .parse = FALSE
  )
}

#' @describeIn create_function Create the 'Parameters'
#'   Block of the Stan Model Code
#' @noRd
create_parameters <- function(dformula, idt, vars) {
  spline_defs <- attr(dformula, "splines")
  splinetext <- ifelse_(
    is.null(spline_defs),
    "",
    paste_rows(
      onlyif(
        spline_defs$shrinkage,
        "vector<lower=0>[D - 1] xi; // Common shrinkage for splines"
      ),
      .indent = idt(1)
    )
  )
  randomtext <- ifelse_(
    identical(length(attr(dformula, "random")$responses), 0L),
    "",
    paste_rows(
      "// Random intercepts",
      onlyif(
        attr(dformula, "random")$correlated,
        "cholesky_factor_corr[M] L; // Cholesky for correlated intercepts"
      ),
      ifelse_(
        attr(dformula, "random")$noncentered,
        "matrix[N, M] nu_raw;", "vector[M] nu_raw[N];"
      ),
      .indent = idt(c(1, 1, 1))
    )
  )
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
  nus <- attr(dformula, "random")$responses
  M <- length(nus)
  if (M > 0) {
    if (attr(dformula, "random")$noncentered) {
      randomtext <- ifelse_(
        attr(dformula, "random")$correlated,
        paste_rows(
          "matrix[N, M] nu = nu_raw * L';",
          glue::glue("vector[N] nu_{nus} = sigma_nu_{nus} * nu[, {1:M}];"),
          .indent = idt(1)
        ),
        paste_rows(
          glue::glue("vector[N] nu_{nus} = sigma_nu_{nus} * nu_raw[, {1:M}];"),
          .indent = idt(1)
        )
      )
    } else {
      randomtext <-
        paste_rows(
          glue::glue("vector[N] nu_{nus} = to_vector(nu_raw[, {1:M}]);"),
          .indent = idt(1)
        )
    }
  }
  n <- length(dformula)
  declarations <- character(n)
  statements <- character(n)
  for (i in seq_len(n)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
    tr_pars <- lines_wrap("transformed_parameters", family, line_args)
    declarations[i] <- tr_pars$declarations
    statements[i] <- tr_pars$statements
  }
  paste_rows(
    "transformed parameters {",
    randomtext,
    declarations,
    statements,
    "}",
    .parse = FALSE
  )
}

#' @describeIn create_function Create the 'Model' Block of the Stan Model Code
#' @noRd
create_model <- function(dformula, idt, vars, backend) {
  spline_defs <- attr(dformula, "splines")
  splinetext <- ""
  if (!is.null(spline_defs) && spline_defs$shrinkage) {
    xi_prior <- attr(vars, "common_priors")
    xi_prior <- xi_prior[xi_prior$parameter == "xi", "prior"]
    splinetext <- paste_rows("xi ~ {xi_prior};", .indent = idt(1))
  }
  randomtext <- ""
  has_nu <- length(attr(dformula, "random")$responses) > 0
  if (has_nu) {
    if (attr(dformula, "random")$correlated) {
      L_prior <- attr(vars, "common_priors")
      L_prior <- L_prior[L_prior$parameter == "L", "prior"]
      if (attr(dformula, "random")$noncentered) {
        randomtext <- paste_rows(
          "to_vector(nu_raw) ~ std_normal();",
          "L ~ {L_prior};",
          .indent = idt(c(1, 1))
        )
      } else {
        randomtext <- paste_rows(
          "nu_raw ~ multi_normal_cholesky(0, diag_pre_multiply(sigma_nu, L));",
          onlyif(attr(dformula, "random")$correlated, "L ~ {L_prior};"),
          .indent = idt(c(1, 1))
        )
      }
    } else {
      if (attr(dformula, "random")$noncentered) {
        randomtext <- paste_rows(
          "to_vector(nu_raw) ~ std_normal();",
          .indent = idt(1)
        )
      } else {
        nus <- attr(dformula, "random")$responses
        M <- length(nus)
        randomtext <- paste_rows(
          glue::glue("nu_raw[, {1:M}] ~ normal(0, sigma_nu_{nus});"),
          .indent = idt(1)
        )
      }
    }
  }
  mod <- character(length(dformula))
  for (i in seq_along(dformula)) {
    family <- dformula[[i]]$family$name
    line_args <- c(
      list(y = vars[[i]]$resp, idt = idt, backend = backend),
      vars[[i]]
    )
    mod[i] <- lines_wrap("model", family, line_args)
  }
  paste_rows("model {", splinetext, randomtext, mod, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Generated Quantities'
#'   Block of the Stan Model Code
#' @noRd
create_generated_quantities <- function(dformula, idt, vars) {
  gen <- ""
  M <- length(attr(dformula, "random")$responses)
  if (M > 0 && attr(dformula, "random")$correlated) {
    # evaluate number of corrs to avoid Stan warning about integer division
    gen <- paste_rows(
      "corr_matrix[M] corr_matrix_nu = multiply_lower_tri_self_transpose(L);",
      "vector<lower=-1,upper=1>[{(M * (M - 1L)) %/% 2L}] corr_nu;",
      "for (k in 1:M) {{",
      "for (j in 1:(k - 1)) {{",
      "corr_nu[choose(k - 1, 2) + j] = corr_matrix_nu[j, k];",
      "}}",
      "}}",
      .indent = idt(c(1, 1, 1, 2, 3, 2, 1))
    )
  }
  if (any(nzchar(gen))) {
    paste_rows("generated quantities {", gen, "}", .parse = FALSE)
  } else {
    NULL
  }
  # uncomment if needed in the future
  # gen <- character(length(dformula))
  # for (i in seq_along(dformula)) {
  #  family <- dformula[[i]]$family$name
  #  line_args <- c(list(y = vars[[i]]$resp, idt = idt), vars[[i]])
  #  gen[i] <- lines_wrap("generated_quantities", family, line_args)
  # }
}

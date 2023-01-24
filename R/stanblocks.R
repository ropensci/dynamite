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
    any(unlist(lapply(vars, "[[", "has_varying_intercept"))) ||
    any(unlist(lapply(vars, "[[", "has_lfactor")))
  mtext <- paste_rows(
    "int<lower=1> T; // number of time points",
    "int<lower=1> N; // number of individuals",
    "int<lower=0> K; // total number of covariates across all channels",
    "matrix[N, K] X[T]; // covariates as an array of N x K matrices",
    "row_vector[K] X_m; // Means of all covariates at first time point",
    onlyif(has_splines, "int<lower=1> D; // number of B-splines"),
    onlyif(has_splines, "matrix[D, T] Bs; // B-spline basis matrix"),
    "int<lower=0> M; // number group-level effects (including intercepts)",
    "int<lower=0> P; // number of channels with latent factor",
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

  declare_AQR <- "matrix[N, N - 1] A_qr;"
  state_AQR <- paste_rows(
    "{",
    "matrix[N, N] A = diag_matrix(rep_vector(1, N));",
    "for (i in 1:N - 1) A[N, i] = -1;",
    "A[N, N] = 0;",
    "A_qr = qr_Q(A)[ , 1:(N - 1)];",
    "}",
    .indent = idt(c(0, 2, 2, 2, 2, 1)),
    .parse = FALSE
  )
  has_lfactor <- any(unlist(lapply(vars, "[[", "has_lfactor")))
  paste_rows(
    "transformed data {",
    "// Parameters for vectorized priors",
    declarations,
    onlyif(has_lfactor, declare_AQR),
    statements,
    onlyif(has_lfactor, state_AQR),
    "}",
    .indent = idt(c(0, 1, 0, 1, 0, 1, 0)),
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
    attr(vars, "random_defs")$M > 0L,
    paste_rows(
      "// Random group-level effects",
      onlyif(
        attr(vars, "random_defs")$correlated,
        "cholesky_factor_corr[M] L_nu; // Cholesky for correlated random effects"
      ),
      "vector<lower=0>[M] sigma_nu; // standard deviations of random effects",
      "matrix[N, M] nu_raw;",
      .indent = idt(c(1, 1, 1, 1))
    ),
    ""
  )

  lfactortext <- ifelse_(
    identical(length(attr(vars, "lfactor_defs")$responses), 0L),
    "",
    paste_rows(
      "// Latent factor splines",
      onlyif(
        attr(vars, "lfactor_defs")$correlated,
        "cholesky_factor_corr[P] L_lf; // Cholesky for correlated factors"
      ),
      ifelse_(
        attr(vars, "lfactor_defs")$noncentered_psi,
        "matrix[P, D - 1] omega_raw_psi;", "matrix[P, D] omega_raw_psi;"
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
  paste_rows("parameters {", splinetext, randomtext, lfactortext, pars, "}",
    .parse = FALSE)
}

#' @describeIn create_function Create the 'Transformed Parameters'
#'   Block of the Stan Model Code
#' @noRd
create_transformed_parameters <- function(dformula, idt, vars) {

  randomtext <- ""
  M <- attr(vars, "random_defs")$M
  if (M > 0L) {
    Ks <- unlist(lapply(vars, "[[", "K_random"))
    Ks <- Ks[Ks > 0]
    y <- names(Ks)
    cKs1 <- cumsum(c(1, Ks[-length(Ks)]))
    cKs2 <- cumsum(Ks)
    if (attr(vars, "random_defs")$noncentered) {
      randomtext <- ifelse_(
        attr(vars, "random_defs")$correlated,
        paste_rows(
          "matrix[N, M] nu = nu_raw * diag_pre_multiply(sigma_nu, L_nu)';",
          glue::glue("matrix[N, {Ks}] nu_{y} = nu[, {cKs1}:{cKs2}];"),
          .indent = idt(1)
        ),
        paste_rows(
          glue::glue("matrix[N, {Ks}] nu_{y} = diag_post_multiply(nu_raw[, {cKs1}:{cKs2}], sigma_nu_{y});"),
          .indent = idt(1)
        )
      )
    } else {
      randomtext <-
        paste_rows(
          glue::glue("matrix[N, {Ks}] nu_{y} = nu_raw[, {cKs1}:{cKs2}];"),
          .indent = idt(1)
        )
    }
    randomtext <- paste_rows(
      glue::glue("vector[{Ks}] sigma_nu_{y} = sigma_nu[{cKs1}:{cKs2}];"),
      randomtext,
      .indent = idt(c(1, 0))
    )
  }

  lfactortext <- ""
  psis <- attr(vars, "lfactor_defs")$responses
  P <- length(psis)
  if (P > 0) {
    if (attr(vars, "lfactor_defs")$noncentered_psi) {
      tau_psi <- ifelse(attr(vars, "lfactor_defs")$nonzero_lambda,
        paste0("tau_psi_", psis, " * "), "")
      # create noise terms for latent factors
      lfactortext <- ifelse_(
        attr(vars, "lfactor_defs")$correlated,
        paste_rows(
          "matrix[P, D - 1] omega_psi = L_lf * omega_raw_psi;",
          glue::glue(
            "row_vector[D] omega_psi_{psis} = \\
            append_col(0, {tau_psi}omega_psi[{1:P}, ]);"),
          .indent = idt(1)
        ),
        paste_rows(
          glue::glue("row_vector[D] omega_psi_{psis} = \\
            append_col(0, {tau_psi}omega_raw_psi[{1:P}, ]);"),
          .indent = idt(1)
        )
      )
    } else {
      lfactortext <-
        paste_rows(
          glue::glue("row_vector[D] omega_psi_{psis} = \\
            omega_raw_psi[{1:P}, ];"),
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
    lfactortext,
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
    splinetext <- paste_rows("xi[1] ~ {xi_prior};", .indent = idt(1))
  }
  randomtext <- ""
  if (attr(vars, "random_defs")$M > 0L) {
    if (attr(vars, "random_defs")$correlated) {
      L_prior <- attr(vars, "common_priors")
      L_prior <- L_prior[L_prior$parameter == "L_nu", "prior"]
      randomtext <- ifelse_(
        attr(vars, "random_defs")$noncentered,
        paste_rows(
          "to_vector(nu_raw) ~ std_normal();",
          "L_nu ~ {L_prior};",
          .indent = idt(c(1, 1))
        ),
        paste_rows(
          "{{",
          "row_vector[M] nu_tmp[N];",
          "for (i in 1:N) {{",
            "nu_tmp[i] = nu_raw[i, ];",
          "}}",
          "nu_tmp ~ multi_normal_cholesky(rep_vector(0, M), diag_pre_multiply(sigma_nu, L_nu));",
          "}}",
          "L_nu ~ {L_prior};",
          .indent = idt(c(1, 2, 2, 3, 2, 2, 1, 1))
        )
      )
    } else {
      M <- attr(vars, "random_defs")$M
      Ks <- unlist(lapply(vars, "[[", "K_random"))
      y <- names(Ks[Ks > 0])
      randomtext <- ifelse_(
        attr(vars, "random_defs")$noncentered,
        paste_rows(
          "to_vector(nu_raw) ~ std_normal();",
          .indent = idt(1)
        ),
        paste_rows(
          "for (i in 1:M) {{",
            "nu_raw[, i] ~ normal(0, sigma_nu[i]);",
          "}}",
          .indent = idt(c(1, 2, 1))
        )
      )
    }
  }

  lfactortext <- ""
  psis <- attr(vars, "lfactor_defs")$responses
  P <- length(psis)
  if (P > 0L) {
    if (attr(vars, "lfactor_defs")$correlated) {
      L_prior <- attr(vars, "common_priors")
      L_prior <- L_prior[L_prior$parameter == "L_lf", "prior"]
      if (attr(vars, "lfactor_defs")$noncentered_psi) {
        lfactortext <- paste_rows(
          "to_vector(omega_raw_psi) ~ std_normal();",
          "L_lf ~ {L_prior};",
          .indent = idt(c(1, 1))
        )
      } else {
        if (any(attr(vars, "lfactor_defs")$nonzero_lambda)) {
          tau <- paste0(
            ifelse(
              attr(vars, "lfactor_defs")$nonzero_lambda,
              paste0("tau_psi_", psis),
              "1"
            ),
            collapse = ", "
          )
          lfactortext <- paste_rows(
            "L_lf ~ {L_prior};",
            "{{",
            "vector[P] tau_psi = [{tau}]';",
            "matrix[P, P] Ltau = diag_pre_multiply(tau_psi, L_lf);",
            "for (i in 2:D) {{",
            paste0(
              "omega_raw_psi[, i] ~ ",
              "multi_normal_cholesky(omega_raw_psi[, i - 1], Ltau);"
            ),
            "}}",
            "}}",
            .indent = idt(c(1, 1, 2, 2, 2, 3, 2, 1))
          )
        } else {
          lfactortext <- paste_rows(
            "L_lf ~ {L_prior};",
            "for (i in 2:D) {{",
            paste0(
              "omega_raw_psi[, i] ~ ",
              "multi_normal_cholesky(omega_raw_psi[, i - 1], L_lf);"
            ),
            "}}",
            .indent = idt(c(1, 1, 2, 1))
          )
        }
      }
    } else {
      if (attr(vars, "lfactor_defs")$noncentered_psi) {
        lfactortext <- paste_rows(
          "to_vector(omega_raw_psi) ~ std_normal();",
          .indent = idt(1)
        )
      } else {
        if (any(attr(vars, "lfactor_defs")$nonzero_lambda)) {
          psis <- attr(vars, "lfactor_defs")$responses
          P <- length(psis)
          tau <- paste0(ifelse(attr(vars, "lfactor_defs")$nonzero_lambda,
            paste0("tau_psi_", psis), "1"), collapse = ", ")
          lfactortext <- paste_rows(
            "{{",
            "vector[P] tau_psi = [{tau}]';",
            "for (i in 2:D) {{",
            "omega_raw_psi[, i] ~ normal(omega_raw_psi[, i - 1], tau_psi);",
            "}}",
            "}}",
            .indent = idt(c(1, 2, 2, 3, 2, 1))
          )
        } else {
          lfactortext <- paste_rows(
            "for (i in 2:D) {{",
            "omega_raw_psi[, i] ~ normal(omega_raw_psi[, i - 1], 1);",
            "}}",
            .indent = idt(c(1, 2, 1))
          )
        }
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
  paste_rows(
    "model {",
    splinetext,
    randomtext,
    lfactortext,
    mod,
    "}",
    .parse = FALSE
  )
}

#' @describeIn create_function Create the 'Generated Quantities'
#'   Block of the Stan Model Code
#' @noRd
create_generated_quantities <- function(dformula, idt, vars) {
  gen_nu <- ""
  M <- attr(vars, "random_defs")$M
  if (M > 1 && attr(vars, "random_defs")$correlated) {
    # evaluate number of corrs to avoid Stan warning about integer division
    gen_nu <- paste_rows(
      "corr_matrix[M] corr_matrix_nu = multiply_lower_tri_self_transpose(L_nu);",
      "vector<lower=-1,upper=1>[{(M * (M - 1L)) %/% 2L}] corr_nu;",
      "for (k in 1:M) {{",
      "for (j in 1:(k - 1)) {{",
      "corr_nu[choose(k - 1, 2) + j] = corr_matrix_nu[j, k];",
      "}}",
      "}}",
      .indent = idt(c(1, 1, 1, 2, 3, 2, 1))
    )
  }
  gen_psi <- ""
  psis <- attr(vars, "lfactor_defs")$responses
  P <- attr(vars, "lfactor_defs")$P
  if (P > 0 && attr(vars, "lfactor_defs")$correlated) {
    # evaluate number of corrs to avoid Stan warning about integer division
    tau <- paste0(ifelse(attr(vars, "lfactor_defs")$nonzero_lambda,
      paste0("tau_psi_", psis), "1"), collapse = ", ")
    gen_psi <- paste_rows(
      "vector[P] tau_psi = [{tau}]';",
      "corr_matrix[P] corr_matrix_psi =
       multiply_lower_tri_self_transpose(diag_pre_multiply(tau_psi, L_lf));",
      "vector<lower=-1,upper=1>[{(P * (P - 1L)) %/% 2L}] corr_psi;",
      "for (k in 1:P) {{",
      "for (j in 1:(k - 1)) {{",
      "corr_psi[choose(k - 1, 2) + j] = corr_matrix_psi[j, k];",
      "}}",
      "}}",
      .indent = idt(c(1, 1, 1, 1, 2, 3, 2, 1))
    )
  }
  if (any(nzchar(gen_nu)) || any(nzchar(gen_psi))) {
    paste_rows("generated quantities {", gen_nu, gen_psi, "}", .parse = FALSE)
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

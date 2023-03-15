#' Create Stan Blocks
#'
#' @param indent \[`integer(1)`] How many units of indentation to use for the
#'   code generation. One unit is equal to one space.
#' @param cvars \[`list()`]\cr The `channel_vars` component of
#' @param cgvars \[`list()`]\cr The `channel_group_vars` component of
#'   [prepare_stan_input()] output.
#' @param cg \[`integer()`]\cr The `"channel_groups"` attribute of the `dformula`
#'   for stochastic channels
#' @param ... Not used.
#' @noRd
create_blocks <- function(indent = 2L, cvars, cgvars, cg, backend) {
  idt <- indenter_(indent)
  paste_rows(
    create_functions(idt, cvars, cgvars, cg),
    create_data(idt, cvars, cgvars, cg, backend),
    create_transformed_data(idt, cvars, cgvars, cg, backend),
    create_parameters(idt, cvars, cgvars, cg, backend),
    create_transformed_parameters(idt, cvars, cgvars, cg, backend),
    create_model(idt, cvars, cgvars, cg, backend),
    create_generated_quantities(idt, cvars, cgvars, cg, backend),
    .parse = FALSE
  )
}

#' Create the 'Functions' Block of the Stan Model Code
#'
#' @inheritParams create_blocks
#' @param idt \[`function`]\cr An indentation function created by [indenter_()]
#' @noRd
create_functions <- function(idt, cvars, cgvars, cg) {
  functions_psi <- ""
  psis <- attr(cvars, "lfactor_def")$responses
  P <- attr(cvars, "lfactor_def")$P
  # From Stan forums https://tinyurl.com/2spznmyv
  if (P > 0) {
    functions_psi <- paste_rows(
      "vector create_Q(int N) {{",
      "vector[2 * N] Q;",
      "for (i in 1:N) {{",
      "Q[i] = -sqrt((N - i)/(N - i + 1.0));",
      "Q[i + N] = inv_sqrt((N - i) * (N - i + 1));",
      "}}",
      "return Q;",
      "}}",
      "vector sum_to_zero(vector x_raw, vector Q) {{",
      "int N = num_elements(x_raw) + 1;",
      "vector[N] x;",
      "real x_aux = 0;",
      "for (i in 1:(N - 1)){{",
      "x[i] = x_aux + x_raw[i] * Q[i];",
      "x_aux = x_aux + x_raw[i] * Q[i + N];",
      "}}",
      "x[N] = x_aux;",
      "return x;",
      "}}",
      .indent = idt(c(1, 2, 2, 3, 3, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, 2, 2, 2, 1))
    )
  }
  paste_rows(
    "functions {",
    functions_psi,
    "}",
    .parse = FALSE
  )
}

#' @describeIn create_function Create The 'Data' Block of the Stan Model Code
#' @noRd
create_data <- function(idt, cvars, cgvars, cg, backend) {
  has_splines <- any(vapply(cvars, "[[", logical(1L), "has_varying")) ||
    any(vapply(cvars, "[[", logical(1L), "has_varying_intercept")) ||
    any(vapply(cvars, "[[", logical(1L), "has_lfactor"))
  K <- sum(vapply(cvars, "[[", integer(1), "K"))
  M <- attr(cvars, "random_def")$M
  P <- attr(cvars, "lfactor_def")$P
  init_text <- paste_rows(
    "int<lower=1> T; // number of time points",
    "int<lower=1> N; // number of individuals",
    onlyif(
      K > 0,
      "int<lower=0> K; // total number of covariates across all channels"
    ),
    onlyif(
      K > 0,
      stan_array(
        backend,
        "matrix",
        "X",
        "T",
        dims = "N, K",
        comment = "covariates as an array of N x K matrices"
      )
    ),
    onlyif(
      K > 0,
      "row_vector[K] X_m; // Means of all covariates at first time point"
    ),
    onlyif(has_splines, "int<lower=1> D; // number of B-splines"),
    onlyif(has_splines, "matrix[D, T] Bs; // B-spline basis matrix"),
    onlyif(
      M > 0,
      "int<lower=0> M; // number group-level effects (including intercepts)"
    ),
    onlyif(
      P > 0,
      "int<lower=0> P; // number of channels with latent factor"
    ),
    .indent = idt(1),
    .parse = FALSE
  )
  n_cg <- n_unique(cg)
  data_text <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    if (is_multivariate(cgvars[[i]]$family)) {
      if (has_univariate(cgvars[[i]]$family)) {
        cgvars[[i]]$default <- vapply(
          cvars[cg_idx],
          function(x) {
            lines_wrap("data", "default", x, idt, backend)
          },
          character(1L)
        )
      }
      data_text[i] <- lines_wrap(
        "data",  cgvars[[i]]$family, cgvars[[i]], idt, backend
      )
    } else {
      j <- cg_idx[1L]
      cvars[[j]]$default <- lines_wrap(
        "data", "default", cvars[[j]], idt, backend
      )
      data_text[i] <- lines_wrap(
        "data", cvars[[j]]$family, cvars[[j]], idt, backend
      )
    }
  }
  paste_rows("data {", init_text, data_text, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Transformed Data'
#'   Block of the Stan Model Code
#' @noRd
create_transformed_data <- function(idt, cvars, cgvars, cg, backend) {
  n_cg <- n_unique(cg)
  declarations <- character(n_cg)
  statements <- character(n_cg)
  tr_data <- list()
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    if (is_multivariate(cgvars[[i]]$family)) {
      cgvars[[i]]$default <- lapply(
        cvars[cg_idx],
        function(x) {
          lines_wrap("transformed_data", "default", x, idt, backend)
        }
      )
      tr_data <- lines_wrap(
        "transformed_data", cgvars[[i]]$family, cgvars[[i]], idt, backend
      )
    } else {
      j <- cg_idx[1L]
      cvars[[j]]$default <- lines_wrap(
        "transformed_data", "default", cvars[[j]], idt, backend
      )
      tr_data <- lines_wrap(
        "transformed_data", cvars[[j]]$family, cvars[[j]], idt, backend
      )
    }
    declarations[i] <- tr_data$declarations
    statements[i] <- tr_data$statements
  }
  declare_QR <- "vector[2 * N] QR_Q = create_Q(N);"
  has_lfactor <- any(vapply(cvars, "[[", logical(1L), "has_lfactor"))
  paste_rows(
    "transformed data {",
    declarations,
    onlyif(has_lfactor, declare_QR),
    statements,
    "}",
    .indent = idt(c(0, 0, 1, 0, 0)),
    .parse = FALSE
  )
}

#' @describeIn create_function Create the 'Parameters'
#'   Block of the Stan Model Code
#' @noRd
create_parameters <- function(idt, cvars, cgvars, cg, backend) {
  spline_def <- attr(cvars, "spline_def")
  spline_text <- ifelse_(
    is.null(spline_def),
    "",
    paste_rows(
      onlyif(
        spline_def$shrinkage,
        "vector<lower=0>[D - 1] xi; // Common shrinkage for splines"
      ),
      .indent = idt(1)
    )
  )
  random_text <- ifelse_(
    attr(cvars, "random_def")$M > 0L,
    paste_rows(
      "// Random group-level effects",
      onlyif(
        attr(cvars, "random_def")$correlated,
        paste0(
          "cholesky_factor_corr[M] L_nu; ",
          "// Cholesky for correlated random effects"
        )
      ),
      "vector<lower=0>[M] sigma_nu; // standard deviations of random effects",
      "matrix[N, M] nu_raw;",
      .indent = idt(c(1, 1, 1, 1))
    ),
    ""
  )
  lfactor_text <- ifelse_(
    identical(length(attr(cvars, "lfactor_def")$responses), 0L),
    "",
    paste_rows(
      "// Latent factor splines",
      onlyif(
        attr(cvars, "lfactor_def")$correlated,
        "cholesky_factor_corr[P] L_lf; // Cholesky for correlated factors"
      ),
      "matrix[P, D - 1] omega_raw_psi;",
      .indent = idt(c(1, 1, 1))
    )
  )
  n_cg <- n_unique(cg)
  parameters_text <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    if (is_multivariate(cgvars[[i]]$family)) {
      univariate <- ""
      for (j in cg_idx) {
        cvars[[j]]$default <- lines_wrap(
          "parameters", "default", cvars[[j]], idt, backend
        )
        univariate <- ifelse_(
          has_univariate(cgvars[[i]]$family),
          paste_rows(
            univariate,
            lines_wrap(
              "parameters",
              get_univariate(cvars[[j]]$family),
              cvars[[j]],
              idt,
              backend
            )
          ),
          paste_rows(univariate, cvars[[j]]$default)
        )
      }
      cgvars[[i]]$univariate <- univariate
      parameters_text[i] <- lines_wrap(
        "parameters", cgvars[[i]]$family, cgvars[[i]], idt, backend
      )
    } else {
      j <- cg_idx[1L]
      cvars[[j]]$default <- lines_wrap(
        "parameters", "default", cvars[[j]], idt, backend
      )
      parameters_text[i] <- lines_wrap(
        "parameters", cvars[[j]]$family, cvars[[j]], idt, backend
      )
    }
  }
  paste_rows(
    "parameters {",
    spline_text,
    random_text,
    lfactor_text,
    parameters_text,
    "}",
    .parse = FALSE
  )
}

#' @describeIn create_function Create the 'Transformed Parameters'
#'   Block of the Stan Model Code
#' @noRd
create_transformed_parameters <- function(idt, cvars, cgvars, cg, backend) {
  random_text <- ""
  M <- attr(cvars, "random_def")$M
  if (M > 0L) {
    Ks <- vapply(cvars, "[[", integer(1L), "K_random")
    Ks <- Ks[Ks > 0]
    y <- names(Ks)
    cKs1 <- cumsum(c(1, Ks[-length(Ks)]))
    cKs2 <- cumsum(Ks)
    if (attr(cvars, "random_def")$noncentered) {
      random_text <- ifelse_(
        attr(cvars, "random_def")$correlated,
        paste_rows(
          "matrix[N, M] nu = nu_raw * diag_pre_multiply(sigma_nu, L_nu)';",
          "matrix[N, {Ks}] nu_{y} = nu[, {cKs1}:{cKs2}];",
          .indent = idt(1)
        ),
        paste_rows(
          paste0(
            "matrix[N, {Ks}] nu_{y} = ",
            "diag_post_multiply(nu_raw[, {cKs1}:{cKs2}], sigma_nu_{y});"
          ),
          .indent = idt(1)
        )
      )
    } else {
      random_text <-
        paste_rows(
          "matrix[N, {Ks}] nu_{y} = nu_raw[, {cKs1}:{cKs2}];",
          .indent = idt(1)
        )
    }
    random_text <- paste_rows(
      "vector[{Ks}] sigma_nu_{y} = sigma_nu[{cKs1}:{cKs2}];",
      random_text,
      .indent = idt(c(1, 0))
    )
  }
  lfactor_text <- ""
  psis <- attr(cvars, "lfactor_def")$responses
  P <- length(psis)
  if (P > 0) {
    if (attr(cvars, "lfactor_def")$noncentered_psi) {
      tau_psi <- ifelse(
        attr(cvars, "lfactor_def")$nonzero_lambda,
        paste0("tau_psi_", psis, " * "),
        ""
      )
      # create noise terms for latent factors
      lfactor_text <- ifelse_(
        attr(cvars, "lfactor_def")$correlated,
        paste_rows(
          "matrix[P, D - 1] omega_psi = L_lf * omega_raw_psi;",
          paste0(
            "row_vector[D] omega_psi_{psis} = ",
            "append_col(0, {tau_psi}omega_psi[{1:P}, ]);"
          ),
          .indent = idt(1)
        ),
        paste_rows(
          paste0(
            "row_vector[D] omega_psi_{psis} = ",
            "append_col(0, {tau_psi}omega_raw_psi[{1:P}, ]);"
          ),
          .indent = idt(1)
        )
      )
    } else {
      lfactor_text <-
        paste_rows(
          paste0(
            "row_vector[D] omega_psi_{psis} = ",
            "append_col(0, omega_raw_psi[{1:P}, ]);"
          ),
          .indent = idt(1)
        )
    }
  }
  n_cg <- n_unique(cg)
  declarations <- character(n_cg)
  statements <- character(n_cg)
  tr_pars <- list()
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    if (is_multivariate(cgvars[[i]]$family)) {
      cgvars[[i]]$default <- lapply(
        cvars[cg_idx],
        function(x) {
          lines_wrap("transformed_parameters", "default", x, idt, backend)
        }
      )
      tr_pars <- lines_wrap(
        "transformed_parameters", cgvars[[i]]$family, cgvars[[i]], idt, backend
      )
    } else {
      j <- cg_idx[1L]
      cvars[[j]]$default <- lines_wrap(
        "transformed_parameters", "default", cvars[[j]], idt, backend
      )
      tr_pars <- lines_wrap(
        "transformed_parameters", cvars[[j]]$family, cvars[[j]], idt, backend
      )
    }
    declarations[i] <- tr_pars$declarations
    statements[i] <-  tr_pars$statements
  }
  paste_rows(
    "transformed parameters {",
    random_text,
    lfactor_text,
    declarations,
    statements,
    "}",
    .parse = FALSE
  )
}

#' @describeIn create_function Create the 'Model' Block of the Stan Model Code
#' @noRd
create_model <- function(idt, cvars, cgvars, cg, backend) {
  spline_def <- attr(cvars, "spline_def")
  spline_text <- ""
  if (!is.null(spline_def) && spline_def$shrinkage) {
    xi_prior <- attr(cvars, "common_priors")
    xi_prior <- xi_prior[xi_prior$parameter == "xi", "prior"]
    spline_text <- paste_rows("xi[1] ~ {xi_prior};", .indent = idt(1))
  }
  random_text <- ""
  if (attr(cvars, "random_def")$M > 0L) {
    if (attr(cvars, "random_def")$correlated) {
      L_prior <- attr(cvars, "common_priors")
      L_prior <- L_prior[L_prior$parameter == "L_nu", "prior"]
      random_text <- ifelse_(
        attr(cvars, "random_def")$noncentered,
        paste_rows(
          "to_vector(nu_raw) ~ std_normal();",
          "L_nu ~ {L_prior};",
          .indent = idt(c(1, 1))
        ),
        paste_rows(
          "{{",
          stan_array(backend, "row_vector", "nu_tmp", "N", dims = "M"),
          "for (i in 1:N) {{",
          "nu_tmp[i] = nu_raw[i, ];",
          "}}",
          paste0(
            "nu_tmp ~ multi_normal_cholesky(rep_vector(0, M), ",
            "diag_pre_multiply(sigma_nu, L_nu));"
          ),
          "}}",
          "L_nu ~ {L_prior};",
          .indent = idt(c(1, 2, 2, 3, 2, 2, 1, 1))
        )
      )
    } else {
      M <- attr(cvars, "random_def")$M
      Ks <- vapply(cvars, "[[", integer(1L), "K_random")
      y <- names(Ks[Ks > 0])
      random_text <- ifelse_(
        attr(cvars, "random_def")$noncentered,
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
  lfactor_text <- ""
  psis <- attr(cvars, "lfactor_def")$responses
  P <- length(psis)
  if (P > 0L) {
    omega1 <- paste0("omega_raw_psi_1_", psis, collapse = ", ")
    if (attr(cvars, "lfactor_def")$correlated) {
      L_prior <- attr(cvars, "common_priors")
      L_prior <- L_prior[L_prior$parameter == "L_lf", "prior"]
      if (attr(cvars, "lfactor_def")$noncentered_psi) {
        lfactor_text <- paste_rows(
          "to_vector(omega_raw_psi) ~ std_normal();",
          "L_lf ~ {L_prior};",
          .indent = idt(c(1, 1))
        )
      } else {
        if (any(attr(cvars, "lfactor_def")$nonzero_lambda)) {
          tau <- paste0(
            ifelse(
              attr(cvars, "lfactor_def")$nonzero_lambda,
              paste0("tau_psi_", psis),
              "1"
            ),
            collapse = ", "
          )
          lfactor_text <- paste_rows(
            "L_lf ~ {L_prior};",
            "{{",
            "vector[P] tau_psi = [{tau}]';",
            "matrix[P, P] Ltau = diag_pre_multiply(tau_psi, L_lf);",
            "vector[P] omega1 = [{omega1}]';",
            "omega_raw_psi[, 2] ~ multi_normal_cholesky(omega1, Ltau);",
            "for (i in 3:(D - 1)) {{",
            paste0(
              "omega_raw_psi[, i] ~ ",
              "multi_normal_cholesky(omega_raw_psi[, i - 1], Ltau);"
            ),
            "}}",
            "}}",
            .indent = idt(c(1, 1, 2, 2, 2, 2, 2, 3, 2, 1))
          )
        } else {
          lfactor_text <- paste_rows(
            "L_lf ~ {L_prior};",
            "vector[P] omega1 = [{omega1}]';",
            "omega_raw_psi[, 2] ~ multi_normal_cholesky(omega1, L_lf);",
            "for (i in 3:(D - 1)) {{",
            paste0(
              "omega_raw_psi[, i] ~ ",
              "multi_normal_cholesky(omega_raw_psi[, i - 1], L_lf);"
            ),
            "}}",
            .indent = idt(c(1, 1, 1, 1, 2, 1))
          )
        }
      }
    } else {
      if (attr(cvars, "lfactor_def")$noncentered_psi) {
        lfactor_text <- paste_rows(
          "to_vector(omega_raw_psi) ~ std_normal();",
          .indent = idt(1)
        )
      } else {
        if (any(attr(cvars, "lfactor_def")$nonzero_lambda)) {
          tau <- paste0(
            ifelse(
              attr(cvars, "lfactor_def")$nonzero_lambda,
              paste0("tau_psi_", psis),
              "1"
            ),
            collapse = ", "
          )
          lfactor_text <- paste_rows(
            "{{",
            "vector[P] tau_psi = [{tau}]';",
            "vector[P] omega1 = [{omega1}]';",
            "omega_raw_psi[, 2] ~ normal(omega1, tau_psi);",
            "for (i in 3:(D - 1)) {{",
            "omega_raw_psi[, i] ~ normal(omega_raw_psi[, i - 1], tau_psi);",
            "}}",
            "}}",
            .indent = idt(c(1, 2, 2, 2, 2, 3, 2, 1))
          )
        } else {
          lfactor_text <- paste_rows(
            "vector[P] omega1 = [{omega1}]';",
            "omega_raw_psi[, 2] ~ normal(omega1, 1);",
            "for (i in 3:(D - 1)) {{",
            "omega_raw_psi[, i] ~ normal(omega_raw_psi[, i - 1], 1);",
            "}}",
            .indent = idt(c(1, 1, 1, 2, 1))
          )
        }
      }
    }
  }
  n_cg <- n_unique(cg)
  model_text <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    if (is_multivariate(cgvars[[i]]$family)) {
      cgvars[[i]]$backend <- backend
      model_text[i] <- lines_wrap(
        "model", cgvars[[i]]$family,
        list(cvars = cvars[cg_idx], cgvars = cgvars[[i]]),
        idt,
        backend
      )
    } else {
      j <- cg_idx[1L]
      cvars[[j]]$backend <- backend
      cvars[[j]]$priors <- do.call(prior_lines, c(cvars[[j]], idt = idt))
      cvars[[j]]$intercept <- do.call(intercept_lines, cvars[[j]])
      model_text[i] <- lines_wrap(
        "model", cvars[[j]]$family, cvars[[j]], idt, backend
      )
    }
  }
  paste_rows(
    "model {",
    spline_text,
    random_text,
    lfactor_text,
    model_text,
    "}",
    .parse = FALSE
  )
}

#' @describeIn create_function Create the 'Generated Quantities'
#'   Block of the Stan Model Code
#' @noRd
create_generated_quantities <- function(idt, cvars, cgvars, cg, backend) {
  gen_nu <- ""
  M <- attr(cvars, "random_def")$M
  if (M > 1 && attr(cvars, "random_def")$correlated) {
    # evaluate number of corrs to avoid Stan warning about integer division
    gen_nu <- paste_rows(
      paste0(
        "matrix[M, M] corr_matrix_nu = ",
        "multiply_lower_tri_self_transpose(L_nu);"
      ),
      "vector[{(M * (M - 1L)) %/% 2L}] corr_nu;",
      "for (k in 1:M) {{",
      "for (j in 1:(k - 1)) {{",
      "corr_nu[choose(k - 1, 2) + j] = corr_matrix_nu[j, k];",
      "}}",
      "}}",
      .indent = idt(c(1, 1, 1, 2, 3, 2, 1))
    )
  }
  gen_psi <- ""
  psis <- attr(cvars, "lfactor_def")$responses
  P <- attr(cvars, "lfactor_def")$P
  if (P > 0 && attr(cvars, "lfactor_def")$correlated) {
    # evaluate number of corrs to avoid Stan warning about integer division
    tau <- paste0(
      ifelse(
        attr(cvars, "lfactor_def")$nonzero_lambda,
        paste0("tau_psi_", psis),
        "1"
      ),
      collapse = ", "
    )
    gen_psi <- paste_rows(
      paste0(
        "matrix[P, P] corr_matrix_psi = ",
        "multiply_lower_tri_self_transpose(L_lf);"
      ),
      "vector[{(P * (P - 1L)) %/% 2L}] corr_psi;",
      "for (k in 1:P) {{",
      "for (j in 1:(k - 1)) {{",
      "corr_psi[choose(k - 1, 2) + j] = corr_matrix_psi[j, k];",
      "}}",
      "}}",
      .indent = idt(c(1, 1, 1, 2, 3, 2, 1))
    )
  }
  n_cg <- n_unique(cg)
  generated_quantities_text <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    j <- cg_idx[1L]
    generated_quantities_text[i] <- ifelse_(
      is_multivariate(cgvars[[i]]$family),
      lines_wrap(
        "generated_quantities", cgvars[[i]]$family, cgvars[[i]], idt, backend
      ),
      lines_wrap(
        "generated_quantities", cvars[[j]]$family, cvars[[j]], idt, backend
      )
    )
  }
  paste_rows(
    "generated quantities {",
    gen_nu,
    gen_psi,
    generated_quantities_text,
    "}",
    .parse = FALSE
  )
}

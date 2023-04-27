#' Create Stan Blocks
#'
#' @param indent \[`integer(1)`]\cr How many units of indentation to use for
#'   the code generation. One unit is equal to one space.
#' @param backend \[`character(1)`]\cr Either `"rstan"` or `"cmdstanr"`.
#' @param cg \[`integer()`]\cr The `"channel_groups"` attribute of the
#'   `dformula` for stochastic channels.
#' @param cvars \[`list()`]\cr The `channel_vars` component of
#'   [prepare_stan_input()] output.
#' @param cgvars \[`list()`]\cr The `channel_group_vars` component of
#'   [prepare_stan_input()] output.
#' @param mvars \[`list()`]\cr The `model_vars` component of
#'   [prepare_stan_input()] output.
#' @noRd
create_blocks <- function(indent = 2L, backend, cg, cvars, cgvars, mvars,
                          threading) {
  idt <- indenter_(indent)
  paste_rows(
    create_functions(idt, backend, cg, cvars, cgvars, mvars, threading),
    create_data(idt, backend, cg, cvars, cgvars, mvars, threading),
    create_transformed_data(idt, backend, cg, cvars, cgvars, mvars),
    create_parameters(idt, backend, cg, cvars, cgvars, mvars),
    create_transformed_parameters(idt, backend, cg, cvars, cgvars, mvars),
    create_model(idt, backend, cg, cvars, cgvars, mvars, threading),
    create_generated_quantities(idt, backend, cg, cvars, cgvars, mvars),
    .parse = FALSE
  )
}

#' Create the 'Functions' Block of the Stan Model Code
#'
#' @inheritParams create_blocks
#' @param idt \[`function`]\cr
#'   An indentation function created by [indenter_()].
#' @noRd
create_functions <- function(idt, backend, cg, cvars, cgvars, mvars, threading) {
  functions_psi <- ""
  psis <- mvars$lfactor_def$responses
  P <- mvars$lfactor_def$P
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

  n_cg <- n_unique(cg)
  likelihood_functions_text <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    likelihood_functions_text[i] <- create_functions_lines(
      idt, backend, cvars[cg_idx], cgvars[[i]], threading
    )
  }
  paste_rows(
    "functions {",
    functions_psi,
    likelihood_functions_text,
    "}",
    .parse = FALSE
  )
}

#' Create Functions Lines for a Distribution
#'
#' @noRd
create_functions_lines <- function(idt, backend, cvars, cgvars, threading) {
  family <- cgvars$family
  if (is_multivariate(family)) {
    lines_wrap("functions", family, idt, backend,
               list(cvars = cvars, cgvars = cgvars, threading = threading))
  } else {
    if (is_categorical(family)) {
      cvars[[1L]]$threading <- threading
    } else {
      cvars[[1L]]$default <- lines_wrap(
        "functions", "default", idt, backend,
        c(cvars[[1L]], threading = threading)
      )
    }
    lines_wrap("functions", family, idt, backend, cvars[[1L]])
  }
}
#' @describeIn create_function Create The 'Data' Block of the Stan Model Code
#' @noRd
create_data <- function(idt, backend, cg, cvars, cgvars, mvars, threading) {
  has_splines <- any(vapply(cvars, "[[", logical(1L), "has_varying")) ||
    any(vapply(cvars, "[[", logical(1L), "has_varying_intercept")) ||
    any(vapply(cvars, "[[", logical(1L), "has_lfactor"))
  K <- mvars$K
  M <- mvars$random_def$M
  P <- mvars$lfactor_def$P
  init_text <- paste_rows(
    "int<lower=1> T; // number of time points",
    "int<lower=1> N; // number of individuals",
    onlyif(
      K > 0L,
      "int<lower=0> K; // total number of covariates across all channels"
    ),
    onlyif(
      K > 0L,
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
      K > 0L,
      "row_vector[K] X_m; // Means of all covariates at first time point"
    ),
    onlyif(has_splines, "int<lower=1> D; // number of B-splines"),
    onlyif(has_splines, "matrix[D, T] Bs; // B-spline basis matrix"),
    onlyif(
      M > 0L,
      "int<lower=0> M; // number group-level effects (including intercepts)"
    ),
    onlyif(
      P > 0L,
      "int<lower=0> P; // number of channels with latent factor"
    ),
    onlyif(threading, "int<lower=1> grainsize;"),
    .indent = idt(1),
    .parse = FALSE
  )
  n_cg <- n_unique(cg)
  data_text <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    data_text[i] <- create_data_lines(idt, backend, cvars[cg_idx], cgvars[[i]])
  }
  paste_rows("data {", init_text, data_text, "}", .parse = FALSE)
}

#' Create Data Lines for a Distribution
#'
#' @noRd
create_data_lines <- function(idt, backend, cvars, cgvars) {
  family <- cgvars$family
  if (is_multivariate(family)) {
    if (has_univariate(family)) {
      cgvars$default <- vapply(
        cvars,
        function(x) {
          paste_rows(
            lines_wrap("data", "default", idt, backend, x),
            do.call("prior_data_lines", c(idt = idt, x)),
            .indent = idt(0),
            .parse = FALSE
          )
        },
        character(1L)
      )
    }
    lines_wrap("data", family, idt, backend, cgvars)
  } else {
    cvars[[1L]]$default <- lines_wrap(
      "data", "default", idt, backend, cvars[[1L]]
    )
    lines_wrap("data", family, idt, backend, cvars[[1L]])
  }
}

#' @describeIn create_function Create the 'Transformed Data'
#'   Block of the Stan Model Code
#' @noRd
create_transformed_data <- function(idt, backend, cg, cvars, cgvars, mvars) {
  n_cg <- n_unique(cg)
  declarations <- character(n_cg)
  statements <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    tr_data <- create_transformed_data_lines(
      idt, backend, cvars[cg_idx], cgvars[[i]]
    )
    declarations[i] <- tr_data$declarations
    statements[i] <- tr_data$statements
  }
  has_lfactor <- any(vapply(cvars, "[[", logical(1L), "has_lfactor"))
  paste_rows(
    "transformed data {",
    declarations,
    onlyif(has_lfactor, "vector[2 * N] QR_Q = create_Q(N);"),
    ifelse_(
      stan_supports_array_keyword(backend),
      "array[T] int seq1T = linspaced_int_array(T, 1, T);",
      "int seq1T[T] = linspaced_int_array(T, 1, T);"
    ),
    statements,
    "}",
    .indent = idt(c(0, 0, 1, 1, 0, 0)),
    .parse = FALSE
  )
}

#' Create Transformed Data Lines for a Distribution
#'
#' @noRd
create_transformed_data_lines <- function(idt, backend, cvars, cgvars) {
  family <- cgvars$family
  if (is_multivariate(family)) {
    cgvars$default <- lapply(
      cvars,
      function(x) {
        lines_wrap("transformed_data", "default", idt, backend, x)
      }
    )
    lines_wrap("transformed_data", family, idt, backend, cgvars)
  } else {
    cvars[[1L]]$default <- lines_wrap(
      "transformed_data", "default", idt, backend, cvars[[1L]]
    )
    lines_wrap("transformed_data", family, idt, backend, cvars[[1L]])
  }
}

#' @describeIn create_function Create the 'Parameters'
#'   Block of the Stan Model Code
#' @noRd
create_parameters <- function(idt, backend, cg, cvars, cgvars, mvars) {
  spline_def <- mvars$spline_def
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
    mvars$random_def$M > 0L,
    paste_rows(
      "// Random group-level effects",
      onlyif(
        mvars$random_def$correlated,
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
    identical(length(mvars$lfactor_def$responses), 0L),
    "",
    paste_rows(
      "// Latent factor splines",
      onlyif(
        mvars$lfactor_def$correlated,
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
    parameters_text[i] <- create_parameters_lines(
      idt, backend, cvars[cg_idx], cgvars[[i]]
    )
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

#' Create Parameters Lines for a Distribution
#'
#' @noRd
create_parameters_lines <- function(idt, backend, cvars, cgvars) {
  family <- cgvars$family
  if (is_multivariate(family)) {
    if (is_multinomial(family)) {
      cgvars$univariate <- ulapply(
        cgvars$y[-1L],
        function(s) {
          cvars[[1L]]$ydim <- cgvars$y_cg
          cvars[[1L]]$y <- s
          lines_wrap(
            "parameters", "default", idt, backend, cvars[[1L]]
          )
        }
      )
    } else {
      univariate <- ""
      for (j in seq_along(cvars)) {
        cvars[[j]]$default <- lines_wrap(
          "parameters", "default", idt, backend, cvars[[j]]
        )
        univariate <- ifelse_(
          has_univariate(family),
          paste_rows(
            univariate,
            lines_wrap(
              "parameters",
              get_univariate(family),
              idt,
              backend,
              cvars[[j]]
            )
          ),
          paste_rows(univariate, cvars[[j]]$default)
        )
      }
      cgvars$univariate <- univariate
    }
    lines_wrap("parameters", family, idt, backend, cgvars)
  } else {
    if (is_categorical(family)) {
      cvars[[1L]]$default <- lapply(
        cvars[[1L]]$categories[-1L],
        function(s) {
          cvars[[1L]]$ydim <- cvars[[1L]]$y
          cvars[[1L]]$y <- paste0(cvars[[1L]]$y, "_", s)
          lines_wrap(
            "parameters", "default", idt, backend, cvars[[1L]]
          )
        }
      )
    } else {
      cvars[[1L]]$default <- lines_wrap(
        "parameters", "default", idt, backend, cvars[[1L]]
      )
    }
    lines_wrap("parameters", family, idt, backend, cvars[[1L]])
  }
}

#' @describeIn create_function Create the 'Transformed Parameters'
#'   Block of the Stan Model Code
#' @noRd
create_transformed_parameters <- function(idt, backend,
                                          cg, cvars, cgvars, mvars) {
  random_text <- ""
  M <- mvars$random_def$M
  if (M > 0L) {
    Ks <- mvars$Ks[mvars$Ks > 0L]
    y <- names(Ks)
    cKs1 <- cumsum(c(1L, Ks[-length(Ks)]))
    cKs2 <- cumsum(Ks)
    if (mvars$random_def$noncentered) {
      random_text <- ifelse_(
        mvars$random_def$correlated,
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
      random_text <- paste_rows(
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
  psis <- mvars$lfactor_def$responses
  P <- length(psis)
  if (P > 0L) {
    if (mvars$lfactor_def$noncentered_psi) {
      tau_psi <- ifelse(
        mvars$lfactor_def$nonzero_lambda,
        paste0("tau_psi_", psis, " * "),
        ""
      )
      # create noise terms for latent factors
      lfactor_text <- ifelse_(
        mvars$lfactor_def$correlated,
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
      lfactor_text <- paste_rows(
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
    tr_pars <- create_transformed_parameters_lines(
      idt, backend, cvars[cg_idx], cgvars[[i]]
    )
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

#' Create Transformed Parameters Lines for a Distribution
#'
#' @noRd
create_transformed_parameters_lines <- function(idt, backend, cvars, cgvars) {
  family <- cgvars$family
  if (is_multivariate(family)) {
    if (is_multinomial(family)) {
      cgvars$default <- lapply(
        cgvars$y[-1L],
        function(s) {
          cvars[[1L]]$ydim <- cgvars$y_cg
          cvars[[1L]]$y <- s
          lines_wrap(
            "transformed_parameters", "default", idt, backend, cvars[[1L]]
          )
        }
      )
    } else {
      cgvars$default <- lapply(
        cvars,
        function(x) {
          lines_wrap("transformed_parameters", "default", idt, backend, x)
        }
      )
    }
    lines_wrap("transformed_parameters", family, idt, backend, cgvars)
  } else {
    if (is_categorical(family)) {
      cvars[[1L]]$default <- lapply(
        cvars[[1L]]$categories[-1L],
        function(s) {
          cvars[[1L]]$ydim <- cvars[[1L]]$y
          cvars[[1L]]$y <- paste0(cvars[[1L]]$y, "_", s)
          lines_wrap(
            "transformed_parameters", "default", idt, backend, cvars[[1L]]
          )
        }
      )
    } else {
      cvars[[1L]]$default <- lines_wrap(
        "transformed_parameters", "default", idt, backend, cvars[[1L]]
      )
    }
    lines_wrap("transformed_parameters", family, idt, backend, cvars[[1L]])
  }
}

#' @describeIn create_function Create the 'Model' Block of the Stan Model Code
#' @noRd
create_model <- function(idt, backend, cg, cvars, cgvars, mvars, threading) {
  spline_def <- mvars$spline_def
  spline_text <- ""
  if (!is.null(spline_def) && spline_def$shrinkage) {
    xi_prior <- mvars$common_priors
    xi_prior <- xi_prior[xi_prior$parameter == "xi", "prior"]
    spline_text <- paste_rows("xi[1] ~ {xi_prior};", .indent = idt(1))
  }
  random_text <- ""
  if (mvars$random_def$M > 0L) {
    if (mvars$random_def$correlated) {
      L_prior <- mvars$common_priors
      L_prior <- L_prior[L_prior$parameter == "L_nu", "prior"]
      random_text <- ifelse_(
        mvars$random_def$noncentered,
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
      M <- mvars$random_def$M
      random_text <- ifelse_(
        mvars$random_def$noncentered,
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
  psis <- mvars$lfactor_def$responses
  P <- length(psis)
  if (P > 0L) {
    omega1 <- paste0("omega_raw_psi_1_", psis, collapse = ", ")
    if (mvars$lfactor_def$correlated) {
      L_prior <- mvars$common_priors
      L_prior <- L_prior[L_prior$parameter == "L_lf", "prior"]
      if (mvars$lfactor_def$noncentered_psi) {
        lfactor_text <- paste_rows(
          "to_vector(omega_raw_psi) ~ std_normal();",
          "L_lf ~ {L_prior};",
          .indent = idt(c(1, 1))
        )
      } else {
        if (any(mvars$lfactor_def$nonzero_lambda)) {
          tau <- paste0(
            ifelse(
              mvars$lfactor_def$nonzero_lambda,
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
            "omega_raw_psi[, 1] ~ multi_normal_cholesky(omega1, Ltau);",
            "for (i in 2:(D - 1)) {{",
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
            "omega_raw_psi[, 1] ~ multi_normal_cholesky(omega1, L_lf);",
            "for (i in 2:(D - 1)) {{",
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
      if (mvars$lfactor_def$noncentered_psi) {
        lfactor_text <- paste_rows(
          "to_vector(omega_raw_psi) ~ std_normal();",
          .indent = idt(1)
        )
      } else {
        if (any(mvars$lfactor_def$nonzero_lambda)) {
          tau <- paste0(
            ifelse(
              mvars$lfactor_def$nonzero_lambda,
              paste0("tau_psi_", psis),
              "1"
            ),
            collapse = ", "
          )
          lfactor_text <- paste_rows(
            "{{",
            "vector[P] tau_psi = [{tau}]';",
            "vector[P] omega1 = [{omega1}]';",
            "omega_raw_psi[, 1] ~ normal(omega1, tau_psi);",
            "for (i in 2:(D - 1)) {{",
            "omega_raw_psi[, i] ~ normal(omega_raw_psi[, i - 1], tau_psi);",
            "}}",
            "}}",
            .indent = idt(c(1, 2, 2, 2, 2, 3, 2, 1))
          )
        } else {
          lfactor_text <- paste_rows(
            "vector[P] omega1 = [{omega1}]';",
            "omega_raw_psi[, 1] ~ normal(omega1, 1);",
            "for (i in 2:(D - 1)) {{",
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
    model_text[i] <- create_model_lines(
      idt, backend, cvars[cg_idx], cgvars[[i]], threading = threading
    )
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

#' Create Model Lines for a Distribution
#'
#' @noRd
create_model_lines <- function(idt, backend, cvars, cgvars, mvars, threading) {
  family <- cgvars$family
  if (is_multivariate(family)) {
    cgvars$backend <- backend
    lines_wrap(
      "model",
      family,
      idt,
      backend,
      list(cvars = cvars, cgvars = cgvars, threading = threading)
    )
  } else {
    cvars[[1L]]$backend <- backend
    if (is_categorical(family)) {
      cvars[[1L]]$priors <- lapply(
        cvars[[1L]]$categories[-1L],
        function(s) {
          cvars[[1L]]$y <- paste0(cvars[[1L]]$y, "_", s)
          cvars[[1L]]$prior_distr <- cvars[[1L]]$prior_distr[[s]]
          do.call(prior_lines, c(cvars[[1L]], idt = idt))
        }
      )
      cvars[[1L]]$backend <- backend
    } else {
      cvars[[1L]]$priors <- do.call(prior_lines, c(cvars[[1L]], idt = idt))
    }
    cvars[[1L]]$threading <- threading
    lines_wrap("model", family, idt, backend, cvars[[1L]])
  }
}

#' @describeIn create_function Create the 'Generated Quantities'
#'   Block of the Stan Model Code
#' @noRd
create_generated_quantities <- function(idt, backend,
                                        cg, cvars, cgvars, mvars) {
  gen_nu <- ""
  M <- mvars$random_def$M
  if (M > 1L && mvars$random_def$correlated) {
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
  psis <- mvars$lfactor_def$responses
  P <- mvars$lfactor_def$P
  if (P > 0L && mvars$lfactor_def$correlated) {
    # evaluate number of corrs to avoid Stan warning about integer division
    tau <- paste0(
      ifelse(
        cvars$lfactor_def$nonzero_lambda,
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
    generated_quantities_text[i] <- create_generated_quantities_lines(
      idt, backend, cvars[cg_idx], cgvars[[i]]
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

create_generated_quantities_lines <- function(idt, backend, cvars, cgvars) {
  family <- cgvars$family
  if (is_multivariate(family)) {
    lines_wrap("generated_quantities", family, idt, backend, cgvars)
  } else {
    lines_wrap("generated_quantities", family, idt, backend, cvars[[1L]])
  }
}

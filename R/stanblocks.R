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
#' @param cvars \[`list()`]\cr The `channel_vars` component of
#' @param cgvars \[`list()`]\cr The `channel_group_vars` component of
#'   [prepare_stan_input()] output.
#' @param ... Not used.
#' @noRd
create_blocks.default <- function(dformula, indent = 2L, cvars, cgvars,
                                  backend, ...) {
  idt <- indenter_(indent)
  paste_rows(
    create_functions(dformula, idt, cvars, cgvars),
    create_data(dformula, idt, cvars, cgvars),
    create_transformed_data(dformula, idt, cvars, cgvars),
    create_parameters(dformula, idt, cvars, cgvars),
    create_transformed_parameters(dformula, idt, cvars, cgvars),
    create_model(dformula, idt, cvars, cgvars, backend),
    create_generated_quantities(dformula, idt, cvars, cgvars),
    .parse = FALSE
  )
}

#' Create the 'Functions' Block of the Stan Model Code
#'
#' @inheritParams create_blocks.default
#' @param idt \[`function`] An indentation function created by [indenter_()]
#' @noRd
create_functions <- function(dformula, idt, cvars, cgvars) {
  NULL
}

#' @describeIn create_function Create The 'Data' Block of the Stan Model Code
#' @noRd
create_data <- function(dformula, idt, cvars, cgvars) {
  has_splines <- any(vapply(cvars, "[[", logical(1L), "has_varying")) ||
    any(vapply(cvars, "[[", logical(1L), "has_varying_intercept")) ||
    any(vapply(cvars, "[[", logical(1L), "has_lfactor"))
  K <- sum(vapply(cvars, "[[", integer(1), "K"))
  M <- attr(cvars, "random_def")$M
  P <- attr(cvars, "lfactor_def")$P
  mtext <- paste_rows(
    "int<lower=1> T; // number of time points",
    "int<lower=1> N; // number of individuals",
    onlyif(K > 0,
      "int<lower=0> K; // total number of covariates across all channels"),
    onlyif(K > 0,
      "matrix[N, K] X[T]; // covariates as an array of N x K matrices"),
    onlyif(K > 0,
      "row_vector[K] X_m; // Means of all covariates at first time point"),
    onlyif(has_splines, "int<lower=1> D; // number of B-splines"),
    onlyif(has_splines, "matrix[D, T] Bs; // B-spline basis matrix"),
    onlyif(M > 0,
      "int<lower=0> M; // number group-level effects (including intercepts)"),
    onlyif(P > 0,
      "int<lower=0> P; // number of channels with latent factor"),
    .indent = idt(1),
    .parse = FALSE
  )
  cg <- attr(dformula, "channel_groups")
  n_cg <- length(unique(cg))
  datatext <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    j <- cg_idx[1L]
    family <- dformula[[j]]$family$name
    if (is_multivariate(dformula[[j]]$family)) {
      datatext[i] <- ""
      for (k in cg_idx) {
        line_args <- c(
          list(y = dformula[[cg_idx[k]]]$response, idt = idt),
          cvars[[j]]
        )
        datatext[i] <- paste_rows(
          datatext[i],
          lines_wrap("data", "default", line_args)
        )
      }
      y <- paste(get_responses(dformula[cg_idx]), collapse = "_")
      line_args <- c(list(y = y, idt = idt), cvars[[j]])
      datatext[i] <- paste_rows(
        datatext[i],
        lines_wrap("data", family, line_args)
      )
    } else {
      line_args <- c(list(y = dformula[[j]]$response, idt = idt), cvars[[j]])
      datatext[i] <- lines_wrap("data", family, line_args)
    }
  }

  paste_rows("data {", mtext, datatext, "}", .parse = FALSE)
}

#' @describeIn create_function Create the 'Transformed Data'
#'   Block of the Stan Model Code
#' @noRd
create_transformed_data <- function(dformula, idt, cvars, cgvars) {
  n <- length(dformula)
  declarations <- character(n)
  statements <- character(n)
  for (i in seq_len(n)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = cvars[[i]]$resp, idt = idt), cvars[[i]])
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
  has_lfactor <- any(vapply(cvars, "[[", logical(1L), "has_lfactor"))
  paste_rows(
    "transformed data {",
    declarations,
    onlyif(has_lfactor, declare_AQR),
    statements,
    onlyif(has_lfactor, state_AQR),
    "}",
    .indent = idt(c(0, 0, 1, 0, 1, 0)),
    .parse = FALSE
  )
}

#' @describeIn create_function Create the 'Parameters'
#'   Block of the Stan Model Code
#' @noRd
create_parameters <- function(dformula, idt, cvars, cgvars) {
  spline_def <- attr(dformula, "splines")
  splinetext <- ifelse_(
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

  randomtext <- ifelse_(
    attr(cvars, "random_def")$M > 0L,
    paste_rows(
      "// Random group-level effects",
      onlyif(
        attr(cvars, "random_def")$correlated,
        "cholesky_factor_corr[M] L_nu; // Cholesky for correlated random effects"
      ),
      "vector<lower=0>[M] sigma_nu; // standard deviations of random effects",
      "matrix[N, M] nu_raw;",
      .indent = idt(c(1, 1, 1, 1))
    ),
    ""
  )

  lfactortext <- ifelse_(
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

  cg <- attr(dformula, "channel_groups")
  n_cg <- length(unique(cg))
  pars <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    j <- cg_idx[1L]
    family <- dformula[[j]]$family$name
    if (is_multivariate(dformula[[j]]$family)) {
      pars[i] <- ""
      for (k in cg_idx) {
        line_args <- c(
          list(y = dformula[[cg_idx[k]]]$response, idt = idt),
          cvars[[j]]
        )
        pars[i] <- paste_rows(
          pars[i],
          lines_wrap(
            "parameters",
            get_univariate(dformula[[j]]$family),
            line_args
          )
        )
      }
      y <- paste(get_responses(dformula[cg_idx]), collapse = "_")
      line_args <- c(list(y = y, idt = idt), cvars[[j]])
      pars[i] <- paste_rows(
        pars[i],
        lines_wrap("parameters", family, line_args)
      )
    } else {
      family <- dformula[[j]]$family$name
      line_args <- c(list(y = dformula[[j]]$response, idt = idt), cvars[[j]])
      pars[i] <- lines_wrap("parameters", family, line_args)
    }
  }
  paste_rows("parameters {", splinetext, randomtext, lfactortext, pars, "}",
    .parse = FALSE
  )
}

#' @describeIn create_function Create the 'Transformed Parameters'
#'   Block of the Stan Model Code
#' @noRd
create_transformed_parameters <- function(dformula, idt, cvars, cgvars) {
  randomtext <- ""
  M <- attr(cvars, "random_def")$M
  if (M > 0L) {
    Ks <- vapply(cvars, "[[", integer(1L), "K_random")
    Ks <- Ks[Ks > 0]
    y <- names(Ks)
    cKs1 <- cumsum(c(1, Ks[-length(Ks)]))
    cKs2 <- cumsum(Ks)
    if (attr(cvars, "random_def")$noncentered) {
      randomtext <- ifelse_(
        attr(cvars, "random_def")$correlated,
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
      lfactortext <- ifelse_(
        attr(cvars, "lfactor_def")$correlated,
        paste_rows(
          "matrix[P, D - 1] omega_psi = L_lf * omega_raw_psi;",
          glue::glue(
            "row_vector[D] omega_psi_{psis} = \\
            append_col(0, {tau_psi}omega_psi[{1:P}, ]);"
          ),
          .indent = idt(1)
        ),
        paste_rows(
          glue::glue(
            "row_vector[D] omega_psi_{psis} = \\
            append_col(0, {tau_psi}omega_raw_psi[{1:P}, ]);"
          ),
          .indent = idt(1)
        )
      )
    } else {
      lfactortext <-
        paste_rows(
          glue::glue(
            "row_vector[D] omega_psi_{psis} = \\
             append_col(0, omega_raw_psi[{1:P}, ]);"
          ),
          .indent = idt(1)
        )
    }
  }

  n <- length(dformula)
  declarations <- character(n)
  statements <- character(n)
  for (i in seq_len(n)) {
    family <- dformula[[i]]$family$name
    line_args <- c(list(y = cvars[[i]]$resp, idt = idt), cvars[[i]])
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
create_model <- function(dformula, idt, cvars, cgvars, backend) {
  spline_def <- attr(dformula, "splines")
  splinetext <- ""
  if (!is.null(spline_def) && spline_def$shrinkage) {
    xi_prior <- attr(cvars, "common_priors")
    xi_prior <- xi_prior[xi_prior$parameter == "xi", "prior"]
    splinetext <- paste_rows("xi[1] ~ {xi_prior};", .indent = idt(1))
  }
  randomtext <- ""
  if (attr(cvars, "random_def")$M > 0L) {
    if (attr(cvars, "random_def")$correlated) {
      L_prior <- attr(cvars, "common_priors")
      L_prior <- L_prior[L_prior$parameter == "L_nu", "prior"]
      randomtext <- ifelse_(
        attr(cvars, "random_def")$noncentered,
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
      randomtext <- ifelse_(
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

  lfactortext <- ""
  psis <- attr(cvars, "lfactor_def")$responses
  P <- length(psis)
  if (P > 0L) {
    omega1 <- paste0("omega_raw_psi_1_", psis,
      collapse = ", "
    )
    if (attr(cvars, "lfactor_def")$correlated) {
      L_prior <- attr(cvars, "common_priors")
      L_prior <- L_prior[L_prior$parameter == "L_lf", "prior"]
      if (attr(cvars, "lfactor_def")$noncentered_psi) {
        lfactortext <- paste_rows(
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
          lfactortext <- paste_rows(
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
          lfactortext <- paste_rows(
            "L_lf ~ {L_prior};",
            "vector[P] omega1 = [{omega1}]';",
            "omega_raw_psi[, 2] ~ multi_normal_cholesky(omega1, Ltau);",
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
        lfactortext <- paste_rows(
          "to_vector(omega_raw_psi) ~ std_normal();",
          .indent = idt(1)
        )
      } else {
        if (any(attr(cvars, "lfactor_def")$nonzero_lambda)) {
          tau <- paste0(ifelse(attr(cvars, "lfactor_def")$nonzero_lambda,
            paste0("tau_psi_", psis), "1"
          ), collapse = ", ")
          lfactortext <- paste_rows(
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
          lfactortext <- paste_rows(
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
  cg <- attr(dformula, "channel_groups")
  n_cg <- length(unique(cg))
  mod <- character(n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    j <- cg_idx[1L]
    family <- dformula[[j]]$family$name
    if (is_multivariate(dformula[[j]]$family)) {
      mod[i] <- ""
      # for (k in cg_idx) {
      #   line_args <- c(
      #     list(y = dformula[[cg_idx[k]]]$response, idt = idt, backend = backend),
      #     cvars[[j]]
      #   )
      #   mod[i] <- paste_rows(
      #     mod[i],
      #     lines_wrap("model", "default", line_args)
      #   )
      # }
     # y <- paste(get_responses(dformula[cg_idx]), collapse = "_")
      line_args <- c(
        list(y = get_responses(dformula[cg_idx]), idt = idt),
        cgvars[[i]]
      )
      mod[i] <- paste_rows(
        mod[i],
        lines_wrap("model", family, line_args), .parse = FALSE
      )
    } else {
      family <- dformula[[j]]$family$name
      line_args <- c(
        list(y = cvars[[j]]$resp, idt = idt, backend = backend),
        cvars[[j]]
      )
      mod[i] <- lines_wrap("model", family, line_args)
    }
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
create_generated_quantities <- function(dformula, idt, cvars, cgvars) {
  gen_nu <- ""
  M <- attr(cvars, "random_def")$M
  if (M > 1 && attr(cvars, "random_def")$correlated) {
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
      paste0("corr_matrix[P] corr_matrix_psi = ",
        "multiply_lower_tri_self_transpose(L_lf);"),
      "vector<lower=-1,upper=1>[{(P * (P - 1L)) %/% 2L}] corr_psi;",
      "for (k in 1:P) {{",
      "for (j in 1:(k - 1)) {{",
      "corr_psi[choose(k - 1, 2) + j] = corr_matrix_psi[j, k];",
      "}}",
      "}}",
      .indent = idt(c(1, 1, 1, 2, 3, 2, 1))
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

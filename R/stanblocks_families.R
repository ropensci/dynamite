#' Wrapper to parse Stan model blocks based on the family of the response
#'
#' @param prefix A character string indicating the Stan model block name
#' @param family A supported family name
#' @param args Channel specific component of `model_vars`
#'   (see [create_blocks()])
#'
#' @noRd
lines_wrap <- function(prefix, family, args) {
  lines_expr <- paste0(prefix, "_lines_", family)
  lines_env <- list2env(args)
  eval(eval(as.name(lines_expr)), envir = lines_env)
}

#' Checks if a prior definition is vectorizable
#'
#' @param x A character string
#'
#' @noRd
vectorizable_prior <- function(x) {
  length(x) == 1 && !grepl("\\(", x)
}

# Data block --------------------------------------------------------------

data_lines_default <- quote({
  paste_rows(
    "// Data for response {y}",
    onlyif(has_missing, "int<lower=0> obs_{y}[N, T];"),
    onlyif(has_missing, "int<lower=0> n_obs_{y}[T];"),
    .indent = idt(1)
  )
})

data_lines_categorical <- quote({
  dtext_def <- eval(data_lines_default)
  dtext <- paste_rows(
    "int<lower=1> {y}[T, N];",
    .indent = idt(1)
  )
  paste_rows(dtext_def, dtext, .parse = FALSE)
})

data_lines_gaussian <- quote({
  dtext_def <- eval(data_lines_default)
  dtext <- paste_rows(
    "matrix[N, T] {y};",
    .indent = idt(1)
  )
  paste_rows(dtext_def, dtext, .parse = FALSE)
})

data_lines_binomial <- quote({
  dtext_def <- eval(data_lines_default)
  dtext <- paste_rows(
    "int<lower=0> {y}[T, N];",
    "int<lower=1> trials_{y}[T, N];",
    .indent = idt(1)
  )
  paste_rows(dtext_def, dtext, .parse = FALSE)
})

data_lines_bernoulli <- quote({
  dtext_def <- eval(data_lines_default)
  dtext <- paste_rows(
    "int<lower=0,upper=1> {y}[T, N];",
    .indent = idt(1)
  )
  paste_rows(dtext_def, dtext, .parse = FALSE)
})

data_lines_poisson <- quote({
  dtext_def <- eval(data_lines_default)
  dtext <- paste_rows(
    "int<lower=0> {y}[T, N];",
    onlyif(has_offset, "real offset_{y}[T, N];"),
    .indent = idt(1)
  )
  paste_rows(dtext_def, dtext, .parse = FALSE)
})

data_lines_negbin <- quote({
  dtext_def <- eval(data_lines_default)
  dtext <- paste_rows(
    "int<lower=0> {y}[T, N];",
    .indent = idt(1)
  )
  paste_rows(dtext_def, dtext, .parse = FALSE)
})

# Transformed data block --------------------------------------------------

transformed_data_lines_default <- quote({
  mtext_fixed <- ""
  mtext_varying <- ""
  mtext_tau <- ""
  if (write_beta) {
    i <- rep(1:K_fixed, beta_prior_npars)
    j <- rep(1:beta_prior_npars, each = K_fixed)
    mtext_fixed <- paste_rows(
      "matrix[{K_fixed}, {beta_prior_npars}] beta_prior_pars_{y};",
      "beta_prior_pars_{y}[{i},{j}] = {beta_prior_pars};",
      .indent = idt(1)
    )
  }
  if (write_delta) {
    i <- rep(1:K_varying, delta_prior_npars)
    j <- rep(1:delta_prior_npars, each = K_varying)
    mtext_varying <- paste_rows(
      "matrix[{K_varying}, {delta_prior_npars}] delta_prior_pars_{y};",
      "delta_prior_pars_{y}[{i},{j}] = {delta_prior_pars};",
      .indent = idt(1)
    )
  }
  if (write_tau) {
    i <- rep(1:K_varying, tau_prior_npars)
    j <- rep(1:tau_prior_npars, each = K_varying)
    mtext_tau <- paste_rows(
      "matrix[{K_varying}, {tau_prior_npars}] tau_prior_pars_{y};",
      "tau_prior_pars_{y}[{i},{j}] = {tau_prior_pars};",
      .indent = idt(1)
    )
  }
  paste_rows(mtext_fixed, mtext_varying, mtext_tau, .parse = FALSE)
})

transformed_data_lines_categorical <- quote({
  mtext <- paste_rows(
    "vector[{K}] zeros_K_{y} = rep_vector(0, {K});",
    "vector[{S}] zeros_S_{y} = rep_vector(0, {S});",
    .indent = idt(1)
  )
  mtext_alpha <- ""
  mtext_fixed <- ""
  mtext_varying <- ""
  mtext_tau <- ""
  if (write_alpha) {
    k <- S - 1
    i <- rep(1:k, alpha_prior_npars)
    j <- rep(1:alpha_prior_npars, each = k)
    mtext_alpha <- paste_rows(
      "matrix[{k}, {alpha_prior_npars}] alpha_prior_pars_{y};",
      "alpha_prior_pars_{y}[{i},{j}] = {alpha_prior_pars};",
      .indent = idt(1)
    )
  }
  if (write_beta) {
    k <- (K_fixed * (S - 1))
    i <- rep(1:k, beta_prior_npars)
    j <- rep(1:beta_prior_npars, each = k)
    mtext_fixed <- paste_rows(
      "matrix[{k}, {beta_prior_npars}] beta_prior_pars_{y};",
      "beta_prior_pars_{y}[{i},{j}] = {beta_prior_pars};",
      .indent = idt(1)
    )
  }
  if (write_delta) {
    k <- (K_varying * (S - 1))
    i <- rep(1:k, delta_prior_npars)
    j <- rep(1:delta_prior_npars, each = k)
    mtext_varying <- paste_rows(
      "matrix[{k}, {delta_prior_npars}] delta_prior_pars_{y};",
      "delta_prior_pars_{y}[{i},{j}] = {delta_prior_pars};",
      .indent = idt(1)
    )
  }
  if (write_tau) {
    i <- rep(1:K_varying, tau_prior_npars)
    j <- rep(1:tau_prior_npars, each = K_varying)
    mtext_tau <- paste_rows(
      "matrix[{K_varying}, {tau_prior_npars}] tau_prior_pars_{y};",
      "tau_prior_pars_{y}[{i},{j}] = {tau_prior_pars};",
      .indent = idt(1)
    )
  }
  paste_rows(mtext, mtext_alpha, mtext_fixed, mtext_varying, mtext_tau,
             .parse = FALSE)
})

transformed_data_lines_gaussian <- quote({
  eval(transformed_data_lines_default)
})

transformed_data_lines_binomial <- quote({
  eval(transformed_data_lines_default)
})

transformed_data_lines_bernoulli <- quote({
  eval(transformed_data_lines_default)
})

transformed_data_lines_poisson <- quote({
  eval(transformed_data_lines_default)
})

transformed_data_lines_negbin <- quote({
  eval(transformed_data_lines_default)
})

# Parameters block --------------------------------------------------------

parameters_lines_default <- quote({

  intercept <- ""
  if (has_fixed_intercept || has_varying_intercept) {
    intercept <- "real a_{y};"
  }
  oname <- ifelse_(noncentered, "omega_raw_", "omega_")
  paste_rows(
    onlyif(has_fixed, "vector[{K_fixed}] beta_{y};"),
    onlyif(has_varying, "matrix[{K_varying}, D] {oname}{y};"),
    onlyif(has_varying, "vector<lower={lb}>[{K_varying}] tau_{y};"),
    intercept,
    onlyif(has_varying_intercept, "row_vector[D - 1] omega_raw_alpha_{y};"),
    onlyif(has_varying_intercept, "real<lower={lb}> tau_alpha_{y};"),
    .indent = idt(1)
  )
})

parameters_lines_categorical <- quote({

  intercept <- ""
  if (has_fixed_intercept || has_varying_intercept) {
    intercept <- "vector[{S - 1}] a_{y};"
  }
  oname <- ifelse_(noncentered, "omega_raw_", "omega_")
  paste_rows(
    onlyif(has_fixed, "matrix[{K_fixed}, {S - 1}] beta_{y};"),
    onlyif(has_varying, "matrix[{K_varying}, D] {oname}{y}[{S - 1}];"),
    onlyif(has_varying, "vector<lower={lb}>[{K_varying}] tau_{y};"),
    intercept,
    onlyif(has_varying_intercept, "row_vector[D - 1] omega_raw_alpha_{y}[{S - 1}];"),
    onlyif(has_varying_intercept, "real<lower={lb}> tau_alpha_{y};"),
    .indent = idt(1)
  )
})

parameters_lines_gaussian <- quote({
  paste_rows(
    eval(parameters_lines_default),
    "real<lower=0> sigma_{y};",
    .indent = idt(c(0, 1))
  )
})

parameters_lines_binomial <- quote({
  eval(parameters_lines_default)
})

parameters_lines_bernoulli <- quote({
  eval(parameters_lines_default)
})

parameters_lines_poisson <- quote({
  eval(parameters_lines_default)
})

parameters_lines_negbin <- quote({
  paste_rows(
    eval(parameters_lines_default),
    "real<lower=0> phi_{y};",
    .indent = idt(c(0, 1))
  )
})


# Transformed parameters --------------------------------------------------

transformed_parameters_lines_default <- quote({
  mtext_varying_noncentered <- ""
  mtext_varying <- ""
  mtext_intercept <- ""
  if (has_varying) {
    if (noncentered) {
      lambda_term <- ifelse_(shrinkage, " * lambda[i - 1];", ";")
      mtext_varying_noncentered <- paste_rows(
        "matrix[{K_varying}, D] omega_{y};",
        "omega_{y}[, 1] = omega_raw_{y}[, 1];",
        "for (i in 2:D) {{",
          "omega_{y}[, i] = omega_{y}[, i - 1] + omega_raw_{y}[, i] .* tau_{y}{lambda_term}",
        "}}",
        .indent = idt(c(1, 1, 1, 2, 1)),
        .parse = FALSE
      )
    }
    mtext_varying <- paste_rows(
      "vector[{K_varying}] delta_{y}[T];",
      "for (t in 1:T) {{",
        "delta_{y}[t] = omega_{y} * Bs[, t];",
      "}}",
      .indent = idt(c(1, 1, 2, 1)),
      .parse = FALSE
    )
  }
  if (has_varying_intercept) {
    mtext_omega_alpha_1 <- paste_rows(
      "real alpha_{y}[T];",
      "real omega_alpha_1_{y};",
      "{{",
        "vector[{K}] gamma_{y};",
        onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y};"),
        onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[1];"),
        "omega_alpha_1_{y} = a_{y} - X_m[{{{cs(J)}}}] * gamma_{y};",
      "}}",
      .indent = idt(c(1, 1, 1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
    if (noncentered) {
      lambda_term <- ifelse_(shrinkage, " * lambda[i - 1];", ";")
      mtext_omega_alpha <- paste_rows(
        "row_vector[D] omega_alpha_{y};",
        "omega_alpha_{y}[1] = omega_alpha_1_{y};",
        "for (i in 2:D) {{",
          "omega_alpha_{y}[i] = omega_alpha_{y}[i - 1] + omega_raw_alpha_{y}[i - 1] * tau_alpha_{y}{lambda_term}",
        "}}",
        .indent = idt(c(1, 1, 1, 2, 1)),
        .parse = FALSE
      )
    } else {
      mtext_omega_alpha <- paste_rows(
        "row_vector[D] omega_alpha_{y};",
        "omega_alpha_{y}[1] = omega_alpha_1_{y};",
        "omega_alpha_{y}[2:D] = omega_raw_alpha_{y};",
        .indent = idt(c(1, 1, 1)),
        .parse = FALSE
      )
    }
    mtext_intercept <- paste_rows(
      mtext_omega_alpha_1,
      mtext_omega_alpha,
      "for (t in 1:T) {{",
        "alpha_{y}[t] = omega_alpha_{y} * Bs[, t];",
      "}}",
      .indent = idt(c(0, 0, 1, 2, 1)),
      .parse = FALSE
    )
  }
  if (has_fixed_intercept) {
    mtext_intercept <- paste_rows(
      "real alpha_{y};",
      "{{",
        "vector[{K}] gamma_{y};",
        onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y};"),
        onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[1];"),
        "alpha_{y} = a_{y} - X_m[{{{cs(J)}}}] * gamma_{y};",
      "}}",
      .indent = idt(c(1, 1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
  }
  paste_rows(mtext_varying_noncentered, mtext_varying, mtext_intercept,
    .indent = idt(c(1, 0, 0))
  )
})

transformed_parameters_lines_categorical <- quote({
  mtext_varying_noncentered <- ""
  mtext_varying <- ""
  mtext_intercept <- ""
  if (has_varying) {
    if (noncentered) {
      lambda_term <- ifelse_(shrinkage, " * lambda[i - 1];", ";")
      mtext_varying_noncentered <- paste_rows(
        "matrix[{K_varying}, D] omega_{y}[{S - 1}];",
        "for (s in 1:{S - 1}) {{",
          "omega_{y}[s, , 1] = omega_raw_{y}[s, , 1];",
          "for (i in 2:D) {{",
            "omega_{y}[s, , i] = omega_{y}[s, , i - 1] + omega_raw_{y}[s, , i] .* tau_{y}{lambda_term}",
          "}}",
        "}}",
        .indent = idt(c(1, 1, 2, 2, 3, 2, 1)),
        .parse = FALSE
      )
    }
    mtext_varying <- paste_rows(
      "matrix[{K_varying}, {S - 1}] delta_{y}[T];",
      "for (s in 1:{S - 1}) {{",
        "for (t in 1:T) {{",
          "delta_{y}[t, , s] = omega_{y}[s] * Bs[, t];",
        "}}",
      "}}",
      .indent = idt(c(1, 1, 2, 3, 2, 1)),
      .parse = FALSE
    )
  }
  if (has_varying_intercept) {
    mtext_omega_alpha_1 <- paste_rows(
      "vector[{S - 1}] alpha_{y}[T];",
      "real omega_alpha_1_{y}[{S - 1}];",
      "for (s in 1:{S - 1}) {{",
        "vector[{K}] gamma_{y};",
        onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y}[, s];"),
        onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[1, , s];"),
        "omega_alpha_1_{y}[s] = a_{y}[s] - X_m[{{{cs(J)}}}] * gamma_{y};",
      "}}",
      .indent = idt(c(1, 1, 1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
    if (noncentered) {
      lambda_term <- ifelse_(shrinkage, " * lambda[i - 1];", ";")
      mtext_omega_alpha <- paste_rows(
        "row_vector[D] omega_alpha_{y}[{S - 1}];",
        "for (s in 1:{S - 1}) {{",
          "omega_alpha_{y}[s, 1] = omega_alpha_1_{y}[s];",
          "for (i in 2:D) {{",
            "omega_alpha_{y}[s, i] = omega_alpha_{y}[s, i - 1] + omega_raw_alpha_{y}[s, i - 1] * tau_alpha_{y}{lambda_term}",
          "}}",
        "}}",
        .indent = idt(c(1, 1, 2, 2, 3, 2, 1)),
        .parse = FALSE
      )
    } else {
      mtext_omega_alpha <- paste_rows(
        "row_vector[D] omega_alpha_{y}[{S - 1}];",
        "for (s in {S - 1}) {{",
          "omega_alpha_{y}[s, 1] = omega_alpha_1_{y}[s];",
          "omega_alpha_{y}[s, 2:D] = omega_raw_alpha_{y}[s];",
        "}}",
        .indent = idt(c(1, 1, 2, 2, 1)),
        .parse = FALSE
      )
    }
    mtext_intercept <- paste_rows(
      mtext_omega_alpha_1,
      mtext_omega_alpha,
      "for (t in 1:T) {{",
        "alpha_{y}[t] = omega_alpha_{y} * Bs[, t];",
      "}}",
      .indent = idt(c(0, 0, 1, 2, 1)),
      .parse = FALSE
    )
  }
  if (has_fixed_intercept) {
    mtext_intercept <- paste_rows(
      "vector[{S - 1}] alpha_{y};",
      "for (s in 1:{S - 1}) {{",
        "vector[{K}] gamma_{y};",
        onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y}[, s];"),
        onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[1, , s];"),
        "alpha_{y}[s] = a_{y}[s] - X_m[{{{cs(J)}}}] * gamma_{y};",
      "}}",
      .indent = idt(c(1, 1, 2, 2, 2, 2, 1)),
      .parse = FALSE
    )
  }
  paste_rows(mtext_varying_noncentered, mtext_varying, mtext_intercept,
    .indent = idt(c(1, 0, 0))
  )
})

transformed_parameters_lines_binomial <- quote({
  eval(transformed_parameters_lines_default)
})

transformed_parameters_lines_gaussian <- quote({
  eval(transformed_parameters_lines_default)
})

transformed_parameters_lines_bernoulli <- quote({
  eval(transformed_parameters_lines_default)
})

transformed_parameters_lines_poisson <- quote({
  eval(transformed_parameters_lines_default)
})

transformed_parameters_lines_negbin <- quote({
  eval(transformed_parameters_lines_default)
})

# Model block -------------------------------------------------------------

model_lines_default <- quote({
  mtext_intercept <- ""
  mtext_fixed <- ""
  mtext_varying <- ""
  mtext_tau <- ""

  if (has_fixed_intercept || has_varying_intercept) {
    mtext_intercept <- "a_{y} ~ {alpha_prior_distr};"
    if (has_varying_intercept) {
      if (noncentered) {
        mtext_omega <- "omega_raw_alpha_{y} ~ std_normal();"
      } else {
        lambda_term <- ifelse_(shrinkage, " * lambda[i - 1]", "")
        mtext_omega <- paste_rows(
          "omega_raw_alpha_{y}[2] ~ normal(omega_alpha_1_{y}, tau_alpha_{y}{lambda_term});",
          "for (i in 3:(D - 1)) {{",
            "omega_raw_alpha_{y}[i] ~ normal(omega_raw_alpha_{y}[i - 1], tau_alpha_{y}{lambda_term});",
          "}}",
          .indent = idt(c(1, 1, 2, 1)),
          .parse = FALSE
        )
      }
      mtext_intercept <- paste_rows(
        mtext_intercept,
        mtext_omega,
        "tau_alpha_{y} ~ {tau_alpha_prior_distr};",
        .indent = idt(c(0, 0, 1)),
        .parse = FALSE
      )
    }
  }
  if (has_fixed) {
    if (vectorizable_prior(beta_prior_distr)) {
      np <- beta_prior_npars
      dpars_fixed <- paste0("beta_prior_pars_", y, "[, ", 1:np, "]",
                            collapse = ", ")
      mtext_fixed <- "beta_{y} ~ {beta_prior_distr}({dpars_fixed});"
    } else {
      mtext_fixed <- "beta_{y}[{{{cs(1:K_fixed)}}}] ~ {beta_prior_distr};"
    }
  }
  if (has_varying) {
    if (noncentered) {
      if (vectorizable_prior(delta_prior_distr)) {
        np <- delta_prior_npars
        dpars_varying <- paste0("delta_prior_pars_", y, "[, ", 1:np, "]",
                                collapse = ", ")
        mtext_varying <- "omega_raw_{y}[, 1] ~ {delta_prior_distr}({dpars_varying});"
      } else {
        mtext_varying <- "omega_raw_{y}[{{{cs(1:K_varying)}}}, 1] ~ {delta_prior_distr};"
      }
      mtext_varying <- paste_rows(
        mtext_varying,
        "to_vector(omega_raw_{y}[, 2:D]) ~ std_normal();",
        .indent = idt(c(0, 1)),
        .parse = FALSE
      )
    } else {
      if (vectorizable_prior(delta_prior_distr)) {
        np <- delta_prior_npars
        dpars_varying <- paste0("delta_prior_pars_", y, "[, ", 1:np, "]",
                                collapse = ", ")
        mtext_varying <- "omega_{y}[, 1] ~ {delta_prior_distr}({dpars_varying});"
      } else {
        mtext_varying <- "omega_{y}[{{{cs(1:K_varying)}}}, 1] ~ {delta_prior_distr};"
      }
      mtext_varying <- paste_rows(
        mtext_varying,
        "for (i in 2:D) {{",
        ifelse_(
          shrinkage,
          "omega_{y}[, i] ~ normal(omega_{y}[, i - 1], lambda[i - 1] * tau_{y});",
          "omega_{y}[, i] ~ normal(omega_{y}[, i - 1], tau_{y});"
        ),
        "}}",
        .indent = idt(c(0, 1, 2, 1)),
        .parse = FALSE
      )
    }
    if (vectorizable_prior(tau_prior_distr)) {
      np <- tau_prior_npars
      dpars_tau <- paste0("tau_prior_pars_", y, "[, ", 1:np, "]", collapse = ", ")
      mtext_tau <- "tau_{y} ~ {tau_prior_distr}({dpars_tau});"
    } else {
      mtext_tau <- "tau_{y}[{{{cs(1:K_varying)}}}] ~ {d};"
    }
  }
  paste_rows(mtext_intercept, mtext_fixed, mtext_varying, mtext_tau,
             .indent = idt(c(1, 1, 1, 1)))
})

model_lines_categorical <- quote({
  mtext_intercept <- ""
  mtext_fixed <- ""
  mtext_varying <- ""
  mtext_tau <- ""
  if (has_fixed_intercept || has_varying_intercept) {
    if (vectorizable_prior(alpha_prior_distr)) {
      np <- alpha_prior_npars
      dpars_alpha <- paste0("alpha_prior_pars_", y, "[, ", 1:np, "]",
                            collapse = ", ")
      mtext_intercept <- "to_vector(a_{y}) ~ {alpha_prior_distr}({dpars_alpha});"
    } else {
      s <- 1:(S - 1)
      mtext_intercept <- "a_{y}[{s}] ~ {alpha_prior_distr};"
    }
    if (has_varying_intercept) {
      if (noncentered) {
        mtext_omega <- paste_rows(
          "for (s in 1:{S - 1}) {{",
            "omega_raw_alpha_{y}[s] ~ std_normal();",
          "}}",
          .indent = idt(c(1, 2, 1)),
          .parse = FALSE
        )
      } else {
        lambda_term <- ifelse_(shrinkage, " * lambda[i - 1]", "")
        mtext_omega <- paste_rows(
          "for (s in 1:{S - 1}) {{",
            "omega_raw_alpha_{y}[s, 2] ~ normal(omega_alpha_1_{y}[s], tau_alpha_{y}{lambda_term});",
            "for (i in 3:(D - 1)) {{",
              "omega_raw_alpha_{y}[s, i] ~ normal(omega_raw_alpha_{y}[s, i - 1], tau_alpha_{y}{lambda_term});",
            "}}",
          "}}",
          .indent = idt(c(1, 2, 2, 3, 2, 1)),
          .parse = FALSE
        )
      }
      mtext_intercept <- paste_rows(
        mtext_intercept,
        mtext_omega,
        "tau_alpha_{y} ~ {tau_alpha_prior_distr};",
        .indent = idt(c(0, 0, 1)),
        .parse = FALSE
      )
    }
  }
  if (has_fixed) {
    if (vectorizable_prior(beta_prior_distr)) {
      np <- beta_prior_npars
      dpars_fixed <- paste0("beta_prior_pars_", y, "[, ", 1:np, "]",
                            collapse = ", ")
      mtext_fixed <- "to_vector(beta_{y}) ~ {beta_prior_distr}({dpars_fixed});"
    } else {
      k <- rep(1:K_fixed, S - 1)
      s <- rep(1:(S - 1), each = K)
      mtext_fixed <- "beta_{y}[{k},{s}] ~ {beta_prior_distr};"
    }
  }
  if (has_varying) {
    if (noncentered) {
      if (vectorizable_prior(delta_prior_distr)) {
        np <- delta_prior_npars
        dpars_varying <- paste0("delta_prior_pars_", y, "[, ", 1:np, "]",
                                collapse = ", ")
        # can't convert vector[] to vector
        # mtext_varying <- c(idt(1), "to_vector(to_matrix(a_raw_", i, "[, , 1])') ~ ", d, "(", paste0("delta_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
        # TODO: this works but might give a false warning about missing jacobian which needs to be suppressed?
        mtext_varying <- "to_vector(delta_{y}[1]) ~ {delta_prior_distr}({dpars_varying});"
      } else {
        k <- rep(1:K_fixed, S - 1)
        s <- rep(1:(S - 1), each = K_fixed)
        mtext_varying <- "omega_raw_{y}[{s},{k},1] ~ {delta_prior_distr};"
      }
      mtext_varying <- paste_rows(
        mtext_varying,
        "for (s in 1:{S - 1}) {{",
          "to_vector(omega_raw_{y}[s, ,2:D]) ~ std_normal();",
        "}}",
        .indent = idt(c(0, 1, 2, 1)),
        .parse = FALSE
      )
    } else {
      if (vectorizable_prior(delta_prior_distr)) {
        np <- delta_prior_npars
        dpars_varying <- paste0("delta_prior_pars_", y, "[, ", 1:np, "]",
                                collapse = ", ")
        # can't convert vector[] to vector
        # mtext_varying <- c(idt(1), "to_vector(to_matrix(omega_raw_", i, "[, , 1])') ~ ", d, "(", paste0("delta_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
        # TODO: this works but might give a false warning about missing jacobian which needs to be suppressed?
        mtext_varying <- "to_vector(delta_{y}[1]) ~ {delta_prior_distr}({dpars_varying});"
      } else {
        k <- rep(1:K_fixed, S - 1)
        s <- rep(1:(S - 1), each = K_varying)
        mtext_varying <- "omega_{y}[{s}, {k}, 1] ~ {delta_prior_distr};"
      }
      mtext_varying <- paste_rows(
        mtext_varying,
        "for (s in 1:{S - 1}) {{",
          "for (i in 2:D) {{",
          ifelse_(
           shrinkage,
           "omega_{y}[s, , i] ~ normal(omega_{y}[s, , i - 1], lambda[i - 1] * tau_{y});",
           "omega_{y}[s, , i] ~ normal(omega_{y}[s, , i - 1], tau_{y});"
          ),
          "}}",
        "}}",
        .indent = idt(c(0, 1, 2, 3, 2, 1)),
        .parse = FALSE
      )
    }
    if (vectorizable_prior(tau_prior_distr)) {
      np <- tau_prior_npars
      dpars_tau <- paste0("tau_prior_pars_", y, "[, ", 1:np, "]",
                          collapse = ", ")
      mtext_tau <- "tau_{y} ~ {tau_prior_distr}({dpars_tau});"
    } else {
      mtext_tau <- "tau_{y}[{{{cs(1:K_varying)}}}] ~ {tau_prior_distr};"
    }
  }
  intercept <- glue::glue("zeros_S_{y}")
  if (has_fixed_intercept) intercept <-  glue::glue("append_row(0, alpha_{y})")
  if (has_varying_intercept) intercept <- glue::glue("append_row(0, alpha_{y}[t])")
  likelihood_term <- "{y}[t, {obs}] ~ categorical_logit_glm(X[t][{obs}, {{{cs(J)}}}], {intercept}, append_col(zeros_K_{y}, gamma_{y}));"
  paste_rows(
    mtext_intercept, mtext_fixed, mtext_varying, mtext_tau,
    "{{",
      "matrix[{K}, {S-1}] gamma_{y};",
      "for (t in 1:T) {{",
        onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y};"),
        onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[t];"),
        likelihood_term,
      "}}",
    "}}",
    .indent = idt(c(1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 2, 1))
  )
})

model_lines_gaussian <- quote({
  mtext_def <- eval(model_lines_default)
  d <- sigma_prior_distr
  sigma_term <- "sigma_{y} ~ {d};"
  intercept <- "0"
  if (has_fixed_intercept) intercept <-  glue::glue("alpha_{y}")
  if (has_varying_intercept) intercept <- glue::glue("alpha_{y}[t]")
  likelihood_term <- "{y}[{obs}, t] ~ normal_id_glm(X[t][{obs}, {{{cs(J)}}}], {intercept}, gamma_{y}, sigma_{y});"
  mtext <- paste_rows(
    sigma_term,
    "{{",
      "vector[{K}] gamma_{y};",
      "for (t in 1:T) {{",
        onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y};"),
        onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[t];"),
        likelihood_term,
      "}}",
    "}}",
    .indent = idt(c(1, 1, 2, 2, 3, 3, 3, 2, 1))
  )
  paste_rows(mtext_def, mtext, .parse = FALSE)
})

model_lines_binomial <- quote({
  mtext_def <- eval(model_lines_default)
  fixed_term <- onlyif(has_fixed, glue::glue("X[t][{obs}, {{{cs(J_fixed)}}}] * beta_{y}"))
  varying_term <- onlyif(has_varying, glue::glue("X[t][{obs}, {{{cs(J_varying)}}}] * delta_{y}[t]"))
  plus <- onlyif(has_fixed && has_varying, " + ")
  intercept <- ""
  if (has_fixed_intercept) intercept <-  glue::glue("alpha_{y} + ")
  if (has_varying_intercept) intercept <- glue::glue("alpha_{y}[t] + ")
  likelihood_term <- "{y}[t, {obs}] ~ binomial_logit(trials_{y}[t, {obs}], {intercept}{fixed_term}{plus}{varying_term});"
  mtext <- paste_rows("for (t in 1:T) {{", likelihood_term, "}}",
                      .indent = idt(c(1, 2, 1)))
  paste_rows(mtext_def, mtext, .parse = FALSE)
})

model_lines_bernoulli <- quote({
  mtext_def <- eval(model_lines_default)
  intercept <- "0"
  if (has_fixed_intercept) intercept <-  glue::glue("alpha_{y}")
  if (has_varying_intercept) intercept <- glue::glue("alpha_{y}[t]")
  likelihood_term <- "{y}[t, {obs}] ~ bernoulli_logit_glm(X[t][{obs}, {{{cs(J)}}}], {intercept}, gamma_{y});"
  mtext <- paste_rows(
    "{{",
    "vector[{K}] gamma_{y};",
    "for (t in 1:T) {{",
    onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y};"),
    onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[t];"),
    likelihood_term,
    "}}",
    "}}",
    .indent = idt(c(1, 2, 2, 3, 3, 3, 2, 1))
  )
  paste_rows(mtext_def, mtext, .parse = FALSE)
})

model_lines_poisson <- quote({
  mtext_def <- eval(model_lines_default)
  offset_term <- ifelse_(has_offset, glue::glue("to_vector(offset_{y}[t, {obs}])"), "0")
  intercept <- ""
  if (has_fixed_intercept) intercept <-  glue::glue("alpha_{y} + ")
  if (has_varying_intercept) intercept <- glue::glue("alpha_{y}[t] + ")
  likelihood_term <- "{y}[t, {obs}] ~ poisson_log_glm(X[t][{obs}, {{{cs(J)}}}], {intercept}{offset_term}, gamma_{y});"
  mtext <- paste_rows(
    "{{",
    "vector[{K}] gamma_{y};",
    "for (t in 1:T) {{",
    onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y};"),
    onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[t];"),
    likelihood_term,
    "}}",
    "}}",
    .indent = idt(c(1, 2, 2, 3, 3, 3, 2, 1))
  )
  paste_rows(mtext_def, mtext, .parse = FALSE)
})

model_lines_negbin <- quote({
  mtext_def <- eval(model_lines_default)
  d <- phi_prior_distr
  phi_term <- "phi_{y} ~ {d};"
  offset_term <- ifelse_(has_offset, glue::glue("to_vector(offset_{y}[t, {obs}])"), "0")
  intercept <- ""
  if (has_fixed_intercept) intercept <-  glue::glue("alpha_{y} + ")
  if (has_varying_intercept) intercept <- glue::glue("alpha_{y}[t] + ")
  likelihood_term <- "{y}[t, {obs}] ~ neg_binomial_2_log_glm(X[t][{obs}, {{{cs(J)}}}], {intercept}{offset_term}, gamma_{y}, phi_{y});"
  mtext <- paste_rows(
    phi_term,
    "{{",
    "vector[{K}] gamma_{y};",
    "for (t in 1:T) {{",
    onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y};"),
    onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[t];"),
    likelihood_term,
    "}}",
    "}}",
    .indent = idt(c(1, 1, 2, 2, 3, 3, 3, 2, 1))
  )
  paste_rows(mtext_def, mtext, .parse = FALSE)
})

# Generated quantities block ----------------------------------------------

generated_quantities_lines_default <- quote({
    ""
})

generated_quantities_lines_categorical <- quote({
  eval(generated_quantities_lines_default)
})

generated_quantities_lines_gaussian <- quote({
  eval(generated_quantities_lines_default)
})

generated_quantities_lines_binomial <- quote({
  eval(generated_quantities_lines_default)
})

generated_quantities_lines_bernoulli <- quote({
  eval(generated_quantities_lines_default)
})

generated_quantities_lines_poisson <- quote({
  eval(generated_quantities_lines_default)
})

generated_quantities_lines_negbin <- quote({
  eval(generated_quantities_lines_default)
})

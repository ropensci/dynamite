#' Wrapper to parse Stan model blocks based on the family of the response
#'
#' @param prefix A character string indicating the Stan model block name
#' @param formula A `dynamiteformula` object
#' @param args Channel specific component of `model_vars`
#'   (see [create_blocks()])
#'
#' @noRd
lines_wrap <- function(prefix, formula, args) {
  lines_expr <- paste0(prefix, "_lines_", formula$family)
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
  mtext_fixed <- ""
  mtext_varying <- ""
  mtext_tau <- ""
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
  paste_rows(mtext, mtext_fixed, mtext_varying, mtext_tau, .parse = FALSE)
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

parameters_lines_default <- quote({
  aname <- ifelse_(noncentered, "alpha_raw_", "alpha_")
  paste_rows(
    onlyif(has_fixed, "vector[{K_fixed}] beta_{y};"),
    onlyif(has_varying, "matrix[{K_varying}, D] {aname}{y};"),
    onlyif(has_varying, "vector<lower={lb}>[{K_varying}] tau_{y};"),
    .indent = idt(1)
  )
})

parameters_lines_categorical <- quote({
  aname <- ifelse_(noncentered, "alpha_raw_", "alpha_")
  paste_rows(
    onlyif(has_fixed, "matrix[{K_fixed}, {S - 1}] beta_{y};"),
    onlyif(has_varying, "matrix[{K_varying}, D] {aname}{y}[{S - 1}];"),
    onlyif(has_varying, "vector<lower={lb}>[{K_varying}] tau_{y};"),
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

transformed_parameters_lines_default <- quote({
  mtext <- ""
  mtext_varying_noncentered <- ""
  mtext_varying <- ""
  if (has_varying) {
    mtext <- "vector[{K_varying}] delta_{y}[T];"
    if (noncentered) {
      lambda_term <- ifelse_(shrinkage, " * lambda[i - 1];", ";")
      mtext_varying_noncentered <- paste_rows(
        "matrix[{K_varying}, D] alpha_{y};",
        "alpha_{y}[, 1] = alpha_raw_{y}[, 1];",
        "for (i in 2:D) {{",
          "alpha_{y}[, i] = alpha_{y}[, i - 1] + alpha_raw_{y}[, i] .* tau_{y}{lambda_term}",
        "}}",
        .indent = idt(c(1, 1, 1, 2, 1)),
        .parse = FALSE
      )
    }
    mtext_varying <- paste_rows(
      "for (t in 1:T) {{",
        "delta_{y}[t] = alpha_{y} * Bs[, t];",
      "}}",
      .indent = idt(c(1, 2, 1)),
      .parse = FALSE
    )
  }
  paste_rows(mtext, mtext_varying_noncentered, mtext_varying,
    .indent = idt(c(1, 0, 0))
  )
})

transformed_parameters_lines_categorical <- quote({
  mtext <- ""
  mtext_varying_noncentered <- ""
  mtext_varying <- ""
  if (has_varying) {
    mtext <- "matrix[{K_varying}, {S - 1}] delta_{y}[T];"
    if (noncentered) {
      lambda_term <- ifelse_(shrinkage, " * lambda[i - 1];", ";")
      mtext_varying_noncentered <- paste_rows(
        "matrix[{K_varying}, D] alpha_{y}[{S - 1}];",
        "for (s in 1:{S - 1}) {{",
          "alpha_{y}[s, , 1] = alpha_raw_{y}[s, , 1];",
          "for (i in 2:D) {{",
            "alpha_{y}[s, , i] = alpha_{y}[s, , i - 1] + alpha_raw_{y}[s, , i] .* tau_{y}{lambda_term}",
          "}}",
        "}}",
        .indent = idt(c(1, 1, 2, 2, 3, 2, 1)),
        .parse = FALSE
      )
    }
    mtext_varying <- paste_rows(
      "for (s in 1:{S - 1}) {{",
        "for (t in 1:T) {{",
          "delta_{y}[t, , s] = alpha_{y}[s] * Bs[, t];",
        "}}",
      "}}",
      .indent = idt(c(1, 2, 3, 2, 1)),
      .parse = FALSE
    )
  }
  paste_rows(mtext, mtext_varying_noncentered, mtext_varying,
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

model_lines_default <- quote({
  mtext_fixed <- ""
  mtext_varying <- ""
  mtext_tau <- ""
  if (has_fixed) {
    if (vectorizable_prior(beta_prior_distr)) {
      np <- beta_prior_npars
      dpars_fixed <- paste0("beta_prior_pars_", y, "[, ", 1:np, "]",
                            collapse = ", ")
      mtext_fixed <- "beta_{y} ~ {beta_prior_distr}({dpars_fixed});"
    } else {
      mtext_fixed <- "beta_{y}[k] ~ {beta_prior_distr};"
    }
  }
  if (has_varying) {
    if (noncentered) {
      if (vectorizable_prior(delta_prior_distr)) {
        np <- delta_prior_npars
        dpars_varying <- paste0("delta_prior_pars_", y, "[, ", 1:np, "]",
                                collapse = ", ")
        mtext_varying <- "alpha_raw_{y}[, 1] ~ {delta_prior_distr}({dpars_varying});"
      } else {
        mtext_varying <- "alpha_raw_{y}[{{{cs(1:K_varying)}}}, 1] ~ {delta_prior_distr};"
      }
      mtext_varying <- paste_rows(mtext_varying,
        "to_vector(alpha_raw_{y}[, 2:D]) ~ std_normal();",
        .indent = idt(c(0, 1)),
        .parse = FALSE
      )
    } else {
      if (vectorizable_prior(delta_prior_distr)) {
        np <- delta_prior_npars
        dpars_varying <- paste0("delta_prior_pars_", y, "[, ", 1:np, "]",
                                collapse = ", ")
        mtext_varying <- "alpha_{y}[, 1] ~ {delta_prior_distr}({dpars_varying});"
      } else {
        mtext_varying <- "alpha_{y}[{{{cs(1:K_varying)}}}, 1] ~ {delta_prior_distr};"
      }
      mtext_varying <- paste_rows(
        mtext_varying,
        "for(i in 2:D) {{",
        ifelse_(
          shrinkage,
          "alpha_{y}[, i] ~ normal(alpha_{y}[, i - 1], lambda[i - 1] * tau_{y});",
          "alpha_{y}[, i] ~ normal(alpha_{y}[, i - 1], tau_{y});"
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
  paste_rows(mtext_fixed, mtext_varying, mtext_tau, .indent = idt(c(1, 1, 1)))
})

model_lines_categorical <- quote({
  mtext_fixed <- ""
  mtext_varying <- ""
  mtext_tau <- ""
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
        mtext_varying <- "alpha_raw_{y}[{s},{k},1] ~ {delta_prior_distr};"
      }
      mtext_varying <- paste_rows(
        mtext_varying,
        "for (s in 1:{S - 1}) {{",
          "to_vector(alpha_raw_{y}[s, ,2:D]) ~ std_normal();",
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
        # mtext_varying <- c(idt(1), "to_vector(to_matrix(alpha_raw_", i, "[, , 1])') ~ ", d, "(", paste0("delta_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
        # TODO: this works but might give a false warning about missing jacobian which needs to be suppressed?
        mtext_varying <- "to_vector(delta_{y}[1]) ~ {delta_prior_distr}({dpars_varying});"
      } else {
        k <- rep(1:K_fixed, S - 1)
        s <- rep(1:(S - 1), each = K_varying)
        mtext_varying <- "alpha_{y}[{s}, {k}, 1] ~ {delta_prior_distr};"
      }
      mtext_varying <- paste_rows(
        mtext_varying,
        "for (s in 1:{S - 1}) {{",
          "for(i in 2:D) {{",
          ifelse_(
           shrinkage,
           "alpha_{y}[s, , i] ~ normal(alpha_{y}[s, , i - 1], lambda[i - 1] * tau_{y});",
           "alpha_{y}[s, , i] ~ normal(alpha_{y}[s, , i - 1], tau_{y});"
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
  likelihood_term <- "{y}[t, {obs}] ~ categorical_logit_glm(X[t][{obs}, {{{cs(J)}}}], zeros_S_{y}, append_col(zeros_K_{y}, gamma_{y}));"
  paste_rows(
    mtext_fixed, mtext_varying, mtext_tau,
    "{{",
      "matrix[{K}, {S-1}] gamma_{y};",
      "for (t in 1:T) {{",
        onlyif(has_fixed, "gamma_{y}[{{{cs(L_fixed)}}}] = beta_{y};"),
        onlyif(has_varying, "gamma_{y}[{{{cs(L_varying)}}}] = delta_{y}[t];"),
        likelihood_term,
      "}}",
    "}}",
    .indent = idt(c(1, 1, 1, 1, 2, 2, 3, 3, 3, 2, 1))
  )
})

model_lines_gaussian <- quote({
  mtext_def <- eval(model_lines_default)
  d <- sigma_prior_distr
  sigma_term <- "sigma_{y} ~ {d};"
  likelihood_term <- "{y}[{obs}, t] ~ normal_id_glm(X[t][{obs}, {{{cs(J)}}}], 0, gamma_{y}, sigma_{y});"
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
  likelihood_term <- "{y}[t, {obs}] ~ binomial_logit(trials_{y}[t, {obs}], {fixed_term}{plus}{varying_term});"
  mtext <- paste_rows("for (t in 1:T) {{", likelihood_term, "}}",
                      .indent = idt(c(1, 2, 1)))
  paste_rows(mtext_def, mtext, .parse = FALSE)
})

model_lines_bernoulli <- quote({
  mtext_def <- eval(model_lines_default)
  likelihood_term <- "{y}[t, {obs}] ~ bernoulli_logit_glm(X[t][{obs}, {{{cs(J)}}}], 0, gamma_{y});"
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
  likelihood_term <- "{y}[t, {obs}] ~ poisson_log_glm(X[t][{obs}, {{{cs(J)}}}], {offset_term}, gamma_{y});"
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
  likelihood_term <- "{y}[t, {obs}] ~ neg_binomial_2_log_glm(X[t][{obs}, {{{cs(J)}}}], 0, gamma_{y}, phi_{y});"
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

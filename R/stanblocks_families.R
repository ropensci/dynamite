# NOTE: Added this wrapper, so that *_lines_* functions can have different arguments (including parameters and model)
lines_wrap <- function(prefix, formula, args) {
    lines_expr <- paste0(prefix, "_lines_", formula$family)
    lines_env <- list2env(args)
    eval(eval(as.name(lines_expr)), envir = lines_env)
}

vectorizable_prior <- function(x) length(x) == 1 && !grepl("\\(", x)

# For data block
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

# For transformed data block
transformed_data_lines_default <- quote({
    mtext_fixed <- ""
    mtext_varying <- ""
    mtext_tau <- ""
    if (write_beta_fixed) {
        i <- rep(1:K_fixed, beta_fixed_prior_npars)
        j <- rep(1:beta_fixed_prior_npars, each = K_fixed)
        mtext_fixed <- paste_rows(
            "matrix[{K_fixed}, {beta_fixed_prior_npars}] beta_fixed_prior_pars_{y};",
            "beta_fixed_prior_pars_{y}[{i},{j}] = {beta_fixed_prior_pars};",
            .indent = idt(1)
        )
    }
    if (write_beta_varying) {
        i <- rep(1:K_varying, beta_varying_prior_npars)
        j <- rep(1:beta_varying_prior_npars, each = K_varying)
        mtext_varying <- paste_rows(
            "matrix[{K_varying}, {beta_varying_prior_npars}] beta_varying_prior_pars_{y};",
            "beta_varying_prior_pars_{y}[{i},{j}] = {beta_varying_prior_pars};",
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
    if (write_beta_fixed) {
        k <- (K_fixed * (S - 1))
        i <- rep(1:k, beta_fixed_prior_npars)
        j <- rep(1:beta_fixed_prior_npars, each = k)
        mtext_fixed <- paste_rows(
            "matrix[{k}, {beta_fixed_prior_npars}] beta_fixed_prior_pars_{y};",
            "beta_fixed_prior_pars_{y}[{i},{j}] = {beta_fixed_prior_pars};",
            .indent = idt(1)
        )
    }
    if (write_beta_varying) {
        k <- (K_varying * (S - 1))
        i <- rep(1:k, beta_varying_prior_npars)
        j <- rep(1:beta_varying_prior_npars, each = k)
        mtext_varying <- paste_rows(
            "matrix[{k}, {beta_varying_prior_npars}] beta_varying_prior_pars_{y};",
            "beta_varying_prior_pars_{y}[{i},{j}] = {beta_varying_prior_pars};",
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

# For parameters block
parameters_lines_default <- quote({
    aname <- ifelse_(noncentered, "a_raw_", "a_")
    paste_rows(
        onlyif(has_fixed,   "vector[{K_fixed}] beta_fixed_{y};"),
        onlyif(has_varying, "matrix[{K_varying}, D] {aname}{y};"),
        onlyif(has_varying, "vector<lower={lb}>[{K_varying}] tau_{y};"),
        .indent = idt(1)
    )
})

parameters_lines_categorical <- quote({
    aname <- ifelse_(noncentered, "a_raw_", "a_")
    paste_rows(
        onlyif(has_fixed,   "matrix[{K_fixed}, {S - 1}] beta_fixed_{y};"),
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

# For transformed parameters block
transformed_parameters_lines_default <- quote({
    mtext_fixed <- ""
    mtext_varying <- ""
    mtext_varying_noncentered <- ""
    mtext <- "vector[{K}] beta_{y}[T];"
    if (has_fixed) {
        mtext_fixed <- "beta_{y}[1:T, {{{cs(L_fixed)}}}] = rep_array(beta_fixed_{y}, T);"
    }
    if (has_varying) {
        if (noncentered) {
            lambda_term <- ifelse_(shrinkage, " * lambda[i - 1];", ";")
            mtext_varying_noncentered <- paste_rows(
                "matrix[{K_varying}, D] a_{y};",
                "a_{y}[, 1] = a_raw_{y}[, 1];",
                "for (i in 2:D) {{",
                    "a_{y}[, i] = a_{y}[, i - 1] + a_raw_{y}[, i] .* tau_{y}{lambda_term}",
                "}}",
                .indent = idt(c(1, 1, 1, 2, 1)),
                .parse = FALSE
            )
        }
        mtext_varying <- paste_rows(
            "for (t in 1:T) {{",
                "beta_{y}[t, {{{cs(L_varying)}}}] = a_{y} * Bs[, t];",
            "}}",
            .indent = idt(c(1, 2, 1)),
            .parse = FALSE
        )
    }
    paste_rows(mtext, mtext_fixed, mtext_varying, mtext_varying_noncentered,
               .indent = idt(c(1, 1, 0, 0)))
})

transformed_parameters_lines_categorical <- quote({
    mtext_fixed <- ""
    mtext_varying <- ""
    mtext_varying_noncentered <- ""
    mtext <- "matrix[{K}, {S - 1}] beta_{y}[T];"
    if (has_fixed) {
        mtext_fixed <- "beta_{y}[1:T, {{{cs(L_fixed)}}}, {{{cs(1:(S - 1))}}}] = rep_array(beta_fixed_{y}, T);"
    }
    if (has_varying) {
        if (noncentered) {
            lambda_term <- ifelse_(shrinkage, " * lambda[i - 1];", ";")
            mtext_varying_noncentered <- paste_rows(
                "matrix[{K_varying}, D] a_{y}[{S - 1}];",
                "for (s in 1:{S - 1}) {{",
                    "a_{y}[s, , 1] = a_raw_{y}[s, , 1];",
                    "for (i in 2:D) {{",
                        "a_{y}[s, , i] = a_{y}[s, , i - 1] + a_raw_{y}[s, , i] .* tau_{y}{lambda_term}",
                    "}}",
                "}}",
                .indent = idt(c(1, 1, 2, 2, 3, 2, 1)),
                .parse = FALSE
            )
        }
        mtext_varying <- paste_rows(
            "for (s in 1:{S - 1}) {{",
                "for (t in 1:T) {{",
                    "beta_{y}[t, {{{cs(L_varying)}}}, s] = a_{y}[s] * Bs[, t];",
                 "}}",
            "}}",
            .indent = idt(c(1, 2, 3, 2, 1)),
            .parse = FALSE
        )
    }
    paste_rows(mtext, mtext_fixed, mtext_varying, mtext_varying_noncentered,
               .indent = idt(c(1, 1, 0, 0)))
})

transformed_parameters_lines_binomial <- quote({
    mtext <- ""
    mtext_noncentered <- ""
    if (has_varying) {
        if (noncentered) {
            lambda_term <- ifelse_(shrinkage, " * lambda[i - 1];", ";")
            mtext_noncentered <- paste_rows(
                "matrix[{K_varying}, D] a_{y};",
                "a_{y}[, 1] = a_raw_{y}[, 1];",
                "for (i in 2:D) {{",
                     "a_{y}[, i] = a_{y}[, i - 1] + a_raw_{y}[, i] .* tau_{y}{lambda_term}",
                "}}",
                .indent = idt(c(1, 1, 1, 2, 1))
            )
        }
        mtext <- paste_rows(
            "vector[{K_varying}] beta_varying_{y}[T];",
            "for (t in 1:T) {{",
                "beta_varying_{y}[t] = a_{y} * Bs[, t];",
            "}}",
            .indent = idt(c(1, 1, 2, 1))
        )
    }
    paste_rows(mtext, mtext_noncentered, .parse = FALSE)
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

# For model block
model_lines_default <- quote({
    mtext_fixed <- ""
    mtext_varying <- ""
    mtext_tau <- ""
    if (has_fixed) {
        if (vectorizable_prior(beta_fixed_prior_distr)) {
            np <- beta_fixed_prior_npars
            dpars_fixed <- paste0("beta_fixed_prior_pars_", y, "[, ", 1:np, "]", collapse = ", ")
            mtext_fixed <- "beta_fixed_{y} ~ {beta_fixed_prior_distr}({dpars_fixed});"
        } else {
            mtext_fixed <- "beta_fixed_{y}[k] ~ {beta_fixed_prior_distr};"
        }
    }
    if (has_varying) {
        if (noncentered) {
            if (vectorizable_prior(beta_varying_prior_distr)) {
                np <- beta_varying_prior_npars
                dpars_varying <- paste0("beta_varying_prior_pars_", y, "[, ", 1:np, "]", collapse = ", ")
                mtext_varying <- "a_raw_{y}[, 1] ~ {beta_varying_prior_distr}({dpars_varying});"
            } else {
                mtext_varying <- "a_raw_{y}[{{{cs(1:K_varying)}}}, 1] ~ {beta_varying_prior_distr};"
            }
            mtext_varying <- paste_rows(mtext_varying,
                "to_vector(a_raw_{y}[, 2:D]) ~ std_normal();",
                .indent = idt(c(1, 1)),
                .parse = FALSE
            )
        } else {
            if (vectorizable_prior(beta_varying_prior_distr)) {
                np <- beta_varying_prior_npars
                dpars_varying <- paste0("beta_varying_prior_pars_", y, "[, ", 1:np, "]", collapse = ", ")
                mtext_varying <- "a_{y}[, 1] ~ {beta_varying_prior_distr}({dpars_varying});"
            } else {
                mtext_varying <- "a_{y}[{{{cs(1:K_varying)}}}, 1] ~ {beta_varying_prior_distr};"
            }
            mtext_varying <- paste_rows(
                mtext_varying,
                "for(i in 2:D) {{",
                ifelse_(shrinkage,
                    "a_{y}[, i] ~ normal(a_{y}[, i - 1], lambda[i - 1] * tau_{y});",
                    "a_{y}[, i] ~ normal(a_{y}[, i - 1], tau_{y});"
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
        if (vectorizable_prior(beta_fixed_prior_distr)) {
            np <- beta_fixed_prior_npars
            dpars_fixed <- paste0("beta_fixed_prior_pars_", y, "[, ", 1:np, "]", collapse = ", ")
            mtext_fixed <- "to_vector(beta_fixed_{y}) ~ {beta_fixed_prior_distr}({dpars_fixed});"
        } else {
            k <- rep(1:K_fixed, S - 1)
            s <- rep(1:(S - 1), each = K)
            mtext_fixed <- "beta_fixed_{y}[{k},{s}] ~ {beta_fixed_prior_distr};"
        }
    }
    if (has_varying) {
        if (noncentered) {
            if (vectorizable_prior(beta_varying_prior_distr)) {
                np <- beta_varying_prior_npars
                dpars_varying <- paste0("beta_varying_prior_pars_", y, "[, ", 1:np, "]", collapse = ", ")
                # can't convert vector[] to vector
                # mtext_varying <- c(idt(1), "to_vector(to_matrix(a_raw_", i, "[, , 1])') ~ ", d, "(", paste0("beta_varying_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
                # TODO: this works but might give a false warning about missing jacobian which needs to be suppressed?
                mtext_varying <- "to_vector(beta_{y}[1, {{{cs(L_varying)}}}, ]) ~ {beta_varying_prior_distr}({dpars_varying});"
            } else {
                k <- rep(1:K_fixed, S - 1)
                s <- rep(1:(S - 1), each = K_fixed)
                mtext_varying <- "a_raw_{y}[{s},{k},1] ~ {beta_varying_prior_distr};"
            }
            mtext_varying <- paste_rows(
                mtext_varying,
                "for (s in 1:{S - 1}) {{",
                    "to_vector(a_raw_{y}[s, ,2:D]) ~ std_normal();",
                "}}",
                .indent = idt(c(1, 1, 2, 1)),
                .parse = FALSE
            )
        } else {
            if (vectorizable_prior(beta_varying_prior_distr)) {
                np <- beta_varying_prior_npars
                dpars_varying <- paste0("beta_varying_prior_pars_", y, "[, ", 1:np, "]", collapse = ", ")
                # can't convert vector[] to vector
                # mtext_varying <- c(idt(1), "to_vector(to_matrix(a_raw_", i, "[, , 1])') ~ ", d, "(", paste0("beta_varying_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
                # TODO: this works but might give a false warning about missing jacobian which needs to be suppressed?
                mtext_varying <- "to_vector(beta_{y}[1, {{{cs(L_varying)}}}, ]) ~ {beta_varying_prior_distr}({dpars_varying});"
            } else {
                k <- rep(1:K_fixed, S - 1)
                s <- rep(1:(S - 1), each = K_varying)
                mtext_varying <- "a_{y}[{s}, {k}, 1] ~ {beta_varying_prior_distr};"
            }
            mtext_varying <- paste_rows(
                mtext_varying,
                "for (s in 1:{S - 1}) {{",
                    "for(i in 2:D) {{",
                    ifelse_(shrinkage,
                        "a_{y}[s, , i] ~ normal(a_{y}[s, , i - 1], lambda[i - 1] * tau_{y});",
                        "a_{y}[s, , i] ~ normal(a_{y}[s, , i - 1], tau_{y});"
                    ),
                    "}}",
                "}}",
                .indent = idt(c(1, 1, 2, 3, 2, 1)),
                .parse = FALSE
            )
        }
        if (vectorizable_prior(tau_prior_distr)) {
            np <- tau_prior_npars
            dpars_tau <-  paste0("tau_prior_pars_", y, "[, ", 1:np, "]", collapse = ", ")
            mtext_tau <- "tau_{y} ~ {tau_prior_distr}({dpars_tau});"
        } else {
            mtext_tau <- "tau_{y}[{{{cs(1:K_varying)}}}] ~ {tau_prior_distr};"
        }
    }
    likelihood_term <- "{y}[t, {obs}] ~ categorical_logit_glm(X[t][{obs}, {{{cs(J)}}}], zeros_S_{y}, append_col(zeros_K_{y}, beta_{y}[t]));"
    paste_rows(mtext_fixed, mtext_varying, mtext_tau, "for (t in 1:T) {{", likelihood_term, "}}",
               .indent = idt(c(1, 0, 1, 1, 2, 1)))
})

model_lines_gaussian <- quote({
    mtext_def <- eval(model_lines_default)
    d <- sigma_prior_distr
    sigma_term <- "sigma_{y} ~ {d};"
    likelihood_term <- "{y}[{obs}, t] ~ normal_id_glm(X[t][{obs}, {{{cs(J)}}}], 0, beta_{y}[t], sigma_{y});"
    mtext <- paste_rows(sigma_term, "for (t in 1:T) {{", likelihood_term, "}}", .indent = idt(c(1, 1, 2, 1)))
    paste_rows(mtext_def, mtext, .parse = FALSE)
})

model_lines_binomial <- quote({
    mtext_def <- eval(model_lines_default)
    fixed_term <- onlyif(has_fixed, glue::glue("X[t][{obs}, {{{cs(J_fixed)}}}] * beta_fixed_{y}"))
    varying_term <- onlyif(has_varying, glue::glue("X[t][{obs}, {{{cs(J_varying)}}}] * beta_varying_{y}[t]"))
    plus <- onlyif(has_fixed && has_varying, " + ")
    likelihood_term <- "{y}[t, {obs}] ~ binomial_logit(trials_{y}[t, {obs}], {fixed_term}{plus}{varying_term});"
    mtext <- paste_rows("for (t in 1:T) {{", likelihood_term, "}}", .indent = idt(c(1, 2, 1)))
    paste_rows(mtext_def, mtext, .parse = FALSE)
})

model_lines_bernoulli <- quote({
    mtext_def <- eval(model_lines_default)
    likelihood_term <- "{y}[t, {obs}] ~ bernoulli_logit_glm(X[t][{obs}, {{{cs(J)}}}], 0, beta_{y}[t]);"
    mtext <- paste_rows(mtext, "for (t in 1:T) {{", likelihood_term, "}}", .indent = idt(c(1, 2, 1)))
    paste_rows(mtext_def, mtext, .parse = FALSE)
})

model_lines_poisson <-  quote({
    mtext_def <- eval(model_lines_default)
    offset_term <- ifelse_(has_offset, glue::glue("to_vector(offset_{y}[t, {obs}])"), "0")
    likelihood_term <- "{y}[t, {obs}] ~ poisson_log_glm(X[t][{obs}, {{{cs(J)}}}], {offset_term}, beta_{y}[t]);"
    mtext <- paste_rows("for (t in 1:T) {{", likelihood_term, "}}", .indent = idt(c(1, 2, 1)))
    paste_rows(mtext_def, mtext, .parse = FALSE)
})

model_lines_negbin <- quote({
    mtext_def <- eval(model_lines_default)
    d <- phi_prior_distr
    phi_term <- "phi_{y} ~ {d};"
    likelihood_term <- "{y}[t, {obs}] ~ neg_binomial_2_log_glm(X[t][{obs}, {{{cs(J)}}}], 0, beta_{y}[t], phi_{y});"
    mtext <- paste_rows(phi_term, "for (t in 1:T) {{", likelihood_term, "}}", .indent = idt(c(1, 1, 2, 1)))
    paste_rows(mtext_def, mtext, .parse = FALSE)
})

# For generated quantities block
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
    mtext <- "vector[{K}] beta_{y}[T];"
    mtext_fixed <- ""
    mtext_varying <- ""
    if (has_fixed) {
        mtext_fixed <- "beta_{y}[1:T, {{{cs(L_fixed)}}}] = rep_array(beta_fixed_{y}, T);"
    }
    if (has_varying) {
        mtext_varying <- "beta_{y}[1:T, {{{cs(L_varying)}}}] = beta_varying_{y};"
    }
    paste_rows(mtext, mtext_fixed, mtext_varying, .indent = idt(c(1, 1, 1)))
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

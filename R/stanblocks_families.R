# NOTE: Added this wrapper, so that *_lines_* functions can have different arguments (including parameters and model)
lines_wrap <- function(prefix, formula, args) {
    lines_fun <- paste0(prefix, "_lines_", formula$family)
    # has_args <- names(args) %in% names(formals(get(lines_fun)))
    # args[!has_args] <- NULL
    do.call(lines_fun, args)
}

# For data block
data_lines_default <- function(i, idt, has_fixed, has_varying, ...) {
    paste_rows(onlyif(has_fixed, c(idt(1), "vector[K_fixed_", i, "] beta_fixed_prior_mean_", i, ";")),
               onlyif(has_fixed, c(idt(1), "vector[K_fixed_", i, "] beta_fixed_prior_sd_", i, ";")),
               onlyif(has_varying, c(idt(1), "vector[K_varying_", i, "] beta_varying_prior_mean_", i, ";")),
               onlyif(has_varying, c(idt(1), "vector[K_varying_", i, "] beta_varying_prior_sd_", i, ";")))

}

data_lines_categorical <- function(i, idt, has_fixed, has_varying, ...) {
    paste_rows(c(idt(1), "int<lower=1> ", i, "[T, N];"),
               c(idt(1), "int<lower=0> S_", i, ";"),
               onlyif(has_fixed, c(idt(1), "matrix[K_fixed_", i, ", S_", i, " - 1] beta_fixed_prior_mean_", i, ";")),
               onlyif(has_fixed, c(idt(1), "matrix[K_fixed_", i, ", S_", i, " - 1] beta_fixed_prior_sd_", i, ";")),
               onlyif(has_varying, c(idt(1), "matrix[K_varying_", i, ", S_", i, " - 1] beta_varying_prior_sd_", i, ";")),
               onlyif(has_varying, c(idt(1), "matrix[K_varying_", i, ", S_", i, " - 1] beta_varying_prior_sd_", i, ";")))
}

data_lines_gaussian <- function(i, idt, ...) {
    paste_rows(c(idt(1), "real ", i, "[T, N];"),
               c(idt(1), "real<lower=0> sigma_scale_", i, ";"),
               data_lines_default(i, idt, ...))
}

data_lines_binomial <- function(i, idt, ...) {
    paste_rows(c(idt(1), "int<lower=0> ", i, "[T, N];"),
               c(idt(1), "int<lower=1> trials_", i, "[T, N];"),
               data_lines_default(i, idt, ...))
}

data_lines_bernoulli <- function(i, idt, ...) {
    paste_rows(c(idt(1), "int<lower=0,upper=1> ", i, "[T, N];"),
               data_lines_default(i, idt, ...))
}

data_lines_poisson <- function(i, idt, has_offset, ...) {
    paste_rows(c(idt(1), "int<lower=0> ", i, "[T, N];"),
               onlyif(has_offset, c(idt(1), "real offset_", i, "[T, N];")),
               data_lines_default(i, idt, ...))
}

data_lines_negbin <- function(i, idt, ...) {
    paste_rows(c(idt(1), "int<lower=0> ", i, "[T, N];"),
               c(idt(1), "real<lower=0> phi_scale_", i),
               data_lines_default(i, idt, ...))
}

# For parameters block
parameters_lines_default <- function(i, lb, idt, has_fixed, has_varying, ...) {
    paste_rows(onlyif(has_fixed,   c(idt(1), "vector[K_fixed_", i, "] beta_fixed_", i, ";")),
               onlyif(has_varying, c(idt(1), "row_vector[D] a_", i, "[K_varying_", i, "];")),
               onlyif(has_varying, c(idt(1), "vector<lower=", lb, ">[K_varying_", i, "] tau_", i, ";")))
}

parameters_lines_categorical <- function(i, lb, idt, has_fixed, has_varying, ...) {
    paste_rows(onlyif(has_fixed,   c(idt(1), "matrix[K_fixed_", i, ", S_", i, " - 1] beta_fixed_", i, ";")),
               onlyif(has_varying, c(idt(1), "row_vector[D] a_", i, "[S_", i, " - 1, K_varying_", i, "];")),
               onlyif(has_varying, c(idt(1), "vector<lower=", lb, ">[K_varying_", i, "] tau_", i, ";")))
}

parameters_lines_gaussian <- function(i, lb, idt, ...) {
    paste_rows(parameters_lines_default(i, lb, idt, ...),
               c(idt(1), "real<lower=0> sigma_", i, ";"))
}

parameters_lines_binomial <- function(...) {
    parameters_lines_default(...)
}

parameters_lines_bernoulli <- function(...) {
    parameters_lines_default(...)
}

parameters_lines_poisson <- function(...) {
    parameters_lines_default(...)
}

parameters_lines_negbin <- function(i, lb, idt, ...) {
    paste_rows(parameters_lines_default(i, lb, idt, ...),
               c(idt(1), "real<lower=0 phi_", i))
}

# For transformed parameters block
transformed_parameters_lines_default <- function(i, idt, has_varying, ...) {
    mtext <- ""
    if (has_varying) {
        mtext <- paste_rows(
            c(idt(1), "vector[K_varying_", i, "] beta_varying_", i, "[T];"),
            c(idt(1), "for (k in 1:K_varying_", i, ") {"),
            c(idt(2),     "for (t in 1:T) {"),
            c(idt(3),         "beta_varying_", i, "[t, k] = a_", i, "[k] * Bs[, t];"),
            c(idt(2),     "}"),
            c(idt(1), "}")
        )
    }
    mtext
    #if (attr(formula, "splines")$noncentered) { # TODO: channel-wise?
    #    c(a_terms) <- paste0(idt(1), "row_vector[D] a_", y, "[K_", y, "];")
    #}
}

transformed_parameters_lines_categorical <- function(i, idt, has_fixed, has_varying, ...) {
    mtext_fixed <- ""
    mtext_varying <- ""
    mtext <- paste0(idt(1), "matrix[K_", i, ", S_", i, " - 1] beta_", i, "[T];")
    if (has_fixed) {
        mtext_fixed <- paste_rows(
            c(idt(1), "beta_", i, "[1:T, L_fixed_", i, ", 1:(S_", i, " - 1)] = rep_array(beta_fixed_", i, ", T);")
        )
    }
    #if (attr(formula, "splines")$noncentered) { # TODO: channel-wise?
    #    c(a_terms) <- paste0(idt(1), "row_vector[D] a_", y, "[S_", y, " - 1, K_", y, "];")
    #}
    # varying
    if (has_varying) {
        mtext_varying <- paste_rows(
            c(idt(1), "for (s in 1:(S_", i, " - 1)) {"),
            c(idt(2),     "for (k in 1:K_varying_", i, ") {"),
            c(idt(3),         "for (t in 1:T) {"),
            c(idt(4),             "beta_", i, "[t, L_varying_", i, "[k], s] = a_", i, "[s, k] * Bs[, t];"),
            c(idt(3),          "}"),
            c(idt(2),     "}"),
            c(idt(1), "}")
        )
    }
    paste_rows(mtext, mtext_fixed, mtext_varying)
}

transformed_parameters_lines_gaussian <- function(...) {
    transformed_parameters_lines_default(...)
}

transformed_parameters_lines_binomial <- function(...) {
    transformed_parameters_lines_default(...)
}

transformed_parameters_lines_bernoulli <- function(i, idt, has_fixed, has_varying, ...) {
    mtext_fixed <- ""
    mtext_varying <- ""
    mtext <- paste0(idt(1), "vector[K_", i, "] beta_", i, "[T];")
    if (has_fixed) {
        mtext_fixed <- paste0(idt(1), "beta_", i, "[1:T, L_fixed_", i, "] = rep_array(beta_fixed_", i, ", T);")
    }
    if (has_varying) {
        mtext_varying <- paste_rows(
            c(idt(1), "for (k in 1:K_varying_", i, ") {"),
            c(idt(2),     "for (t in 1:T) {"),
            c(idt(3),         "beta_", i, "[t, L_varying_", i, "[k]] = a_", i, "[k] * Bs[, t];"),
            c(idt(2),     "}"),
            c(idt(1), "}")
        )
    }
    paste_rows(mtext, mtext_fixed, mtext_varying)
}

transformed_parameters_lines_poisson <- function(...) {
    transformed_parameters_lines_bernoulli(...)
}

# For model block
model_lines_default <- function(i, idt, shrinkage, noncentered, has_fixed, has_varying, ...) {
    mtext_fixed <- ""
    mtext_varying <- ""
    if (has_fixed) {
        mtext_fixed <- paste_rows(
            c(idt(1), "for (k in 1:K_fixed_", i, ") {"),
            c(idt(2),     "beta_fixed_", i, "[k] ~ normal(beta_fixed_prior_mean_", i, "[k], beta_fixed_prior_sd_", i, "[k]);"),
            c(idt(1), "}")
        )
    }
    if (has_varying) {
        if (noncentered) {
            mtext_varying <- paste_rows(
                c(idt(1), "for (k in 1:K_varying_", i, ") {"),
                c(idt(2),      "a_raw_", i, "[k] ~ std_normal();"),
                c(idt(1), "}")
            )
        } else {
            # Prior for the first a (beta) is always normal given the RW prior
            mtext_varying <- paste_rows(
                c(idt(1), "for (k in 1:K_varying_", i, ") {"),
                c(idt(2),      "a_", i, "[k, 1] ~ normal(beta_varying_prior_mean_", i, "[k], beta_varying_prior_sd_", i, "[k]);"),
                c(idt(2),      "for(i in 2:D) {"),
                if (shrinkage) {
                    c(idt(3),      "a_", i, "[k, i] ~ normal(a_", i,"[k, i - 1], lambda[i - 1] * tau_", i, "[k]);")
                } else {
                    c(idt(3),      "a_", i, "[k, i] ~ normal(a_", i,"[k, i - 1], tau_", i, "[k]);")
                },
                c(idt(2),      "}"),
                c(idt(1), "}")
            )
        }
        mtext_varying <- paste_rows(mtext_varying, c(idt(1), "tau_", i, " ~ normal(0, 1);"))
    }
    paste_rows(mtext_fixed, mtext_varying)
}

model_lines_categorical <- function(i, idt, shrinkage, noncentered, has_fixed, has_varying, ...) {
    mtext_fixed <- ""
    mtext_varying <- ""
    prior_term <- ""
    if (has_fixed) {
        mtext_fixed <- paste_rows(
            c(idt(2), "for (k in 1:K_fixed_", i, ") {"),
            c(idt(2),     "for (s in 1:(S_", i, " - 1)) {"),
            c(idt(3),          "beta_fixed_", i, "[k, s] ~ normal(beta_fixed_prior_mean_", i, "[k], beta_fixed_prior_sd_", i, "[k]);"),
            c(idt(2),     "}"),
            c(idt(1), "}")
        )
    }
    if (has_varying) {
        if (noncentered) {
            # fixed given noncentered parameterisation
            mtext_varying <- paste_rows(
                c(idt(1), "for (s in 1:(S_", i, " - 1)) {"),
                c(idt(2),     "for (k in 1:K_varying_", i, ") {"),
                c(idt(3),          "a_raw_", i, "[s, k] ~ std_normal();"),
                c(idt(2),     "}"),
                c(idt(1), "}")
            )
        } else {
            # Prior for the first a (beta) is always normal given the RW prior
            mtext_varying <- paste_rows(
                c(idt(1), "for (s in 1:(S_", i, " - 1)) {"),
                c(idt(2),     "for (k in 1:K_varying_", i, ") {"),
                c(idt(3),         "a_", i, "[s, k, 1] ~ normal(beta_varying_prior_mean_", i, "[k, s], beta_varying_prior_sd_", i, "[k, s]);"),
                c(idt(3),         "for(i in 2:D) {"),
                if (shrinkage) {
                    c(idt(4),         "a_", i, "[s, k, i] ~ normal(a_", i,"[s, k, i - 1], lambda[i - 1] * tau_", i, "[k]);")
                } else {
                    c(idt(4),         "a_", i, "[s, k, i] ~ normal(a_", i,"[s, k, i - 1], tau_", i, "[k]);")
                },
                c(idt(3),         "}"),
                c(idt(2),     "}"),
                c(idt(1), "}")
            )
        }
        prior_term <- c(idt(1), "to_vector(tau_", i, ") ~ normal(0, 1);")
    }
    likelihood_term <-
        c(idt(2), i, "[t] ~ categorical_logit_glm(X[t][,J_", i, "], zeros_S_", i, ", append_col(beta_", i, "[t], zeros_K_", i, "));")

    paste_rows(mtext_fixed, mtext_varying, prior_term, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_gaussian <- function(i, idt, has_fixed, has_varying, ...) {
    mtext <- model_lines_default(i, idt, has_fixed, has_varying, ...)
    # TODO user-defined prior for sigma
    sigma_term <- c(idt(1), "sigma_", i, " ~ exponential(sigma_scale_", i, ");")
    fixed_term <- onlyif(has_fixed, paste0("X[t][,J_fixed_", i, "] * beta_fixed_", i))
    varying_term <- onlyif(has_varying, paste0("X[t][,J_varying_", i, "] * beta_varying_", i, "[t]"))
    plus <- onlyif(has_fixed && has_varying, " + ")
    likelihood_term <-  c(idt(2), i, "[t] ~ normal(", fixed_term, plus, varying_term, ", sigma_", i, ");")
    paste_rows(mtext, sigma_term, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_binomial <- function(i, idt, has_fixed, has_varying, ...) {
    mtext <- model_lines_default(i, idt, has_fixed, has_varying, ...)
    fixed_term <- onlyif(has_fixed, paste0("X[t][,J_fixed_", i, "] * beta_fixed_", i))
    varying_term <- onlyif(has_varying, paste0("X[t][,J_varying_", i, "] * beta_varying_", i, "[t]"))
    plus <- onlyif(has_fixed && has_varying, " + ")
    likelihood_term <- c(idt(2), i, "[t] ~ binomial_logit(trials_", i, "[t], ", fixed_term, plus, varying_term, ");")
    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_bernoulli <- function(i, idt, ...) {
    mtext <- model_lines_default(i, idt, ...)
    likelihood_term <- c(idt(2), i, "[t] ~ bernoulli_logit_glm(X[t][,J_", i, "], 0, beta_", i, "[t]);")
    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_poisson <- function(i, idt, has_offset, ...) {
    mtext <- model_lines_default(i, idt, ...)
    offset_term <- ifelse_(has_offset, paste0("to_vector(offset_", i, "[t])"), "0")
    likelihood_term <- c(idt(2), i, "[t] ~ poisson_log_glm(X[t][,J_", i, "], ", offset_term, ", beta_", i, "[t]);")
    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_negbin <- function(i, idt, ...) {
    mtext <- model_lines_default(i, idt, ...)
    # TODO user-defined prior for phi
    mtext <- paste_rows(mtext, c(idt(1), "phi_", i, " ~ exponential(phi_scale_", i, ");"))
    likelihood_term <- c(idt(2), i, "[t] ~ neg_binomial_2_log_glm(X[t][,J_", i, "], 0, beta_", i, "[t], ", "phi_", i, ");")
    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

# For generated quantities block
generated_quantities_lines_default <- function(...) {
    ""
}

generated_quantities_lines_categorical <- function(...) {
    generated_quantities_lines_default(...)
}

generated_quantities_lines_gaussian <- function(i, idt, has_fixed, has_varying, ...) {
    mtext <- paste0(idt(1), "vector[K_", i, "] beta_", i, "[T];")
    mtext_fixed <- ""
    mtext_varying <- ""
    if (has_fixed) {
        mtext_fixed <- paste_rows(
            c(idt(1), "beta_", i, "[1:T, L_fixed_", i, "] = rep_array(beta_fixed_", i, ", T);")
        )
    }
    if (has_varying) {
        mtext_varying <- paste_rows(
            c(idt(1), "beta_", i, "[1:T, L_varying_", i, "[1:K_varying_", i, "]] = beta_varying_", i, ";")
            #c(idt(1), "for (k in 1:K_varying_", i, ") {"),
            #c(idt(2),     "for (t in 1:T) {"),
            #c(idt(3),         "beta_", i, "[t, L_varying_", i, "[k]] = beta_varying_", i, "[t, k];"),
            #c(idt(2),     "}"),
            #c(idt(1), "}")
        )
    }
    paste_rows(mtext, mtext_fixed, mtext_varying)
}

generated_quantities_lines_binomial <- function(...) {
    generated_quantities_lines_gaussian(...)
}

generated_quantities_lines_bernoulli <- function(...) {
    generated_quantities_lines_default(...)
}

generated_quantities_lines_poisson <- function(...) {
    generated_quantities_lines_default(...)
}

generated_quantities_lines_negbin <- function(...) {
    generated_quantities_lines_default(...)
}

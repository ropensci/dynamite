# NOTE: Added this wrapper, so that *_lines_* functions can have different arguments (including parameters and model)
lines_wrap <- function(prefix, formula, args) {
    lines_fun <- paste0(prefix, "_lines_", formula$family)
    # has_args <- names(args) %in% names(formals(get(lines_fun)))
    # args[!has_args] <- NULL
    do.call(lines_fun, args)
}
vectorizable_prior <- function(x) length(x) == 1 && !grepl("\\(", x)

# For data block
data_lines_default <- function(i, idt, has_fixed, has_varying, data, ...) {

    write_beta_fixed <- has_fixed && length(data[[paste0("beta_fixed_prior_distr_", i)]]) == 1
    write_beta_varying <- has_varying && length(data[[paste0("beta_varying_prior_distr_", i)]]) == 1
    write_tau <- has_varying && length(data[[paste0("tau_prior_distr_", i)]]) == 1

    paste_rows(onlyif(write_beta_fixed, c(idt(1), "int beta_fixed_prior_npars_", i, ";")),
               onlyif(write_beta_fixed, c(idt(1), "matrix[K_fixed_", i, ", beta_fixed_prior_npars_", i, "] beta_fixed_prior_pars_", i, ";")),
               onlyif(write_beta_varying, c(idt(1), "int beta_varying_prior_npars_", i, ";")),
               onlyif(write_beta_varying, c(idt(1), "matrix[K_varying_", i, ", beta_varying_prior_npars_", i, "] beta_varying_prior_pars_", i, ";")),
               onlyif(write_tau, c(idt(1), "int tau_prior_npars_", i, ";")),
               onlyif(write_tau, c(idt(1), "matrix[K_varying_", i, ", tau_prior_npars_", i, "] tau_prior_pars_", i, ";")))

}

data_lines_categorical <- function(i, idt, has_fixed, has_varying, data, ...) {

    write_beta_fixed <- has_fixed && length(data[[paste0("beta_fixed_prior_distr_", i)]]) == 1
    write_beta_varying <- has_varying && length(data[[paste0("beta_varying_prior_distr_", i)]]) == 1
    write_tau <- has_varying && length(data[[paste0("tau_prior_distr_", i)]]) == 1

    paste_rows(c(idt(1), "int<lower=1> ", i, "[T, N];"),
               c(idt(1), "int<lower=0> S_", i, ";"),
               onlyif(write_beta_fixed, c(idt(1), "int beta_fixed_prior_npars_", i, ";")),
               onlyif(write_beta_fixed, c(idt(1), "matrix[K_fixed_", i, " * (S_", i, " - 1), beta_fixed_prior_npars_", i, "] beta_fixed_prior_pars_", i, ";")),
               onlyif(write_beta_varying, c(idt(1), "int beta_varying_prior_npars_", i, ";")),
               onlyif(write_beta_varying, c(idt(1), "matrix[K_varying_", i, " * (S_", i, " - 1), beta_varying_prior_npars_", i, "] beta_varying_prior_pars_", i, ";")),
               onlyif(write_tau, c(idt(1), "int tau_prior_npars_", i, ";")),
               onlyif(write_tau, c(idt(1), "matrix[K_varying_", i, ", tau_prior_npars_", i, "] tau_prior_pars_", i, ";")))
}

data_lines_gaussian <- function(i, idt, has_fixed, has_varying, data, ...) {
    paste_rows(c(idt(1), "real ", i, "[T, N];"),
               data_lines_default(i, idt, has_fixed, has_varying, data, ...))
}

data_lines_binomial <- function(i, idt, has_fixed, has_varying, data, ...) {
    paste_rows(c(idt(1), "int<lower=0> ", i, "[T, N];"),
               c(idt(1), "int<lower=1> trials_", i, "[T, N];"),
               data_lines_default(i, idt, has_fixed, has_varying, data, ...))
}

data_lines_bernoulli <- function(i, idt, has_fixed, has_varying, data, ...) {
    paste_rows(c(idt(1), "int<lower=0,upper=1> ", i, "[T, N];"),
               data_lines_default(i, idt, has_fixed, has_varying, data, ...))
}

data_lines_poisson <- function(i, idt, has_offset, has_fixed, has_varying, data, ...) {
    paste_rows(c(idt(1), "int<lower=0> ", i, "[T, N];"),
               onlyif(has_offset, c(idt(1), "real offset_", i, "[T, N];")),
               data_lines_default(i, idt, has_fixed, has_varying, data, ...))
}

data_lines_negbin <- function(i, idt, has_fixed, has_varying, data, ...) {
    paste_rows(c(idt(1), "int<lower=0> ", i, "[T, N];"),
               data_lines_default(i, idt, has_fixed, has_varying, data, ...))
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
               c(idt(1), "real<lower=0> phi_", i, ";"))
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

transformed_parameters_lines_negbin <- function(...) {
    transformed_parameters_lines_bernoulli(...)
}

# For model block
model_lines_default <- function(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...) {
    mtext_fixed <- ""
    mtext_varying <- ""
    mtext_tau <- ""
    if (has_fixed) {
        #d <- priors |> dplyr::filter(type == "beta_fixed")
        # if (d$vectorized[1]) {
        #     distribution <- sub("(.*", "", d$prior[1])
        #     mtext_fixed <- c(idt(1), distribution, "~ (", "beta_fixed_")
        d <- data[[paste0("beta_fixed_prior_distr_", i)]]
        if (vectorizable_prior(d)) {
            np <- data[[paste0("beta_fixed_prior_npars_", i)]]
            mtext_fixed <- c(idt(1), "beta_fixed_", i, " ~ ", d, "(", paste0("beta_fixed_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
        } else {
            mtext_fixed <- c(idt(1), "beta_fixed_", i, "[k] ~ ", d)
        }
    }
    if (has_varying) {
        if (noncentered) {
            mtext_varying <- paste_rows(
                c(idt(1), "for (k in 1:K_varying_", i, ") {"),
                c(idt(2),      "a_raw_", i, "[k] ~ std_normal();"),
                c(idt(1), "}")
            )
        } else {
            d <- data[[paste0("beta_varying_prior_distr_", i)]]
            if (vectorizable_prior(d)) {
                np <- data[[paste0("beta_varying_prior_npars_", i)]]
                mtext_varying <- c(idt(1), "a_", i, "[, 1] ~ ", d, "(", paste0("beta_varying_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
            } else {
                K <- data[[paste0("K_varying_", i)]]
                mtext_varying <- do.call(paste_rows, as.list(paste0(idt(1), "a_", i, "[", 1:K, "] ~ ", d, ";")))
            }
            mtext_varying <- paste_rows(mtext_varying,
                c(idt(1), "for (k in 1:K_varying_", i, ") {"),
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
        d <- data[[paste0("tau_prior_distr_", i)]]
        if (vectorizable_prior(d)) {
            np <- data[[paste0("tau_prior_npars_", i)]]
            mtext_tau <- c(idt(1), "tau_", i, " ~ ", d, "(", paste0("tau_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
        } else {
            K <- data[[paste0("K_varying_", i)]]
            mtext_tau <- do.call(paste_rows, as.list(paste0(idt(1), "tau_", i, "[", 1:K, "] ~ ", d, ";")))
        }
    }
    paste_rows(mtext_fixed, mtext_varying, mtext_tau)
}

model_lines_categorical <- function(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...) {
    mtext_fixed <- ""
    mtext_varying <- ""
    mtext_tau <- ""
    if (has_fixed) {
        d <- data[[paste0("beta_fixed_prior_distr_", i)]]
        if (vectorizable_prior(d)) {
            np <- data[[paste0("beta_fixed_prior_npars_", i)]]
            mtext_fixed <- c(idt(1), "to_vector(beta_fixed_", i, ") ~ ", d, "(", paste0("beta_fixed_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
        } else {
            K <- data[[paste0("K_fixed_", i)]]
            S <- data[[paste0("S_", i)]]
            k <- rep(1:K, S - 1)
            s <- rep(1:(S - 1), each = K)
            mtext_fixed <- do.call(paste_rows, as.list(paste0(idt(1), "beta_fixed_", i, "[", k, ",", s, "] ~ ", d, ";")))
        }

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
            d <- data[[paste0("beta_varying_prior_distr_", i)]]
            if (vectorizable_prior(d)) {
                np <- data[[paste0("beta_varying_prior_npars_", i)]]
                mtext_varying <- c(idt(1), "to_vector(beta_", i, "[1, L_varying_", i, ", ]) ~ ", d, "(", paste0("beta_varying_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
            } else {
                K <- data[[paste0("K_fixed_", i)]]
                S <- data[[paste0("S_", i)]]
                k <- rep(1:K, S - 1)
                s <- rep(1:(S - 1), each = K)
                mtext_varying <- do.call(paste_rows, as.list(paste0(idt(1), "a_", i, "[", k, ",", s, "] ~ ", d, ";")))
            }
            mtext_varying <- paste_rows(mtext_varying,
                c(idt(1), "for (s in 1:(S_", i, " - 1)) {"),
                c(idt(2),     "for (k in 1:K_varying_", i, ") {"),
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
        d <- data[[paste0("tau_prior_distr_", i)]]
        if (vectorizable_prior(d)) {
            np <- data[[paste0("tau_prior_npars_", i)]]
            mtext_tau <- c(idt(1), "tau_", i, " ~ ", d, "(", paste0("tau_prior_pars_", i, "[, ", 1:np, "]", collapse = ", "), ");")
        } else {
            K <- data[[paste0("K_varying_", i)]]
            mtext_tau <- do.call(paste_rows, as.list(paste0(idt(1), "tau_", i, "[", 1:K, "] ~ ", d, ";")))
        }
    }
    likelihood_term <-
        c(idt(2), i, "[t] ~ categorical_logit_glm(X[t][,J_", i, "], zeros_S_", i, ", append_col(beta_", i, "[t], zeros_K_", i, "));")

    paste_rows(mtext_fixed, mtext_varying, mtext_tau, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_gaussian <- function(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...) {

    mtext <- model_lines_default(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...)
    d <- data[[paste0("sigma_prior_distr_", i)]]
    sigma_term <- c(idt(1), "sigma_", i, " ~ ", d, ";")

    fixed_term <- onlyif(has_fixed, paste0("X[t][,J_fixed_", i, "] * beta_fixed_", i))
    varying_term <- onlyif(has_varying, paste0("X[t][,J_varying_", i, "] * beta_varying_", i, "[t]"))
    plus <- onlyif(has_fixed && has_varying, " + ")
    likelihood_term <-  c(idt(2), i, "[t] ~ normal(", fixed_term, plus, varying_term, ", sigma_", i, ");")
    paste_rows(mtext, sigma_term, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_binomial <- function(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...) {
    mtext <- model_lines_default(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...)
    fixed_term <- onlyif(has_fixed, paste0("X[t][,J_fixed_", i, "] * beta_fixed_", i))
    varying_term <- onlyif(has_varying, paste0("X[t][,J_varying_", i, "] * beta_varying_", i, "[t]"))
    plus <- onlyif(has_fixed && has_varying, " + ")
    likelihood_term <- c(idt(2), i, "[t] ~ binomial_logit(trials_", i, "[t], ", fixed_term, plus, varying_term, ");")
    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_bernoulli <- function(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...) {
    mtext <- model_lines_default(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...)
    likelihood_term <- c(idt(2), i, "[t] ~ bernoulli_logit_glm(X[t][,J_", i, "], 0, beta_", i, "[t]);")
    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_poisson <-  function(i, idt, has_offset, shrinkage, noncentered, has_fixed, has_varying, data, ...) {
    mtext <- model_lines_default(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...)
    offset_term <- ifelse_(has_offset, paste0("to_vector(offset_", i, "[t])"), "0")
    likelihood_term <- c(idt(2), i, "[t] ~ poisson_log_glm(X[t][,J_", i, "], ", offset_term, ", beta_", i, "[t]);")
    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_negbin <- function(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...) {
    mtext <- model_lines_default(i, idt, shrinkage, noncentered, has_fixed, has_varying, data, ...)
    d <- data[[paste0("phi_prior_distr_", i)]]
    phi_term <- c(idt(1), "phi_", i, " ~ ", d, ";")
    likelihood_term <- c(idt(2), i, "[t] ~ neg_binomial_2_log_glm(X[t][,J_", i, "], 0, beta_", i, "[t], ", "phi_", i, ");")
    paste_rows(mtext, phi_term, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
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

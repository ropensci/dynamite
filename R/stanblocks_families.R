
# For data block
data_lines_default <- function(i, idt) {
    mtext <- c(idt(1), "vector[K_", i, "] a_prior_mean_", i, ";")
    paste_rows(mtext, c(idt(1), "vector[K_", i, "] a_prior_sd_", i, ";"))
}

data_lines_categorical <- function(i, idt) {
    mtext <- c(idt(1), "int<lower=1> response_", i, "[T, N];")
    mtext <- paste_rows(mtext, c(idt(1), "int<lower=0> S_", i, ";"))
    mtext <- paste_rows(mtext, c(idt(1), "matrix[K_", i, ", S_", i, " - 1] a_prior_mean_", i, ";"))
    paste_rows(mtext, c(idt(1), "matrix[K_", i, ", S_", i, " - 1] a_prior_sd_", i, ";"))
}

data_lines_gaussian <- function(i, idt) {
    mtext <- c(idt(1), "real response_", i, "[T, N];")
    mtext <- paste_rows(mtext, c(idt(1), "real<lower=0> sigma_scale_", i, ";"))
    paste_rows(mtext, data_lines_default(i, idt))
}

data_lines_binomial <- function(i, idt) {
    mtext <- c(idt(1), "int<lower=0> response_", i, "[T, N];")
    mtext <- paste_rows(mtext, c(idt(1), "int<lower=1>[N, T] trials_", i, ";"))
    paste_rows(mtext, data_lines_default(i, idt))
}

data_lines_bernoulli <- function(i, idt) {
    mtext <- c(idt(1), "int<lower=0,upper=1> response_", i, "[T, N];")
    paste_rows(mtext, data_lines_default(i, idt))
}

data_lines_poisson <- function(i, idt) {
    mtext <- c(idt(1), "int<lower=0> response_", i, "[T, N];")
    paste_rows(mtext, data_lines_default(i, idt))
}

data_lines_negbin <- function(i, idt) {
    mtext <- c(idt(1), "int<lower=0> response_", i, "[T, N];")
    mtext <- c(idt(1), "real<lower=0> phi_scale_", i)
    paste_rows(mtext, data_lines_default(i, idt))
}

# For parameters block
parameters_lines_default <- function(i, lb, idt) {
    a_term <- c(idt(1), "row_vector[D] a_", i, "[K_", i, "];")
    tau_term <- c(idt(1), "vector<lower=", lb, ">[K_", i, "] tau_", i, ";")
    paste_rows(a_term, tau_term)
}

parameters_lines_categorical <- function(i, lb, idt) {
    a_term <- c(idt(1), "row_vector[D] a_", i, "[S_", i, "- 1, K_", i, "];")
    tau_term <- c(idt(1), "vector<lower=", lb, ">[K_", i, "] tau_", i, ";")
    paste_rows(a_term, tau_term)
}

parameters_lines_gaussian <- function(i, lb, idt) {
    sigma_term <- c(idt(1), "real<lower=0> sigma_", i, ";")
    paste_rows(parameters_lines_default(i, lb, idt), sigma_term)
}

parameters_lines_binomial <- function(i, lb, idt) {
    parameters_lines_default(i, lb, idt)
}

parameters_lines_bernoulli <- function(i, lb, idt) {
    parameters_lines_default(i, lb, idt)
}

parameters_lines_poisson <- function(i, lb, idt) {
    parameters_lines_default(i, lb, idt)
}

parameters_lines_negbin <- function(i, lb, idt) {
    phi_term <- c(idt(1), "real<lower=0 phi_", i)
    paste_rows(parameters_lines_default(i, lb, idt), phi_term)
}

# For model block
model_lines_default <- function(i, shrinkage, noncentered, idt) {
    if (noncentered) {
        mtext <- paste_rows(
            c(idt(1), "for (k in 1:K_", i, ") {"),
            c(idt(2), "a_raw_", i, "[k] ~ std_normal();"),
            c(idt(1), "}")
        )
    } else {
        # Prior for the first a (beta) is always normal given the RW prior
        mtext <- paste_rows(
            c(idt(1), "for (k in 1:K_", i, ") {"),
            c(idt(2), "a_", i, "[k, 1] ~ normal(a_prior_mean_", i, "[k], a_prior_sd_", i, "[k]);"),
            c(idt(2), "for(i in 2:D) {"),
            if (shrinkage) {
                c(idt(3), "a_", i, "[k, i] ~ normal(a_", i,"[k, i - 1], lambda[i - 1] * tau_", i, "[k]);")
            } else {
                c(idt(3), "a_", i, "[k, i] ~ normal(a_", i,"[k, i - 1], tau_", i, "[k]);")
            },
            c(idt(2), "}"),
            c(idt(1), "}")
        )
    }
    mtext
}

model_lines_categorical <- function(i, shrinkage, noncentered, idt) {
    if(noncentered) {
        # fixed given noncentered parameterisation
        a_term <- paste_rows(
            c(idt(1), "for (s in 1:(S_", i, " - 1)) {"),
            c(idt(2), "for (k in 1:K_", i, ") {"),
            c(idt(3), "a_raw_", i, "[s, k] ~ std_normal();"),
            c(idt(2), "}"),
            c(idt(1), "}"))
    } else {
        # Prior for the first a (beta) is always normal given the RW prior
        a_term <- paste_rows(
            c(idt(1), "for (s in 1:(S_", i, " - 1)) {"),
            c(idt(2), "for (k in 1:K_", i, ") {"),
            c(idt(3), "a_", i, "[s, k, 1] ~ normal(a_prior_mean_", i, "[k, s], a_prior_sd_", i, "[k, s]);"),
            c(idt(3), "for(i in 2:D) {"),
            if (shrinkage) {
                c(idt(4), "a_", i, "[s, k, i] ~ normal(a_", i,"[s, k, i - 1], lambda[i - 1] * tau_", i, "[k]);")
            } else {
                c(idt(4), "a_", i, "[s, k, i] ~ normal(a_", i,"[s, k, i - 1], tau_", i, "[k]);")
            },
            c(idt(3), "}"),
            c(idt(2), "}"),
            c(idt(1), "}")
        )
    }
    likelihood_term <-
        c(idt(2), "response_", i, "[t] ~ categorical_logit_glm(X[t][,J_", i, "], zeros_S_", i, ", beta_", i, "[t]);")

    paste_rows(a_term, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_gaussian <- function(i, shrinkage, noncentered, idt) {
    mtext <- model_lines_default(i, shrinkage, noncentered, idt)
    # TODO user-defined prior for sigma
    mtext <- paste_rows(mtext, c(idt(1), "sigma_", i, " ~ exponential(sigma_scale_", i, ");"))

    likelihood_term <-
        c(idt(2), "response_", i, "[t] ~ normal(X[t][,J_", i, "] * beta_", i, "[t], ", "sigma_", i, ");")

    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_binomial <- function(i, shrinkage, noncentered, idt) {
    mtext <- model_lines_default(i, shrinkage, noncentered, idt)
    likelihood_term <-
        c(idt(2), "response_", i, "[t] ~ binomial_logit(X[t][,J_", i, "] * beta_", i, "[t], trials_, ", i, "[t]);")

    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_bernoulli <- function(i, shrinkage, noncentered, idt) {
    mtext <- model_lines_default(i, shrinkage, noncentered, idt)
    likelihood_term <-
        c(idt(2), "response_", i, "[t] ~ bernoulli_logit_glm(X[t][,J_", i, "], 0, beta_", i, "[t]);")

    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_poisson <- function(i, shrinkage, noncentered, idt) {
    mtext <- model_lines_default(i, shrinkage, noncentered, idt)
    likelihood_term <-
        c(idt(2), "response_", i, "[t] ~ poisson_log_glm(X[t][,J_", i, "], 0, beta_", i, "[t]);")

    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

model_lines_negbin <- function(i, shrinkage, noncentered, idt) {
    mtext <- model_lines_default(i, shrinkage, noncentered, idt)
    # TODO user-defined prior for phi
    mtext <- paste_rows(mtext, c(idt(1), "phi_", i, " ~ exponential(phi_scale_", i, ");"))

    likelihood_term <-
        c(idt(2), "response_", i, "[t] ~ neg_binomial_2_log_glm(X[t][,J_", i, "], 0, beta_", i, "[t], ", "phi_", i, ");")

    paste_rows(mtext, c(idt(1), "for (t in 1:T) {"), likelihood_term, c(idt(1) ,"}"))
}

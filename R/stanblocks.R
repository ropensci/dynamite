#' Create code blocks for the Stan model
#'
#' @param becauseformula defining model components
#' @export
create_blocks <- function(formula, ...) {
    UseMethod("create_blocks")
}

#' @export
create_blocks.default <- function(formula, indent = 2L, ...) {

    # only need this for the hidden state case
    #functions <- paste("functions {", .loglik_stan, "}", sep = "\n")
    idt <- indenter_(indent)
    functions <- create_functions(formula, idt, ...)
    data <- create_data(formula, idt, ...)
    transformed_data <- create_transformed_data(formula, idt, ...)
    parameters <- create_parameters(formula, idt, ...)
    transformed_parameters <- create_transformed_parameters(formula, idt, ...)
    model <- create_model(formula, idt, ...)
    generated_quantities <- create_generated_quantities(formula, idt, ...)
    # combine above text blocks
    model_code <- paste_rows(functions, data, transformed_data, parameters,
        transformed_parameters, model, generated_quantities)
    model_code
}
#'
#' @export
create_functions <- function(formula, idt, ...) {
    NULL
}
#'
#' @export
create_data <- function(formula, idt, ...) {

    # Maybe not the best way...
    mtext <- paste_rows(
        c(idt(1), "int<lower=1> T; // number of time points"),
        c(idt(1), "int<lower=1> N; // number of individuals"),
        c(idt(1), "int<lower=1> C; // number of channels/response variables"),
        c(idt(1), "int<lower=1> K; // total number of covariates across all channels"),
        c(idt(1), "matrix[N, K] X[T]; // all covariates as an array of N x K matrices"),
        c(idt(1), "int<lower=0> D; // number of B-splines"),
        c(idt(1), "matrix[D, T] Bs; // B-spline basis matrix")
    )

    # loop over channels
    for (i in seq_along(formula)) {
        # Number of covariates in channel i
        mtext <- paste_rows(mtext, c(idt(1), "int<lower=1> K_", i, ";"))
        # index vector of covariates related to channel i
        mtext <- paste_rows(mtext, c(idt(1), "int J_", i, "[K_", i, "];"))
        mtext <- paste_rows(mtext,
            do.call(paste0("data_lines_", formula[[i]]$family), list(i = i, idt = idt)))
    }
    paste_rows("data {", mtext, "}")
}

#'
#' @export
create_transformed_data <- function(formula, idt, ...) {
    transformed_data <- character(0)
    for (i in seq_along(formula)) {
        if (is_categorical(formula[[i]]$family)) {
            zeros_term <- paste_rows(
                c(idt(1), "vector[K_", i, "] zeros_K_", i, " = rep_vector(0, K_", i, ");"),
                c(idt(1), "vector[S_", i, "] zeros_S_", i, " = rep_vector(0, S_", i, ");")
            )
            c(transformed_data) <- zeros_term
        }
    }
    paste_rows("transformed data {", collapse_rows(transformed_data), "}")
}

#'
#' @export
create_parameters <- function(formula, idt, ...) {
    lb <- attr(formula, "splines")$lb_tau
    # TODO channel-wise
    if (attr(formula, "splines")$noncentered) {
        stop("Noncentered parameterisation is currently not supported.")
        mtext <- paste_rows(
            c(idt(1), "// Spline parameters"),
            c(idt(1), "// use noncentered parameterisation for sampling efficiency")
        )
        for (i in seq_along(formula)) {

            if (is_categorical(formula[[i]]$family)) {
                a_raw_term <- c(idt(1), "row_vector[D] a_raw_", i, "[S_", i, "- 1, K_", i, "];")
                # Do we want separate tau for K coefficients or K * (S - 1) coefficients?
                tau_term <- c(idt(1), "vector<lower=", lb, ">[K_", i, "] tau_", i, ";")
            } else {
                a_raw_term <- c(idt(1), "row_vector[D] a_raw_", i, "[K_", i, "];")
                tau_term <- c(idt(1), "vector<lower=", lb, ">[K_", i, "] tau_", i, ";")
            }
            mtext <- paste_rows(mtext, a_raw_term, tau_term)
            if (is_gaussian(formula[[i]]$family)) {
                sigma_term <- c(idt(1), "real<lower=0> sigma_", i, ";")
                mtext <- paste_rows(mtext, sigma_term)
            }
            if (is_gamma(formula[[i]]$family)) {
                # TODO: shape parameter for the gamma distribution, not a priority
                stop("Gamma distribution is not yet supported")
            }

        }
    } else {
        mtext <- paste0(idt(1), "// Spline parameters")
        for (i in seq_along(formula)) {
            mtext <- paste_rows(mtext,
                do.call(paste0("parameters_lines_", formula[[i]]$family), list(i = i, lb = lb, idt = idt)))
        }
    }
    if (attr(formula, "splines")$shrinkage) {
        mtext <- paste_rows(mtext, c(idt(1), "vector<lower=0>[D - 1] lambda; // shrinkage parameter"))
    }
    paste_rows("parameters {", mtext, "}")
}

#'
#' @export
create_transformed_parameters <- function(formula, idt, ...) {

    # define variables
    beta_terms <- character(0)
    a_terms <- character(0)
    for (i in seq_along(formula)) {
        if (is_categorical(formula[[i]]$family)) {
            c(beta_terms) <- paste0(idt(1), "matrix[K_", i, ", S_", i, "] beta_", i, "[T];")
            if (attr(formula, "splines")$noncentered) { # TODO: channel-wise?
                c(a_terms) <- paste0(idt(1), "row_vector[D] a_", i, "[S_", i, " - 1, K_", i, "];")
            }
        } else {
            c(beta_terms) <- paste0(idt(1), "vector[K_", i, "] beta_", i, "[T];")
            if (attr(formula, "splines")$noncentered) { # TODO: channel-wise?
                c(a_terms) <- paste0(idt(1), "row_vector[D] a_", i, "[K_", i, "];")
            }
        }
    }
    beta_terms <- collapse_rows(beta_terms)
    a_terms <- collapse_rows(a_terms)
    mtext <- paste_rows(beta_terms, a_terms)
    for (i in seq_along(formula)) {

        if (attr(formula, "splines")$noncentered) {
            stop("Noncentered parameterisation is currently not supported.")
            # TODO: could separate construction of a and beta so less repetition in the codes
            # TODO: check whether lambda is used
            # TODO use prior definitions
            # i.e. a_i[s, k, 1] = a_prior_mean_i[s, k] + a_prior_sd_i[s, k] * a_raw_i[s, k, 1]
            if (is_categorical(formula[[i]]$family)) {
                spline_term <- paste_rows(
                    # TODO define zeros_K_i in transformed data?
                    c(idt(1),  "for (s in 1:(S_", i, " - 1)) {"),
                    c(idt(2), "for (k in 1:K_", i, ") {"),
                    c(idt(3), "a_", i, "[s, k, 1] = a_raw_", i, "[s, k, 1];"),
                    c(idt(3), "for (i in 2:D) {"),
                    c(idt(4), "a_", i, "[s, k, i] = a_", i, "[s, k, i-1] + a_raw_", i, "[s, k, i] * tau_", i, "[k] * lambda[i - 1];"),
                    c(idt(3), "}"),
                    c(idt(3), "for (t in 1:T) {"),
                    c(idt(4), "beta_", i, "[t, k, s] = a_", i, "[s, k] * Bs[, t];"),
                    c(idt(3), "}"),
                    c(idt(2), "}"),
                    c(idt(1), "}"),
                    c(idt(1), "for (t in 1:T) {"),
                    c(idt(2), "beta_", i, "[t, , S_", i, "] = zeros_K_", i, ";"),
                    c(idt(2), "}")
                )
            } else {
                spline_term <- paste_rows(
                    c(idt(1), "for (k in 1:K_", i, ") {"),
                    c(idt(2), "a_", i, "[k, 1] = a_raw_", i, "[k, 1];"),
                    c(idt(2), "for (i in 2:D) {"),
                    c(idt(3), "a_", i, "[k, i] = a_", i, "[k, i-1] + a_raw_", i, "[k, i] * tau_", i, "[k] * lambda[i - 1];"),
                    c(idt(2), "}"),
                    c(idt(2), "for (t in 1:T) {"),
                    c(idt(3), ", beta_", i, "[t, k] = a_", i, "[k] * Bs[, t];"),
                    c(idt(2), "}"),
                    c(idt(1), "}")
                )
            }
        } else {
            if (is_categorical(formula[[i]]$family)) {
                spline_term <- paste_rows(
                    c(idt(1), "for (s in 1:(S_", i, " - 1)) {"),
                    c(idt(2), "for (k in 1:K_", i, ") {"),
                    c(idt(3), "for (t in 1:T) {"),
                    c(idt(4), "beta_", i, "[t, k, s] = a_", i, "[s, k] * Bs[, t];"),
                    c(idt(3), "}"),
                    c(idt(2), "}"),
                    c(idt(1), "}"),
                    c(idt(1), "for (t in 1:T) {"),
                    c(idt(2), "beta_", i, "[t, , S_", i, "] = zeros_K_", i, ";"),
                    c(idt(1), "}")
                )
            } else {
                spline_term <- paste_rows(
                    c(idt(1), "for (k in 1:K_", i, ") {"),
                    c(idt(2), "for (t in 1:T) {"),
                    c(idt(3), "beta_", i, "[t, k] = a_", i, "[k] * Bs[, t];"),
                    c(idt(2), "}"),
                    c(idt(1), "}")
                )
            }
        }
        mtext <- paste_rows(mtext, spline_term)
    }
    paste_rows("transformed parameters {", mtext, "}")
}

#'
#' @export
create_model <- function(formula, idt, ...) {
    # TODO: Without global shrinkage prior it probably makes sense to use user-defined prior for tau
    # With lambda&tau, need more testing if this is fine or do we need to support other forms
    # e.g. as in https://arxiv.org/abs/1611.01310 and https://www.mdpi.com/2225-1146/8/2/20
    priors <- character(0)
    if (attr(formula, "splines")$shrinkage) {
        c(priors) <- paste0(idt(1), "lambda ~ std_normal();  // prior for shrinkage terms")
    }
    c(priors) <- paste0(idt(1), "tau_", 1:length(formula), " ~ normal(0, 1);")
    mtext <- character(0)
    for (i in seq_along(formula)) {
            mtext <- paste_rows(mtext,
                do.call(paste0("model_lines_", formula[[i]]$family),
                list(i = i, shrinkage = attr(formula, "splines")$shrinkage,
                    noncentered = attr(formula, "splines")$noncentered, idt = idt)))
    }
    mtext <- paste_rows(collapse_rows(priors), mtext)
    mtext <- paste_rows("model {", mtext, "}")
    mtext
}

#'
#' @export
create_generated_quantities <- function(formula, ...) {
    NULL
}

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
create_data <- function(formula, idt, resp, helpers, data, ...) {

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
    channels <- character(length(formula))
    for (i in seq_along(formula)) {
        y <- resp[i]
        h <- helpers[[i]]
        # Number of covariates in channel i
        line_args <- c(list(i = y, idt = idt, data), h)
        channels[i] <- paste_rows(
            # number of covariates (all, fixed, varying)
            c(idt(1), "int<lower=1> K_", y, ";"),
            onlyif(h$has_fixed, c(idt(1), "int<lower=1> K_fixed_", y, ";")),
            onlyif(h$has_varying, c(idt(1), "int<lower=1> K_varying_", y, ";")),
            # index vectors of covariates (all, fixed, varying) related to channel i
            c(idt(1), "int J_", y, "[K_", y, "];"),
            onlyif(h$has_fixed, c(idt(1), "int J_fixed_", y, "[K_fixed_", y, "];")),
            onlyif(h$has_varying, c(idt(1), "int J_varying_", y, "[K_varying_", y, "];")),
            # index vectors of fixed and varying betas
            onlyif(h$has_fixed, c(idt(1), "int L_fixed_", y, "[K_fixed_", y, "];")),
            onlyif(h$has_varying, c(idt(1), "int L_varying_", y, "[K_varying_", y, "];")),
            lines_wrap("data", formula[[i]], line_args)
        )
    }
    paste_rows("data {", mtext, collapse_rows(channels), "}")
}

#'
#' @export
create_transformed_data <- function(formula, idt, resp, ...) {
    transformed_data <- character(length(formula))
    for (i in seq_along(formula)) {
        if (is_categorical(formula[[i]]$family)) {
            y <- resp[i]
            transformed_data[i] <- paste_rows(
                c(idt(1), "vector[K_", y, "] zeros_K_", y, " = rep_vector(0, K_", y, ");"),
                c(idt(1), "vector[S_", y, "] zeros_S_", y, " = rep_vector(0, S_", y, ");")
            )
        }
    }
    paste_rows("transformed data {", collapse_rows(transformed_data), "}")
}

#'
#' @export
create_parameters <- function(formula, idt, resp, helpers, ...) {
    lb <- character(0)
    splinetext <- ""
    if (!is.null(spline_defs <- attr(formula, "splines"))) {
        if (spline_defs$noncentered) {
            stop_("Noncentered parameterisation is currently not supported.")
        }
        lb <- attr(formula, "splines")$lb_tau
        splinetext <- paste0(idt(1), "// Spline parameters")
        if (spline_defs$shrinkage) {
            splinetext <- paste_rows(splinetext, c(idt(1), "vector<lower=0>[D - 1] lambda; // shrinkage parameter"))
        }
        # TODO handle centered case where spline is not defined but user inserts varying(.) terms
    }
    ## TODO channel-wise
    #if (attr(formula, "splines")$noncentered) {
    #    stop("Noncentered parameterisation is currently not supported.")
    #    mtext <- paste_rows(
    #        c(idt(1), "// Spline parameters"),
    #        c(idt(1), "// use noncentered parameterisation for sampling efficiency")
    #    )
    #    for (i in seq_along(formula)) {
    #        if (is_categorical(formula[[i]]$family)) {
    #            a_raw_term <- c(idt(1), "row_vector[D] a_raw_", i, "[S_", i, "- 1, K_", i, "];")
    #            # Do we want separate tau for K coefficients or K * (S - 1) coefficients?
    #            tau_term <- c(idt(1), "vector<lower=", lb, ">[K_", i, "] tau_", i, ";")
    #        } else {
    #            a_raw_term <- c(idt(1), "row_vector[D] a_raw_", i, "[K_", i, "];")
    #            tau_term <- c(idt(1), "vector<lower=", lb, ">[K_", i, "] tau_", i, ";")
    #        }
    #        mtext <- paste_rows(mtext, a_raw_term, tau_term)
    #        if (is_gaussian(formula[[i]]$family)) {
    #            sigma_term <- c(idt(1), "real<lower=0> sigma_", i, ";")
    #            mtext <- paste_rows(mtext, sigma_term)
    #        }
    #        if (is_gamma(formula[[i]]$family)) {
    #            # TODO: shape parameter for the gamma distribution, not a priority
    #            stop("Gamma distribution is not yet supported")
    #        }
    #    }
    #} else {
    pars <- character(length(formula))
    for (i in seq_along(formula)) {
        line_args <- c(list(i = resp[i], lb = lb, idt = idt), helpers[[i]])
        pars[i] <- lines_wrap("parameters", formula[[i]], line_args)
    }
    paste_rows("parameters {", splinetext, collapse_rows(pars), "}")
}

#'
#' @export
create_transformed_parameters <- function(formula, idt, resp, helpers, ...) {
    #if (spline_defs$noncentered) {
    #stop("Noncentered parameterisation is currently not supported.")
    # TODO: could separate construction of a and beta so less repetition in the codes
    # TODO: check whether lambda is used
    # TODO use prior definitions
    # i.e. a_i[s, k, 1] = a_prior_mean_i[s, k] + a_prior_sd_i[s, k] * a_raw_i[s, k, 1]
    # if (is_categorical(formula[[i]]$family)) {
    #     c(varying_terms) <- paste_rows(
    #         c(idt(1), "for (s in 1:(S_", y, " - 1)) {"),
    #         c(idt(2),     "for (k in 1:K_varying_", y, ") {"),
    #         c(idt(3),         "a_", y, "[s, k, 1] = a_raw_", y, "[s, k, 1];"),
    #         c(idt(3),         "for (i in 2:D) {"),
    #         c(idt(4),             "a_", y, "[s, k, i] = a_", y, "[s, k, i-1] + a_raw_", y, "[s, k, i] * tau_", y, "[k] * lambda[i - 1];"),
    #         c(idt(3),         "}"),
    #         c(idt(3),         "for (t in 1:T) {"),
    #         c(idt(4),             "beta_", y, "[t, J_varying_", y, "[k], s] = a_", y, "[s, k] * Bs[, t];"),
    #         c(idt(3),         "}"),
    #         c(idt(2),     "}"),
    #         c(idt(1), "}")
    #     )
    # } else {
    #     c(varying_terms) <- paste_rows(
    #         c(idt(1), "for (k in 1:K_", y, ") {"),
    #         c(idt(2),     "a_", y, "[k, 1] = a_raw_", y, "[k, 1];"),
    #         c(idt(2),     "for (i in 2:D) {"),
    #         c(idt(3),         "a_", y, "[k, i] = a_", y, "[k, i - 1] + a_raw_", y, "[k, i] * tau_", y, "[k] * lambda[i - 1];"),
    #         c(idt(2),     "}"),
    #         c(idt(2),     "for (t in 1:T) {"),
    #         c(idt(3),          "beta_", y, "[t, k] = a_", y, "[k] * Bs[, t];"),
    #         c(idt(2),     "}"),
    #         c(idt(1), "}")
    #     )
    # }
    spline_defs <- attr(formula, "splines")
    transpars <- character(length(formula))
    for (i in seq_along(formula)) {
        line_args <- c(list(i = resp[i],
                            idt = idt,
                            noncentered = spline_defs$noncetered),
                       helpers[[i]])
        transpars[i] <- lines_wrap("transformed_parameters", formula[[i]], line_args)

    }
    paste_rows("transformed parameters {", collapse_rows(transpars), "}")
}

#'
#' @export
create_model <- function(formula, idt, resp, helpers, priors, data, ...) {
    # TODO: Without global shrinkage prior it probably makes sense to use user-defined prior for tau
    # With lambda&tau, need more testing if this is fine or do we need to support other forms
    # e.g. as in https://arxiv.org/abs/1611.01310 and https://www.mdpi.com/2225-1146/8/2/20
    #priors <- character(0)

    # if (!is.null(spline_defs <- attr(formula, "splines"))) {
    #     if (spline_defs$shrinkage) {
    #         priors <- paste0(idt(1), "lambda ~ std_normal();  // prior for shrinkage terms")
    #     }
    # }
    spline_defs <- attr(formula, "splines")
    mod <- character(length(formula))
    for (i in seq_along(formula)) {
        line_args <- c(list(i = resp[i],
                            idt = idt,
                            shrinkage = spline_defs$shrinkage,
                            noncentered = spline_defs$noncentered,
                            data = data),
                       helpers[[i]],
                       priors = priors[[i]])
        mod[i] <- lines_wrap("model", formula[[i]], line_args)
    }
    paste_rows("model {", collapse_rows(mod), "}")
}

#'
#' @export
create_generated_quantities <- function(formula, idt, resp, helpers, ...) {
    gen <- character(length(formula))
    for (i in seq_along(formula)) {
        line_args <- c(list(i = resp[i], idt), helpers[[i]])
        gen[i] <- lines_wrap("generated_quantities", formula[[i]], line_args)
    }
    if (any(nzchar(gen))) {
        paste_rows("generated quantities {", collapse_rows(gen), "}")
    } else {
        NULL
    }
}

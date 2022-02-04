#' Create code blocks for the Stan model
#'
#' @param becauseformula defining model components
#' @export
create_blocks <- function(formula, ...) {
    UseMethod("create_blocks")
}

#' @export
create_blocks.default <- function(formula, ...) {

    functions <- paste("functions {", .loglik_stan, "}", sep = "\n")
    data <- create_data(formula, ....)
    transformed_data <- create_transformed_data(formula, ....)
    parameters <- create_parameters(formula, ....)
    transformed_parameters <- create_transformed_parameters(formula, ....)
    model <- create_model(formula, ....)
    generated_quantities <- create_generated_quantities(formula, ....)
    model_code <- paste(functions, data, transformed_data, parameters,
        transformed_parameters, model, generated_quantities,
        sep = "\n") # combine above text blocks
    model_code
}

#'
#' @export
create_data <- function(formula, ...) {

    continuous_distributions <- c("gaussian", "gamma") # ? Maybe define these globally
    discrete_distributions <- c("categorical", "binomial", "poisson") # ? use proper family objects later

    # Maybe not the best way...
    mtext <- paste(
        "int<lower=1> T; // number of time points",
        "int<lower=1> N; // number of individuals",
        "int<lower=1> C; // number of channels/response variables",
        "int<lower=1> K; // total number of covariates across all channels",
        "matrix[N, K] X[T]; // all covariates as an array of N x K matrices",
        "int<lower=0> D; // number of B-splines",
        "matrix[D, T] Bs; // B-spline basis matrix",
        sep = "\n")

    # loop over channels
    # TODO now using formula$resp and formula$families, should probably have
    # formula$resp[i]$family, formula$resp[i]$name or something
    for (i in seq_along(formula$resp)) {

        if (formula$families[i] %in% continuous_distributions) {
            type <- "real"
        } else {
            type <- "int"
        }

        # create response for channel i as T x N array
        y <- paste0(type, paste0(" response_", i), "[T,N];")
        mtext <- paste(mtext, y, sep = "\n")
        # Number of covariates in channel i
        mtext <- paste(mtext, paste0("int<lower=1> K_", i, ";"))
        # index vector of covariates related to channel i
        mtext <- paste(mtext, paste0("int J_", i, "[K_", i, "];"), sep = "\n")

        # TODO, need to add other distribution-specific components as well
        if (formula$families[i] == "categorical") {
            mtext <- paste(mtext,  paste0("int<lower=0> S_", i, ";"), sep = "\n")
        }
    }
    mtext <- paste("data {", mtext, "}", sep="\n") # TODO indentation?
    mtext
}

#'
#' @export
create_transformed_data <- function(formula, ...) {
    # TODO: Consider transforming X by centering, typically improves sampling efficiency
    # Should center over all time points?
    # TODO: Add vectors of zeros e.g. row_vector[K] zeros_K = rep_row_vector(0, K);
    # At least for categorical sequences
    mtext <- ""
    mtext <- paste("transformed data {", mtext, "}", sep="\n")
    mtext
}

#'
#' @export
create_parameters <- function(formula, ...) {

    mtext <- paste(
        "// Spline parameters",
        "// use noncentered parameterisation for sampling efficiency", sep = "\n")

    # Not the most efficient way but gets the job done
    for (i in seq_along(formula$resp)) {

        # do we want to handle intercept separately or as first beta?
        if (formula$families[i] == "categorical") {
            a_raw_term <- paste0("row_vector[D] a_raw_", i, "[S_", i, "- 1, K_", i, "];")
            # Do we want separate tau for K coefficients or K * (S - 1) coefficients?
            tau_term <- paste0("vector<lower=0>[K_", i, "] tau_", i, ";")
        } else {
            # TODO: need to add distribution-specific parameters i.e sigma for gaussian case
            a_raw_term <- paste0("row_vector[D] a_raw_", i, "[K_", i, "];")
            tau_term <- paste0("vector<lower=0>[K_", i, "] tau_", i, ";")
        }
        mtext <- paste(mtext, a_raw_term, tau_term, sep = "\n")
    }
    if (formula$splines$shrinkage) {
        mtext <- paste(mtext, "vector<lower=0>[D - 1] lambda; // shrinkage parameter", sep = "\n")
    }
    mtext <- paste("parameters {", mtext, "}", sep = "\n") # no indentation... TODO?
    mtext
}

#'
#' @export
create_transformed_parameters <- function(formula, ...) {

    mtext <- NULL
    # define variables
    for (i in seq_along(formula$resp)) {
        if (formula$families[i] == "categorical") { # Fix after definining families
            beta_term <- paste0("matrix[K_", i, ", S_", i, "] beta_", i, "[T];")
            a_term <- paste0("row_vector[D] a_", i, "[S_", i, " - 1, K_", i, "];")
        } else {
            beta_term <- paste0("vector[K_", i, "] beta_", i, "[T];")
            a_term <- paste0("row_vector[D] a_", i, "[K_", i, "];")
        }

        mtext <- paste(mtext, beta_term, a_term, sep = "\n")
    }

    for (i in seq_along(formula$resp)) {
        if (formula$families[i] == "categorical") {
            # N(0, 1) prior for first a, should be user-defined
            spline_term <- paste(
                # TODO define zeros_K_i in transformed data?
                paste0("for(s in 1:(S_", i, " - 1)) {"),
                paste0("  for(k in 1:K_", i, ") {"),
                paste0("    a[s, k, 1] = a_raw_", i, "[s, k, 1];"),
                "    for (i in 2:D) {",
                paste0("      a_", i, "[s, k, i] = a_", i, "[s, k, i-1] + a_raw_", i, "[s, k, i] * tau_", i, "[k] * lambda[i - 1];"),
                "    }",
                "    for(t in 1:T){",
                paste0("      beta_", i, "[t, k, s] = a_", i, "[s, k] * Bs[, t];"),
                "    }",
                "  }",
                "}",
                "for(t in 1:T){",
                paste0("  beta_", i, "[t, , S_", i, "] = zeros_K_", i, ";"),
                "}", sep = "\n")
        } else {

            spline_term <- paste(
                paste0("for(k in 1:K_", i, ") {"),
                paste0("  a_", i, "[k, 1] = a_raw_", i, "[k, 1];"),
                "  for (i in 2:D) {",
                paste0("    a_", i, "[k, i] = a_", i, "[k, i-1] + a_raw_", i, "[k, i] * tau_", i, "[k] * lambda[i - 1];"),
                "  }",
                "  for(t in 1:T){",
                paste0("    beta_", i, "[t, k] = a_", i, "[k] * Bs[, t];"),
                "  }",
                "}", sep = "\n")
        }
        mtext <- paste(mtext, spline_term, sep = "\n")
    }
    mtext <- paste("transformed parameters {", mtext, "}", sep="\n")
    mtext
}


#'
#' @export
create_model <- function(formula, ...) {
    # TODO: user-defined priors
    #
    mtext <- if (formula$splines$shrinkage) {
        "lambda ~ std_normal();  // prior for shrinkage terms"
    } else NULL
    mtext <- paste(mtext, paste0("tau_", 1:length(formula$families), " ~ cauchy(0, 1);", collapse = "\n"), sep = "\n")
    for (i in seq_along(formula$families)) {
        if (formula$families[i] == "categorical") {
            a_term <- paste(
                paste0("for(s in 1:(S_", i, " - 1)) {"),
                paste0("  for(k in 1:K_", i, ") {"),
                paste0("    a_raw_", i, "[s, k] ~ std_normal();"),
                "  }",
                "}", sep = "\n")
        } else {
            a_term <- paste(
                paste0("for(k in 1:K_", i, ") {"),
                paste0("  a_raw_", i, "[k] ~ std_normal();"),
                "}", sep = "\n")
        }
        mtext <- paste(mtext, a_term, sep = "\n")
    }
    for (i in seq_along(formula$families)) {
        if (formula$families[i] == "categorical") {
            likelihood_term <-
                paste0("  target += categorical_logit_glm_lupmf(response_", i, "[t] | X[t][, J_",i , "], zeros_S_", i, ", beta_", i, "[t]);")
        } else {
            stop("Likelihoods not yet defined for non-categorical data")
        }
        mtext <- paste(mtext, "for(t in 1:T) {", likelihood_term, "}", sep = "\n")
    }

    mtext <- paste("model {", mtext, "}", sep = "\n")
    mtext
}

#'
#' @export
create_generated_quantities <- function(formula, ...) {
    # TODO? For simulating predictions, we probably need to do it in R as we
    # need to create new covariate matrix at each time point given previous samples
    # and we might be interested in multiple things based on the results anyway
    mtext <- ""
    mtext <- paste("generated quantities {", mtext, "\n}")
    mtext
}

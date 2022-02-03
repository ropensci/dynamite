## Commented out for a now in order to avoid duplicate function names, needs rework anyway

#' #' Create code blocks for the Stan model
#' #'
#' #' @param becauseformula defining model components
#' #' @export
#' create_blocks <- function(formula, ...) {
#'     UseMethod("create_blocks")
#' }
#'
#' #' @export
#' create_blocks.default <- function(formula, ...) {
#'
#'     functions <- paste("functions {", .loglik_stan, "}", sep = "\n")
#'     data <- create_data(formula, ....)
#'     transformed_data <- create_transformed_data(formula, ....)
#'     parameters <- create_parameters(formula, ....)
#'     transformed_parameters <- create_transformed_parameters(formula, ....)
#'     model <- create_model(formula, ....)
#'     generated_quantities <- create_generated_quantities(formula, ....)
#'     model_code <- paste(functions, data, transformed_data, parameters,
#'         transformed_parameters, model, generated_quantities,
#'         sep = "\n") # combine above text blocks
#'     model_code
#' }
#'
#' #'
#' #' @export
#' create_data <- function(formula, ...) {
#'
#'     continuous_distributions <- c("gaussian", "gamma") # ? Maybe define these globally
#'     discrete_distributions <- c("categorical", "binomial", "poisson") # ? use proper family objects later
#'
#'     # Maybe not the best way...
#'     mtext <- paste(
#'         "int<lower=1> T; // number of time points",
#'         "int<lower=1> N; // number of individuals",
#'         "int<lower=1> M; // number of hidden states",
#'         "int<lower=1> C; // number of channels/response variables",
#'         "int<lower=1> K; // total number of covariates accross all channels",
#'         "int<lower=0, upper=T> T_fixed_start; // fix to state one",
#'         "int<lower=0, upper=T> T_fixed_end; // fix to state M",
#'         "matrix[N, K] X[T] // all covariates as an array of N x K matrices", sep = "\n")
#'
#'     # loop over channels
#'     # TODO now using formula$resp and formula$families, should probably have
#'     # formula$resp[i]$family, formula$resp[i]$name or something
#'     for (i in seq_along(formula$resp)) {
#'
#'         if (formula$families[i] %in% continuous_distributions) {
#'             type <- "real"
#'         } else {
#'             type <- "int"
#'         }
#'
#'         # create response for channel i as T x N array
#'         y <- paste0(type, paste0(" response_", i), "[T,N];")
#'         mtext <- paste(mtext, y, sep = "\n")
#'         # index vector of covariates related to channel i
#'         mtext <- paste(mtext, paste0("int J_", i, "[K_", i, "]"), sep = "\n")
#'
#'         # TODO, need to add other distribution-specific components as well
#'         if (formula$families[i] == "categorical") {
#'             mtext <- paste(mtext,  paste0("int<lower=0> S_", i), sep = "\n")
#'         }
#'     }
#'     mtext <- paste("data {", mtext, "}", sep="\n") # TODO indentation?
#'     mtext
#' }
#'
#' #'
#' #' @export
#' create_transformed_data <- function(formula, ...) {
#'     # Consider transforming X by centering and maybe scaling, or do it in R?
#'     mtext <- ""
#'     mtext <- paste("transformed data {", mtext, "}", sep="\n")
#'     mtext
#' }
#'
#' #'
#' #' @export
#' create_parameters <- function(formula, ...) {
#'
#'     #TODO
#'     mtext <- paste(
#'         "// state-specific parameters",
#'         "// use noncentered parameterisation for sampling efficiency", sep = "\n")
#'
#'     # Not the most efficient way but gets the job done
#'     for (i in seq_along(formula$resp)) {
#'
#'         # do we want to handle intercept separately or as first beta?
#'         if (formula$families[i] == "categorical") {
#'             # last column not sampled but fixed to 0 for identifiability
#'             # intercept_term <-
#'             #    paste0("vector[S_", i, " - 1] intercept_raw_", i, "[M];")
#'             beta_term <- paste0("matrix[K_", i, ", S_", i, " - 1] beta_raw_", i, "[M];")
#'             sigma_term <- paste0("matrix[K_", i, ", S_", i, " - 1] sigma_beta_", i, ";")
#'         } else {
#'             #intercept_term <- paste0("real intercept_raw_", i, "[M];")
#'             beta_term <- paste0("vector[K_", i, "] beta_raw_", i, "[M];")
#'             sigma_term <- paste0("vector[K_", i, "] sigma_beta_", i, ";")
#'         }
#'         # # standard deviations of the random walk priors for parameters
#'         # # here same sigma is used for for all coefficients, probably not optimal if scales differ.
#'         # # and some standardizations are needed for X in this case
#'         # mtext <- paste(mtext, intercept_term, beta_term,
#'         #     paste0("real<lower=0> sigma_intercept_", i, ";"),
#'         #     paste0("real<lower=0> sigma_beta_", i, ";"),
#'         #     sep = "\n")
#'
#'         # unique sigmas for all
#'         mtext <- paste(mtext, beta_term, sigma_term,
#'             sep = "\n")
#'         # TODO, need to add distribution-specific parameters
#'     }
#'     # Strict left-to-right model
#'     mtext <- paste(mtext,
#'         "vector<lower=0,upper=1>[M-1] A; // diagonal of the transition matrix",
#'         sep = "\n")
#'     mtext <- paste("parameters {", mtext, "}", sep = "\n") # no indentation... TODO?
#'     mtext
#' }
#'
#' #'
#' #' @export
#' create_transformed_parameters <- function(formula, ...) {
#'     #TODO
#'     # currently we assume common states for all individuals
#'     # if not, log_py should be defined as matrix[M, T] log_py[N]
#'     # with obvious changes to log_py computations and separate
#'     # forward algorithm run for each individual
#'     mtext <- paste(
#'         "// log-likelihoods p(y_t | z_t)",
#'         "matrix[M, T] log_py = rep_matrix(0, M, T);", sep = "\n")
#'
#'     # define variables
#'     for (i in seq_along(formula$resp)) {
#'
#'         if (formula$families[i] == "categorical") { # Fix after definining families
#'             # intercept_term <-
#'             #    paste0("vector[S_", i, "] intercept_", i, "[M]; // last column as 0 for identifiability")
#'             beta_term <- paste0("matrix[K_", i, ", S_", i, "] beta_", i, "[M];")
#'         } else {
#'             #intercept_term <- paste0("real intercept_", i, "[M];")
#'             beta_term <- paste0("vector[K_", i, "] beta_", i, "[M];")
#'         }
#'
#'         mtext <- paste(mtext, intercept_term, beta_term,
#'             sep = "\n")
#'     }
#'
#'     mtext <- paste(mtext,
#'         "// transition matrix for forward algorithm",
#'         "matrix[M, M] log_A = inf_mat;",
#'         "for(i in 1:(M-1)) {",
#'         "  log_A[i, i] = log(A[i]);",
#'         "  log_A[i, i + 1] = log(1 - A[i]);",
#'         "}",
#'         "log_A[M, M] = 0;",
#'         sep = "\n")
#'
#'     # random walks
#'     for (i in seq_along(formula$resp)) {
#'         if (formula$families[i] == "categorical") {
#'             # N(0, 1) prior for first state, should be user-defined
#'             rw <- paste(
#'                 # paste0("intercept_", i, "[1] = append_row(intercept_raw_", i, "[1], 0);"),
#'                 # TODO define zeros_K_i in transformed data?
#'                 paste0("beta_", i, "[1] = append_col(beta_raw_", i, "[1], zeros_K_", i,");"),
#'                 "for(m in 2:M) {",
#'                 #paste0("  intercept_", i, "[m] = intercept_", i, "[m-1] + ",
#'                 #    "sigma_intercept_", i, " .* append_row(intercept_raw_", i, "[m], 0);"),
#'                 paste0("  beta_", i, "[m] = beta_", i, "[m - 1] + sigma_beta_", i,
#'                     " .* append_col(beta_raw_", i, "[m], zeros_K_", i, ");"),
#'                 "}",
#'                 sep = "\n")
#'         } else {
#'
#'             rw <- paste(
#'                 #paste0("intercept_", i, "[1] = intercept_raw_", i, "[1];"),
#'                 paste0("beta_", i, "[1] = beta_raw_", i, "[1]; "),
#'                     "for(m in 2:M) {",
#'                 #paste0("  intercept_", i, "[m] = intercept_", i, "[m-1] + sigma_intercept_", i, " .* intercept_raw_", i, "[m];"),
#'                 paste0("  beta_", i, "[m] = beta_", i, "[m-1] + sigma_beta_", i, " .* beta_raw_", i, "[m];"),
#'                     "}",
#'                 sep = "\n")
#'         }
#'         mtext <- paste(mtext, rw, sep = "\n")
#'     }
#'     # likelihood terms
#'     distributions <- character(length(formula$resp))
#'     for (i in seq_along(formula$resp)) {
#'         # TODO get_stan_distribution(formula$families[i], i)?
#'
#'         if (formula$families[i] == "categorical") {
#'             distributions[i] <-
#'                # all individuals at once
#'                "categorical_logit_glm_lupmf(response_i[t] | X[t][,J_i], zeros_S_i, beta_i[m])"
#'         } else {
#'             #distribution <- "normal_lpdf(response_i[t] | intercept_i[m] + X_i[t] * beta_i[m], sigma_i)"
#'             stop("Only categorical sequences are currently supported.")
#'         }
#'         ll <- paste(
#'             "for(t in 1:T) {",
#'             "  for(m in 1:M) {",
#'             "    log_py[m, t] += ", paste0(distributions, collapse = " +\n"), ";",
#'             "  }",
#'             "}",
#'             sep = "\n")
#'     }
#'     mtext <- paste(mtext, ll, sep = "\n")
#'
#'     mtext <- paste(mtext,
#'         "// fix first and last states",
#'         "log_py[2:M, 2:T_fixed_start] += negative_infinity();",
#'         "log_py[1:(M - 1), (T - T_fixed_end + 1):T] += negative_infinity();",
#'         sep = "\n")
#'
#'
#'     mtext <- paste("transformed parameters {", mtext, "}", sep="\n")
#'     mtext
#' }
#'
#' #'
#' #' @export
#' create_model <- function(formula, ...) {
#'     #TODO user-defined priors
#'
#'     # Prior for staying in same state
#'     # default prior could be something like this:
#'     mtext <- "A ~ beta(T, M);  // prior average holding time (T + M) / M"
#'
#'     # Some (default) priors for random walk standard deviations
#'     mtext <- paste(mtext,
#'         "for(m in 1:M) {",
#'         paste0("  sigma_beta_", 1:length(formula$resp), " ~ std_normal();", collapse = "\n"),
#'         "}",
#'         sep = "\n")
#'
#'     # random walks
#'     for (i in seq_along(formula$resp)) {
#'         # Fix after definining families
#'         if (formula$families[i] == "categorical") {
#'             rw <- paste(
#'                 paste0("to_vector(beta_raw_", i, "[m]) ~ std_normal();"),
#'                # paste0("intercept_raw_", i, "[m] ~ std_normal();"),
#'                 sep = "\n")
#'
#'         } else {
#'             rw <- paste(
#'                 paste0("beta_raw_", i, "[m] ~ std_normal();"),
#'                # paste0("intercept_raw_", i, "[m] ~ std_normal();"),
#'                 sep = "\n")
#'         }
#'         mtext <- paste(mtext, rw, sep = "\n")
#'     }
#'
#'
#'     # defined in the functions block (see loglik_stan.R)
#'     mtext <- paste(mtext,
#'         "target += loglik(log_py, M, T, log_A);",
#'         sep = "\n")
#'
#'     mtext <- paste("model {", mtext, "}", sep = "\n")
#'     mtext
#' }
#'
#' #'
#' #' @export
#' create_generated_quantities <- function(formula, ...) {
#'     #TODO?
#'     mtext <- ""
#'     mtext <- paste("generated quantities {", mtext, "\n}")
#'     mtext
#' }

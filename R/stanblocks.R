#' Create code blocks for the Stan model
#'
#' @param becauseformula defining model components
#' @export
create_blocks <- function(formula, ...) {
    UseMethod("becauseformula")
}

#'
#'@export
create_blocks.default <- function(formula, ...) {
    data <- create_data(formula, ....)
    transformed_data <- create_transformed_data(formula, ....)
    parameters <- create_parameters(formula, ....)
    transformed_parameters <- create_transformed_parameters(formula, ....)
    model <- create_model(formula, ....)
    generated_quantities <- create_generated_quantities(formula, ....)
    model_code <- paste(data, transformed_data, parameters,
        transformed_parameters, mdoel, generated_quantities,
        sep = "\n") # combine above text blocks
    model_code
}

#'
#'@export
create_data <- function(formula, ...) {

    continuous_distributions <- c("gaussian", "gamma") # ? Maybe define these globally
    discrete_distributions <- c("categorical", "binomial", "poisson") # ?

    # Maybe not the best way...
    mtext <- paste(
        "int<lower=1> T; // number of time points",
        "int<lower=1> N; // number of individuals",
        "int<lower=1> M; // number of hidden states",
        "int<lower=0, upper=T> T_fixed_start;",
        "int<lower=0, upper=T> T_fixed_end;",
        "int<lower=1> C; // number of channels", sep = "\n")

    for(i in seq_along(formula$resp)) {

        if(formula$families[i] %in% continuous_distributions) {
            type <- "real"
        } else {
            type <- "int"
        }

        y <- paste0(type, paste0(" response_", i), "[T,N];")
        mtext <- paste(mtext, y, sep = "\n")
        mtext <- paste(mtext,
            paste0("int<lower=0> K", i),
            paste0("matrix[N, K_", i, "] X_", i, "[T];"),
            sep = "\n")

        # TODO, need to add distribution-specific components as well, e.g, number of symbols
    }
    mtext <- paste("data {", mtext, "\n}") # no indentation... TODO?
    mtext
}

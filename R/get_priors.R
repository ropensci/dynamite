#' @export
get_priors <- function(formula, data, group, time) {
  btvcmfit(formula, data, group, time, debug = list(no_compile = TRUE))$priors
}
#' @export
get_code <- function(formula, data, group, time) {
    btvcmfit(formula, data, group, time, debug = list(no_compile = TRUE, model_code = TRUE))$model_code
}



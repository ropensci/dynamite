#' @export
get_priors <- function(formula, data, group, time) {

    if (!is.character(group)) {
        group <- deparse(substitute(group))
    }
    if (!is.character(time)) {
        time <- deparse(substitute(time))
    }
    btvcmfit(formula, data, group, time, debug = list(no_compile = TRUE))$priors
}
#' @export
get_code <- function(formula, data, group, time) {
    if (!is.character(group)) {
        group <- deparse(substitute(group))
    }
    if (!is.character(time)) {
        time <- deparse(substitute(time))
    }
    btvcmfit(formula, data, group, time, debug = list(no_compile = TRUE, model_code = TRUE))$model_code
}



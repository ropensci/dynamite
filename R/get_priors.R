#' @export
get_priors <- function(formula, data, group, time) {
  out <- do.call(
    "dynamitefit",
    list(
      formula = formula,
      data = data,
      group = substitute(group),
      time = substitute(time),
      debug = list(no_compile = TRUE)
    )
  )
  out$priors
}
#' @export
get_code <- function(formula, data, group, time) {
  out <- do.call(
    "dynamitefit",
    list(
      formula = formula,
      data = data,
      group = substitute(group),
      time = substitute(time),
      debug = list(no_compile = TRUE, model_code = TRUE)
    )
  )
  out$model_code
}
#' @export
get_data <- function(formula, data, group, time) {
  out <- do.call(
    "dynamitefit",
    list(
      formula = formula,
      data = data,
      group = substitute(group),
      time = substitute(time),
      debug = list(no_compile = TRUE, model_data = TRUE)
    )
  )
  out$model_data
}

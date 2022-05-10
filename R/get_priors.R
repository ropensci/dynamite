#' Get prior definitions for a dynamite model
#'
#' TODO description
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula.
#'   See [dynamiteformula()].
#' @param data \[`data.frame`]\cr The data frame containing the variables in
#'   the model.
#' @param group \[`character(1)`]\cr A column name of `data` that denotes the
#'   unique groups.
#' @param time \[`character(1)`]\cr A column name of `data` that denotes the
#'   time axis.
#'
#' @export
get_priors <- function(dformula, data, group, time) {
  out <- do.call(
    "dynamite",
    list(
      dformula = dformula,
      data = data,
      group = substitute(group),
      time = substitute(time),
      debug = list(no_compile = TRUE)
    )
  )
  out$priors
}


#' TODO documentation
#' @inheritParams get_priors
#' @export
get_code <- function(dformula, data, group, time) {
  out <- do.call(
    "dynamite",
    list(
      dformula = dformula,
      data = data,
      group = substitute(group),
      time = substitute(time),
      debug = list(no_compile = TRUE, model_code = TRUE)
    )
  )
  out$model_code
}

#' TODO documentation
#' @inheritParams get_priors
#' @export
get_data <- function(dformula, data, group, time) {
  out <- do.call(
    "dynamite",
    list(
      dformula = dformula,
      data = data,
      group = substitute(group),
      time = substitute(time),
      debug = list(no_compile = TRUE, model_data = TRUE)
    )
  )
  out$model_data
}

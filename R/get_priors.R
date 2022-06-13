#' Get Prior Definitions for a Dynamite Model
#'
#' Extracts the priors used in the dynamite model as a data frame. You
#' can then alter the priors by changing the contents of the `prior` column and
#' supplying this data frame to `dynamite` function using the argument `priors`.
#'
#' Note that the prior for the intercept term `alpha` is actually defined
#' in a centered form, so the prior is related to the `alpha` when the
#' covariates at the first time point are centered around their means. In other
#' words, the prior is defined for `alpha + x_m * gamma` where `x_m` is vector
#' of covariate means and gamma contains the corresponding coefficients (`beta`
#' and `delta_1`). If you want to use prior directly on `alpha`, remove
#' intercept from the formula and add a dummy covariate consisting of ones to
#' the model.
#'
#' @note Warning! Currently the structure of the `prior` argument is not
#' checked in the `dynamite` function, and it is assumed that no rows are
#' deleted and the ordering of the rows is not changed (i.e. only `prior`
#' column is altered).
#'
#' @param x \[`dynamiteformula` or `dynamitefit`]\cr The model formula or
#'   existing `dynamitefit` object. See [dynamiteformula()] and [dynamite()].
#' @param data \[`data.frame`]\cr The data frame containing the variables in
#'   the model.
#' @param group \[`character(1)`]\cr A column name of `data` that denotes the
#'   unique groups.
#' @param time \[`character(1)`]\cr A column name of `data` that denotes the
#'   time axis.
#' @param ... Ignored.
#' @rdname get_priors
#' @export
#' @examples
#'
#' d <- data.frame(y = rnorm(10), x = 1:10, time = 1:10, id = 1)
#' get_priors(obs(y ~ x, family = gaussian()),
#'   data = d, time = "time", group = "id")
#' @srrstats {RE2.3} *Where applicable, Regression Software should enable data to be centred (for example, through converting to zero-mean equivalent values; or to z-scores) or offset (for example, to zero-intercept equivalent values) via additional parameters, with the effects of any such parameters clearly documented and tested.*
#' @srrstats {BS5.2} *Bayesian Software should either return the input function or prior distributional specification in the return object; or enable direct access to such via additional functions which accept the return object as single argument.*
get_priors <- function(x, data, group, time, ...) {
  UseMethod("get_priors", x)
}
#' @method get_priors dynamiteformula
#' @rdname get_priors
#' @export
get_priors.dynamiteformula <- function(x, data, group, time, ...) {
  out <- do.call(
    "dynamite",
    list(
      dformula = x,
      data = data,
      group = substitute(group),
      time = substitute(time),
      debug = list(no_compile = TRUE),
      ...
    )
  )
  out$priors
}
#' @method get_priors dynamitefit
#' @rdname get_priors
#' @export
get_priors.dynamitefit <- function(x, ...) {
  x$priors
}

#' Extract the Stan Code of the Dynamite Model
#'
#' Returns the Stan code of the model. Mostly useful for debugging or for
#' building customized version of the model.
#'
#' @inheritParams get_priors
#' @rdname get_code
#' @export
#' @examples
#'
#' d <- data.frame(y = rnorm(10), x = 1:10, time = 1:10, id = 1)
#' cat(get_code(obs(y ~ x, family = gaussian()),
#'   data = d, time = "time", group = "id"))
get_code <- function(x, data, group, time, ...) {
  UseMethod("get_code", x)
}
#' @rdname get_code
#' @export
get_code.dynamiteformula <- function(x, data, group, time, ...) {
  out <- do.call(
    "dynamite",
    list(
      dformula = x,
      data = data,
      group = substitute(group),
      time = substitute(time),
      debug = list(no_compile = TRUE, model_code = TRUE),
      ...
    )
  )
  out$model_code
}
#' @rdname get_code
#' @export
get_code.dynamitefit <- function(x, ...) {
  x$stanfit@stanmodel
}

#' Extract the Model Data of the Dynamite Model
#'
#' Returns the input data to Stan model. Mostly useful for debugging.
#'
#' @inheritParams get_priors
#' @rdname get_data
#' @export
#' @examples
#'
#' d <- data.frame(y = rnorm(10), x = 1:10, time = 1:10, id = 1)
#' str(get_data(obs(y ~ x, family = gaussian()),
#'   data = d, time = "time", group = "id"))
get_data <- function(x, data, group, time, ...) {
  UseMethod("get_data", x)
}
#' @rdname get_data
#' @export
get_data.dynamiteformula <- function(x, data, group, time, ...) {
  out <- do.call(
    "dynamite",
    list(
      dformula = x,
      data = data,
      group = substitute(group),
      time = substitute(time),
      debug = list(no_compile = TRUE, sampling_vars = TRUE),
      ...
    )
  )
  out$sampling_vars
}
#' @rdname get_data
#' @export
get_data.dynamitefit <- function(x, ...) {

  out <- do.call(
    "dynamite",
    list(
      dformula = x$formula,
      data = x$data,
      group = x$group_var,
      time = x$time_var,
      priors = x$priors,
      debug = list(no_compile = TRUE, sampling_vars = TRUE)
    )
  )
  out$sampling_vars
}

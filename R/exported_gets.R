#' Get Prior Definitions of a Dynamite Model
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
#' @note Only the `prior` column of the output should be altered when defining
#' the user-defined priors for the `dynamite`.
#'
#' @export
#' @rdname get_priors
#' @param x \[`dynamiteformula` or `dynamitefit`]\cr The model formula or an
#'   existing `dynamitefit` object. See [dynamiteformula()] and [dynamite()].
#' @inheritParams dynamite
#' @param ... Ignored.
#' @return A `data.frame` containing the prior definitions.
#' @srrstats {BS5.2} Provides access to the prior definitions of the model.
#' @examples
#' d <- data.frame(y = rnorm(10), x = 1:10, time = 1:10, id = 1)
#' get_priors(obs(y ~ x, family = "gaussian"),
#'   data = d, time = "time", group = "id"
#' )
#'
get_priors <- function(x, ...) {
  UseMethod("get_priors", x)
}

#' @rdname get_priors
#' @export
get_priors.dynamiteformula <- function(x, data, group = NULL, time, ...) {
  out <- do.call(
    "dynamite",
    list(
      dformula = x,
      data = data,
      group = group,
      time = time,
      debug = list(no_compile = TRUE),
      ...
    )
  )
  out$priors
}

#' @rdname get_priors
#' @export
get_priors.dynamitefit <- function(x, ...) {
  x$priors
}

#' Extract the Stan Code of the Dynamite Model
#'
#' Returns the Stan code of the model. Mostly useful for debugging or for
#' building a customized version of the model.
#'
#' @export
#' @rdname get_code
#' @inheritParams get_priors.dynamiteformula
#' @return A Stan model code as a `character` string.
#' @examples
#' d <- data.frame(y = rnorm(10), x = 1:10, time = 1:10, id = 1)
#' cat(get_code(obs(y ~ x, family = "gaussian"),
#'   data = d, time = "time", group = "id"
#' ))
#' # same as
#' cat(dynamite(obs(y ~ x, family = "gaussian"),
#'   data = d, time = "time", group = "id",
#'   debug = list(model_code = TRUE, no_compile = TRUE)
#' )$model_code)
#'
get_code <- function(x, ...) {
  UseMethod("get_code", x)
}

#' @rdname get_code
#' @export
get_code.dynamiteformula <- function(x, data, group = NULL, time, ...) {
  out <- do.call(
    "dynamite",
    list(
      dformula = x,
      data = data,
      group = group,
      time = time,
      debug = list(no_compile = TRUE, model_code = TRUE),
      ...
    )
  )
  out$model_code
}

#' @rdname get_code
#' @export
get_code.dynamitefit <- function(x, ...) {
  x$stanfit@stanmodel@model_code
}

#' Extract the Model Data of the Dynamite Model
#'
#' Returns the input data to the Stan model. Mostly useful for debugging.
#'
#' @export
#' @rdname get_data
#' @inheritParams get_priors.dynamiteformula
#' @return A `list` containing the input data to Stan.
#' @examples
#' d <- data.frame(y = rnorm(10), x = 1:10, time = 1:10, id = 1)
#' str(get_data(obs(y ~ x, family = "gaussian"),
#'   data = d, time = "time", group = "id"
#' ))
#'
get_data <- function(x, ...) {
  UseMethod("get_data", x)
}

#' @rdname get_data
#' @export
get_data.dynamiteformula <- function(x, data, group = NULL, time, ...) {
  out <- do.call(
    "dynamite",
    list(
      dformula = x,
      data = data,
      group = group,
      time = time,
      debug = list(no_compile = TRUE, sampling_vars = TRUE),
      ...
    )
  )
  out$stan$sampling_vars
}

#' @rdname get_data
#' @export
get_data.dynamitefit <- function(x, ...) {
  x$stan$sampling_vars
}

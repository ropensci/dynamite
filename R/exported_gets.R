#' Get Prior Definitions of a Dynamite Model
#'
#' Extracts the priors used in the dynamite model as a data frame. You
#' can then alter the priors by changing the contents of the `prior` column and
#' supplying this data frame to `dynamite` function using the argument
#' `priors`.
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
#' @family fitting
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
get_priors.dynamiteformula <- function(x, data, time, group = NULL, ...) {
  out <- do.call(
    "dynamite",
    list(
      dformula = x,
      data = data,
      time = time,
      group = group,
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
#' @family output
#' @rdname get_code
#' @inheritParams get_priors.dynamiteformula
#' @param blocks \[`character()`]\cr Stan block names to extract. If `NULL`,
#'   extracts the full model code.
#' @return The Stan model blocks as a `character` string.
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
get_code.dynamiteformula <- function(x, data, time,
                                     group = NULL, blocks = NULL, ...) {
  out <- do.call(
    "dynamite",
    list(
      dformula = x,
      data = data,
      time = time,
      group = group,
      debug = list(no_compile = TRUE, model_code = TRUE),
      ...
    )
  )
  get_code_(out$model_code, blocks)
}

#' @rdname get_code
#' @export
get_code.dynamitefit <- function(x, blocks = NULL, ...) {
  stopifnot_(
    !is.null(x$stanfit),
    "No Stan model fit is available."
  )
  get_code_(x$stanfit@stanmodel@model_code, blocks)
}

#' Internal Stan Code Block Extraction
#'
#' @param x \[`character(1L)`]\cr The Stan model code string.
#' @param blocks \[`character`]\cr Stan block names to extract. If `NULL`,
#'   extracts the full model code.
#' @noRd
get_code_ <- function(x, blocks = NULL) {
  if (is.null(blocks)) {
    return(x)
  }
  stopifnot_(
    checkmate::test_character(blocks, null.ok = TRUE),
    "Argument {.arg blocks} must be a {.cls character} vector or NULL."
  )
  block_names <- c(
    "data",
    "transformed data",
    "parameters",
    "transformed parameters",
    "model",
    "generated quantities"
  )
  invalid_blocks <- !blocks %in% block_names
  stopifnot_(
    all(!invalid_blocks),
    c(
      "Invalid Stan blocks provided: {cs(blocks[invalid_blocks])}",
      `i` = "Argument {.arg blocks} must be NULL or a subset of
             {cs(paste0(\"'\", block_names, \"'\"))}."
    )
  )
  x <- strsplit(x, "\n")[[1L]]
  block_rows <- paste0(block_names, " {")
  block_start <- which(x %in% block_rows)
  block_end <- c(block_start[-1L] - 1L, length(x))
  names(block_start) <- names(block_end) <- block_names
  out <- ""
  for (block in blocks) {
    out <- c(
      out,
      x[block_start[block]:block_end[block]]
    )
  }
  paste_rows(out, .parse = FALSE)
}

#' Extract the Model Data of the Dynamite Model
#'
#' Returns the input data to the Stan model. Mostly useful for debugging.
#'
#' @export
#' @family output
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
get_data.dynamiteformula <- function(x, data, time, group = NULL, ...) {
  out <- dynamite(
    dformula = x,
    data = data,
    time = time,
    group = group,
    debug = list(no_compile = TRUE, sampling_vars = TRUE),
    ...
  )
  out$stan$sampling_vars
}

#' @rdname get_data
#' @export
get_data.dynamitefit <- function(x, ...) {
  x$stan$sampling_vars
}

#' Get Parameter Dimensions of the Dynamite Model
#'
#' Extracts the names and dimensions of all parameters used in the
#' `dynamite` model. See also [dynamite::get_parameter_types()] and
#' [dynamite::get_parameter_names()]. The returned dimensions match those of
#' the `stanfit` element of the `dynamitefit` object. When applied to
#' `dynamiteformula` objects, the model is compiled and sampled for 1 iteration
#' to get the parameter dimensions.
#'
#' @param x \[`dynamitefit`, `dynamiteformula`]\cr A model fit object, or
#'   a model formula object.
#' @param ... Ignored
#' @return A named list with all parameter dimensions of the input model.
#' @export
#' @family output
#' @examples
#' get_parameter_dims(multichannel_example_fit)
get_parameter_dims <- function(x, ...) {
  UseMethod("get_parameter_dims", x)
}

#' @rdname get_parameter_dims
#' @export
get_parameter_dims.dynamiteformula <- function(x, data, time,
                                               group = NULL, ...) {
  # use capture.output to make Stan silent, maybe
  utils::capture.output(
    utils::capture.output(
      out <- dynamite(
        dformula = x,
        data = data,
        time = time,
        group = group,
        algorithm = "Fixed_param",
        chains = 1,
        iter = 1,
        refresh = -1,
        backend = "rstan",
        verbose_stan = FALSE,
        ...
      ),
      type = "message"
    ),
    type = "output"
  )
  get_parameter_dims(out)
}

#' @rdname get_parameter_dims
#' @export
get_parameter_dims.dynamitefit <- function(x, ...) {
  stopifnot_(
    !is.null(x$stanfit),
    "No Stan model fit is available."
  )
  draws <- rstan::extract(
    x$stanfit,
    pars = "lp__",
    include = FALSE
  )
  lapply(
    draws,
    function(y) {
      d <- dim(y)
      ifelse_(length(d) == 1L, 1L, d[-1L])
    }
  )
}


#' Get Parameter Types of the Dynamite Model
#'
#' Extracts all parameter types of used in the `dynamitefit` object. See
#' [dynamite::as.data.frame.dynamitefit()] for explanations of different types.
#'
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A `character` vector with all parameter types of the input model.
#' @export
#' @family output
#' @examples
#' get_parameter_types(multichannel_example_fit)
#'
get_parameter_types <- function(x, ...) {
  UseMethod("get_parameter_types", x)
}

#' @rdname get_parameter_types
#' @export
get_parameter_types.dynamitefit <- function(x, ...) {
  types <- c(
    "alpha", "beta", "delta", "tau", "tau_alpha", "xi",
    "sigma_nu", "sigma", "phi", "nu", "lambda", "sigma_lambda",
    "psi", "tau_psi", "corr", "corr_psi", "corr_nu",
    "omega", "omega_alpha", "omega_psi"
  )
  d <- as.data.table(x, types =  types)
  unique(d$type)
}

#' Get Parameter Names of the Dynamite Model
#'
#' Extracts all parameter names of used in the `dynamitefit` object.
#'
#' The naming of parameters generally follows style where the name starts with
#' the parameter type (e.g. beta for time-invariant regression coefficient),
#' followed by underscore and the name of the response variable, and in case of
#' time-invariant, time-varying or random effect, the name of the predictor. An
#' exception to this is spline coefficients omega, which also contain the
#' number denoting the knot number.
#'
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param types \[`character()`]\cr Extract only names of parameter of a
#'   certain type. See [dynamite::get_parameter_types()].
#' @param ... Ignored.
#' @return A `character` vector with parameter names of the input model.
#' @export
#' @family output
#' @examples
#' get_parameter_names(multichannel_example_fit)
#'
get_parameter_names <- function(x, types = NULL, ...) {
  UseMethod("get_parameter_names", x)
}

#' @rdname get_parameter_names
#' @export
get_parameter_names.dynamitefit <- function(x, types = NULL, ...) {
  if (is.null(types)) {
    types <- get_parameter_types(x)
  }
  d <- as.data.table(x, types =  types)
  unique(d$parameter)
}

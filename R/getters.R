#' Get Prior Definitions of a \pkg{dynamite} Model
#'
#' Extracts the priors used in the `dynamite` model as a data frame. You
#' can then alter the priors by changing the contents of the `prior` column and
#' supplying this data frame to `dynamite` function using the argument
#' `priors`. See vignettes for details.
#'
#' @note Only the `prior` column of the output should be altered when defining
#' the user-defined priors for `dynamite`.
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
#' data.table::setDTthreads(1) # For CRAN
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
  out <- dynamite(
    dformula = x,
    data = data,
    time = time,
    group = group,
    debug = list(no_compile = TRUE),
    ...
  )
  out$priors
}

#' @rdname get_priors
#' @export
get_priors.dynamitefit <- function(x, ...) {
  x$priors
}

#' Extract the Stan Code of the \pkg{dynamite} Model
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
#' data.table::setDTthreads(1) # For CRAN
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
  out <- dynamite(
    dformula = x,
    data = data,
    time = time,
    group = group,
    debug = list(no_compile = TRUE, model_code = TRUE),
    ...
  )
  get_code_(out$model_code, blocks)
}

#' @rdname get_code
#' @export
get_code.dynamitefit <- function(x, blocks = NULL, ...) {
  if (is.null(x$stanfit)) {
    out <- dynamite(
      dformula = eval(formula(x)),
      data = x$data,
      time = x$time_var,
      group = x$group_var,
      backend = x$backend,
      debug = list(no_compile = TRUE, model_code = TRUE),
      verbose = FALSE,
      ...
    )$model_code
  } else {
    out <- get_model_code(x$stanfit)
  }
  get_code_(out, blocks)
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

#' Extract the Model Data of the \pkg{dynamite} Model
#'
#' Returns the input data to the Stan model. Mostly useful for debugging.
#'
#' @export
#' @family output
#' @rdname get_data
#' @inheritParams get_priors.dynamiteformula
#' @return A `list` containing the input data to Stan.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
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
    debug = list(no_compile = TRUE, stan_input = TRUE, model_code = FALSE),
    ...
  )
  out$stan_input$sampling_vars
}

#' @rdname get_data
#' @export
get_data.dynamitefit <- function(x, ...) {
  if (!is.null(x$stan_input)) {
    return(x$stan_input$sampling_vars)
  }
  out <- dynamite(
    dformula = eval(formula(x)),
    data = x$data,
    time = x$time_var,
    group = x$group_var,
    priors = x$priors,
    backend = x$backend,
    interval = x$interval,
    debug = list(no_compile = TRUE, stan_input = TRUE, model_code = FALSE),
    verbose = FALSE,
    ...
  )
  out$stan_input$sampling_vars
}

#' Get Parameter Dimensions of the \pkg{dynamite} Model
#'
#' Extracts the names and dimensions of all parameters used in the
#' `dynamite` model. See also [get_parameter_types()] and
#' [get_parameter_names()]. The returned dimensions match those of
#' the `stanfit` element of the `dynamitefit` object. When applied to
#' `dynamiteformula` objects, the model is compiled and sampled for 1 iteration
#' to get the parameter dimensions.
#'
#' @rdname get_parameter_dims
#' @inheritParams get_priors.dynamiteformula
#' @return A named list with all parameter dimensions of the input model.
#' @export
#' @family output
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' get_parameter_dims(multichannel_example_fit)
#'
get_parameter_dims <- function(x, ...) {
  UseMethod("get_parameter_dims", x)
}

#' @rdname get_parameter_dims
#' @export
get_parameter_dims.dynamiteformula <- function(x, data, time,
                                               group = NULL, ...) {
  out <- dynamite(
    dformula = x,
    data = data,
    time = time,
    group = group,
    debug = list(no_compile = TRUE, stan_input = TRUE, model_code = FALSE),
    ...
  )
  get_parameter_dims(out, ...)
}

#' @rdname get_parameter_dims
#' @export
get_parameter_dims.dynamitefit <- function(x, ...) {
  pars_text <- get_code(x, blocks = "parameters")
  pars_text <- strsplit(pars_text, split = "\n")[[1L]]
  pars_text <- pars_text[grepl(";", pars_text)]
  par_regex <- regexec(
    pattern = "^.+\\s([^\\s]+);.*$",
    text = pars_text,
    perl = TRUE
  )
  par_matches <- regmatches(pars_text, par_regex)
  par_names <- vapply(par_matches, "[[", character(1L), 2L)
  dim_regex <- regexec(
    pattern = "^[^\\[]+\\[([^\\]]+)\\].+",
    text = pars_text,
    perl = TRUE
  )
  dim_matches <- regmatches(pars_text, dim_regex)
  dim_names <- lapply(
    dim_matches,
    function(y) {
      if (length(y) > 0L) {
        paste0("c(", y[2L], ")")
      } else {
        "1"
      }
    }
  )
  e <- list2env(get_data(x, ...))
  stats::setNames(
    lapply(
      dim_names,
      function(y) {
        eval(str2lang(y), envir = e)
      }
    ),
    par_names
  )
}

#' Internal Parameter Block Variable Name Extraction
#'
#' @param x \[`character(1L)`]\cr The Stan model code string of the
#'   "Parameters" block.
#' @noRd
get_parameters <- function(x) {
  x <- strsplit(x, split = "\n")[[1L]]
  x <- x[grepl(";", x)]
  par_regex <- regexec(
    pattern = "^.+\\s([^\\s]+);.*$",
    text = x,
    perl = TRUE
  )
  par_matches <- regmatches(x, par_regex)
  vapply(par_matches, "[[", character(1L), 2L)
}

#' Get Parameter Types of the \pkg{dynamite} Model
#'
#' Extracts all parameter types of used in the `dynamitefit` object. See
#' [as.data.frame.dynamitefit()] for explanations of different types.
#'
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param ... Ignored.
#' @return A `character` vector with all parameter types of the input model.
#' @export
#' @family output
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' get_parameter_types(multichannel_example_fit)
#'
get_parameter_types <- function(x, ...) {
  UseMethod("get_parameter_types", x)
}

#' @rdname get_parameter_types
#' @export
get_parameter_types.dynamitefit <- function(x, ...) {
  d <- as.data.table(x, types = all_types)
  unique(d$type)
}

#' Get Parameter Names of the \pkg{dynamite} Model
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
#'   certain type. See [get_parameter_types()].
#' @param ... Ignored.
#' @return A `character` vector with parameter names of the input model.
#' @export
#' @family output
#' @examples
#' data.table::setDTthreads(1) # For CRAN
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

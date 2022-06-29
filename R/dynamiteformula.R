#' Model formula for \pkg{dynamite}
#'
#' TODO description and details of all features
#' TODO document families here as well
#' TODO document lag-conversion to data
#' TODO explain fixed time points
#'
#' @param formula \[`formula`]\cr An R formula describing the model.
#' @param family \[`character(1)`]\cr The family name.
#' @param random_intercept \[`logical(1)`]\cr If `TRUE`, adds
#'   individual-level intercepts to the channel.
#' @export
#' @examples
#' obs(y ~ -1 + varying(~x), family = gaussian()) +
#'   lags(type = "varying") + splines(df = 20)
#'
#' @srrstats {G2.3b} *Either: use `tolower()` or equivalent to ensure input of character parameters is not case dependent; or explicitly document that parameters are strictly case-sensitive.*
#' @srrstats {RE1.0} *Regression Software should enable models to be specified via a formula interface, unless reasons for not doing so are explicitly documented.*
#' @srrstats {RE1.1} *Regression Software should document how formula interfaces are converted to matrix representations of input data.*
#' @srrstats {RE1.4} *Regression Software should document any assumptions made with regard to input data; for example distributional assumptions, or assumptions that predictor data have mean values of zero. Implications of violations of these assumptions should be both documented and tested.*
dynamiteformula <- function(formula, family, random_intercept = FALSE) {
  stopifnot_(
    is.formula(formula),
    "Argument {.arg formula} must be a {.cls formula} object."
  )
  stopifnot_(
    checkmate::test_string(x = family, na.ok = FALSE),
    "Argument {.arg family} must be a single {.cls character} string."
  )
  stopifnot_(
    checkmate::test_flag(x = random_intercept),
    "Argument {.arg random_intercept} must be a single {.cls logical} value."
  )
  family <- tolower(family)
  if (is_supported(family)) {
    family <- do.call(paste0(family, "_"), args = list())
  } else {
    stop_("Family {.val {family}} is not supported.")
  }
  stopifnot_(
    !has_as_is(deparse1(formula)),
    "{.code I(.)} is not supported by {.fun dynamiteformula}."
  )
  x <- dynamiteformula_(formula, formula, family, random_intercept)
  structure(
    list(
      dynamitechannel(
        formula = x$formula,
        original = formula,
        family = x$family,
        response = x$response,
        fixed = x$fixed,
        varying = x$varying,
        specials = x$specials,
        has_fixed_intercept = x$has_fixed_intercept,
        has_varying_intercept = x$has_varying_intercept,
        has_random_intercept = x$has_random_intercept
      )
    ),
    class = "dynamiteformula"
  )
}

#' @describeIn dynamiteformula Internal Version of `dynamiteformula`
#' @noRd
dynamiteformula_ <- function(formula, original, family,
                             random_intercept = FALSE) {
  if (is_deterministic(family)) {
    out <- formula_past(formula)
    resp_parsed <- formula_response(deparse1(formula_lhs(formula)))
    out$specials$resp_type <- resp_parsed$type
    if (!is.null(out$specials$past)) {
      out$specials$past <- do.call(
        what = paste0("as.", resp_parsed$type),
        args = list(out$specials$past)
      )
    }
    out$response <- resp_parsed$resp
  } else {
    out <- formula_specials(formula)
    out$response <- deparse1(formula_lhs(formula))
  }
  out$family <- family
  stopifnot_(
    !(random_intercept && is_categorical(family)),
    "Random intercepts are not yet supported for the categorical family."
  )
  out$original <- original
  out$has_random_intercept <- random_intercept
  out
}

#' Create a Channel for a `dynamiteformula` Directly
#'
#' @inheritParams dynamiteformula
#' @param response \[`character(1)`]\cr Name of the response.
#' @param fixed \[`integer()`]\cr Time-invariant covariate indices.
#' @param varying \[`integer()`]\cr Time-varying covariate indices.
#' @param has_fixed_intercept \[`logical(1)`]\cr Does the channel contain fixed
#'   intercept?
#' @param has_varying_intercept \[`logical(1)`]\cr Does the channel contain
#'   varying intercept?
#' @param has_rand_intercept \[`logical(1)`]\cr Does the channel contain random
#'   individual-level intercept term?
#'@noRd
dynamitechannel <- function(formula, original = NULL, family, response,
                            fixed = integer(0L), varying = integer(0L),
                            specials = list(),
                            has_fixed_intercept = FALSE,
                            has_varying_intercept = FALSE,
                            has_random_intercept = FALSE) {
  list(
    formula = formula,
    original = original,
    family = family,
    response = response,
    fixed = fixed,
    varying = varying,
    specials = specials,
    has_fixed_intercept = has_fixed_intercept,
    has_varying_intercept = has_varying_intercept,
    has_random_intercept = has_random_intercept
  )
}

#' @rdname dynamiteformula
#' @export
obs <- dynamiteformula

#' Is the Argument a `dynamiteformula` Object
#'
#' @param x An \R object.
#' @noRd
is.dynamiteformula <- function(x) {
  inherits(x, "dynamiteformula")
}

#' @describeIn dynamiteformula Prepare a Deterministic Auxiliary Channel
#' @export
aux <- function(formula) {
  dynamiteformula(formula, family = "deterministic")
}

#' Join Two `dynamiteformula` Objects
#'
#' @param e1 An \R object.
#' @param e2 An \R object.
#' @export
#' @examples
#' obs(y ~ x, family = "gaussian") + obs(z ~ w, family = "exponential")
`+.dynamiteformula` <- function(e1, e2) {
  if (is.dynamiteformula(e1)) {
    out <- add_dynamiteformula(e1, e2)
  } else {
    stop_("Method {.fun +.dynamiteformula} is not supported for
          {.cls {class(e1)}} objects.")
  }
  out
}

#' Is the Argument a `formula` Object.
#'
#' @param x An \R object
#' @noRd
is.formula <- function(x) {
  inherits(x, "formula")
}

#' Get All Response Variables of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_responses <- function(x) {
  vapply(x, function(y) y$response, character(1L))
}

#' Get The RHS of All Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_predictors <- function(x) {
  vapply(x, function(y) deparse1(formula_rhs(y$formula)), character(1L))

}

#' Get Terms of All Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_terms <- function(x) {
  lapply(x, function(y) {
    if (is_deterministic(y$family)) {
      character(0L)
    } else {
      attr(terms(y$formula), "term.labels")
    }
  })
}

#' Get All Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_formulas <- function(x) {
  lapply(x, "[[", "formula")
}

#' Get All Original Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_originals <- function(x) {
  lapply(x, "[[", "original")
}

#' Get All Family Objects of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_families <- function(x) {
  lapply(x, "[[", "family")
}

#' Get a Quoted Expression of Deterministic Channel Definitions
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_quoted <- function(x) {
  resp <- get_responses(x)
  if (length(resp) > 0L) {
    expr <- lapply(x, function(x) deparse1(formula_rhs(x$formula)))
    quote_str <- paste0(
      "`:=`(",
      paste0(resp, " = ", expr, collapse = ","),
      ")"
    )
    str2lang(quote_str)
  } else {
    NULL
  }
}

#' Get Indices of Deterministic Channels in a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
which_deterministic <- function(x) {
  which(vapply(x, function(y) is_deterministic(y$family), logical(1L)))
}

#' Get Indices of Stochastic Channels in a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object
#' @noRd
which_stochastic <- function(x) {
  which(vapply(x, function(y) !is_deterministic(y$family), logical(1L)))
}

#' Get Channels with Past Value Definitions of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object
#' @noRd
has_past <- function(x) {
  vapply(x, function(y) length(y$specials$past) > 0L, logical(1L))
}

#' Internal `+.dynamiteformula` For Model Construction
#'
#' @param e1 An \R object.
#' @param e2 An \R object.
#' @noRd
add_dynamiteformula <- function(e1, e2) {
  if (is.dynamiteformula(e2)) {
    out <- join_dynamiteformulas(e1, e2)
  } else if (is.lags(e2)) {
    out <- set_lags(e1, e2)
  } else if (is.splines(e2)) {
    out <- set_splines(e1, e2)
  } else {
    stop_(
      "Unable to add an object of class {.cls {class(e2)}}
      to an object of class {.cls dynamiteformula}."
    )
  }
  out
}

#' Join Two Model Definitions and Verify Compatibility
#'
#' @param e1 A `dynamiteformula` object.
#' @param e2 A `dynamiteformula` object.
#' @noRd
join_dynamiteformulas <- function(e1, e2) {
  out <- c(e1, e2)
  resp_list <- list(
    get_responses(e1),
    get_responses(e2)
  )
  resp_all <- unlist(resp_list)
  resp_duped <- duplicated(resp_all)
  stopifnot_(
    !any(resp_duped),
    "Multiple definitions for response variable{?s}
     {.var {resp_all[resp_duped]}}."
  )
  stopifnot_(
    is.null(attr(e1, "lags")) || is.null(attr(e2, "lags")),
    "Both dynamiteformulas contain a lags definition."
  )
  stopifnot_(
    is.null(attr(e1, "splines")) || is.null(attr(e2, "splines")),
    "Both dynamiteformulas contain a splines definition."
  )
  rhs_list <- list(
    lapply(get_terms(e1), extract_nonlags),
    lapply(get_terms(e2), extract_nonlags)
  )
  stoch_list <- list(
    which_stochastic(e1),
    which_stochastic(e2)
  )
  for (i in 1L:2L) {
    resp_a <- resp_list[[i]]
    resp_b <- resp_list[[3L-i]][stoch_list[[3L-i]]]
    rhs <- rhs_list[[3L-i]][stoch_list[[3L-i]]]
    if (length(rhs) > 0L) {
      for (j in seq_along(resp_a)) {
        simul_lhs <- resp_a[j]
        simul <- vapply(rhs, function(x) simul_lhs %in% x, logical(1L))
        stopifnot_(
          !any(simul),
          c(
            "Simultaneous regression is not supported:",
            `x` = "Response variable {.var {simul_lhs}} appears in
                  the formula of {.var {resp_b[which(simul)[1L]]}}."
          )
        )
      }
    }
  }
  attributes(out) <- c(attributes(e1), attributes(e2))
  class(out) <- "dynamiteformula"
  out
}

#' Set Lags Definitions for All Channels in a `dynamiteformula` Object
#'
#' @param e1 A `dynamiteformula` object,
#' @param e2 A `lags` object,
#' @noRd
set_lags <- function(e1, e2) {
  stopifnot_(
    is.null(attr(e1, "lags")),
    "Multiple definitions for lags."
  )
  attr(e1, "lags") <- e2
  e1
}

#' Set the Regression Coefficient Splines of the Model
#'
#' @param e1 A `dynamiteformula` object.
#' @param e2 A `splines` object.
#' @noRd
set_splines <- function(e1, e2) {
  stopifnot_(
    is.null(attr(e1, "splines")) || attr(e2, "override"),
    "Multiple definitions for splines."
  )
  attr(e1, "splines") <- e2
  e1
}

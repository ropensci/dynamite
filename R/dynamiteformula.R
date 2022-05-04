#' Model formula for \pkg{dynamite}
#'
#' TODO description and details of all features
#'
#' @param formula \[`formula`]\cr An R formula describing the model.
#' @param family \[`call`]\cr A call to a family function, e.g., `gaussian()`.
#' @param ... TODO
#'
#' @export
dynamiteformula <- function(formula, family, ...) {
  if (!is.formula(formula)) {
    stop_("Argument 'formula' is not a formula.")
  }
  family_call <- is_valid_family_call(substitute(family))
  if (!family_call$supported) {
    stop_("Unsupported family object")
  }
  x <- formula_specials(formula)
  structure(
    list(
      list(
        formula = x$formula,
        family = eval(family_call$call),
        response = formula_lhs(x$formula),
        predictors = formula_rhs(x$formula),
        fixed = x$fixed,
        varying = x$varying,
        specials = x$specials
      )
    ),
    class = "dynamiteformula"
  )
}

#' @rdname dynamiteformula
#' @export
obs <- dynamiteformula

#' Checks if argument is a `dynamiteformula` object
#'
#' @param x An \R object
#'
#' @export
is.dynamiteformula <- function(x) {
  inherits(x, "dynamiteformula")
}

#' Join two dynamiteformulas
#'
#' @param e1 An \R object
#' @param e2 An \R object
#'
#' @export
`+.dynamiteformula` <- function(e1, e2) {
  if (is.dynamiteformula(e1)) {
    out <- add_dynamiteformula(e1, e2)
  } else {
    stop_("Method '+.dynamiteformula is not supported for ",
          class(e1), " objects.")
  }
  out
}

#' Checks if argument is a formula
#'
#' @param x An \R object
#'
#' @noRd
is.formula <- function(x) {
  inherits(x, "formula")
}

#' Get all response variables of a dynamiteformula object
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_resp <- function(x) {
  sapply(x, "[[", "response")
}

#' Get all predictor variables of a dynamiteformula object
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_pred <- function(x) {
  lapply(x, "[[", "predictors")
}

#' Get all formulas of a dynamiteformula object
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_form <- function(x) {
  lapply(x, "[[", "formula")
}

#' Check whether a dynamiteformula contains an intercept
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
has_intercept <- function(x) {
  attr(terms(x$formula), "intercept") == 1
}

#' Internal +.dynamiteformula for model constructions
#'
#' @param e1 An \R object
#' @param e2 An \R object
#'
#' @noRd
add_dynamiteformula <- function(e1, e2) {
  if (is.dynamiteformula(e2)) {
    out <- join_dynamiteformulas(e1, e2)
  } else if (is.lags(e2)) {
    out <- set_lags(e1, e2)
  } else if (is.splines(e2)) {
    out <- set_splines(e1, e2)
  } else if (is.modeldata(e2)) {
    out <- set_modeldata(e1, e2)
  } else {
    stop_(
      "Unable to add an object of class ", class(e2),
      " to an object of class 'dynamiteformula'"
    )
  }
  out
}

#' Join two model definitions and verify compatibility
#'
#'
#' @param e1 A `dynamiteformula` object
#' @param e2 A `dynamiteformula` object
#'
#' @noRd
join_dynamiteformulas <- function(e1, e2) {
  out <- c(e1, e2)
  resp_all <- get_resp(out)
  duped <- duplicated(resp_all)
  if (any(duped)) {
    stop_("Multiple definitions for response variables: ", resp_all[duped])
  }
  if (!is.null(attr(e1, "lag_all")) && !is.null(attr(e2, "lag_all"))) {
    stop_("Multiple definitions for lags")
  }
  if (!is.null(attr(e1, "splines")) && !is.null(attr(e2, "splines"))) {
    stop_("Multiple definitions for splines")
  }
  # if (!is.null(e1$data) && !is.null(e2$data)) {
  #     stop_("Multiple definitions for data")
  # }
  attributes(out) <- c(attributes(e1), attributes(e2))
  class(out) <- "dynamiteformula"
  out
}

#' Set lag definitions for all channels
#'
#' @param e1 A `dynamiteformula` object
#' @param e2 A `lags` object
#'
#' @noRd
set_lags <- function(e1, e2) {
  if (!is.null(attr(e1, "lags"))) {
    stop_("Multiple definitions for lags")
  }
  attr(e1, "lags") <- e2
  e1
}

#' Set the regression coefficient splines of the model
#'
#' @param e1 A `dynamiteformula` object
#' @param e2 A `splines` object
#'
#' @noRd
set_splines <- function(e1, e2) {
  if (!is.null(attr(e1, "splines")) && !attr(e2, "override")) {
    stop_("Multiple definitions for splines")
  }
  attr(e1, "splines") <- e2
  e1
}

#' Model formula for \pkg{dynamite}
#'
#' TODO description and details of all features
#'
#' @param formula \[`formula`]\cr An R formula describing the model.
#' @param family \[`call`, `character(1)`]\cr
#'   A call to a family function, e.g., `gaussian()` or the family
#'   name as character, e.g., `"gaussian"`.
#'
#' @export
dynamiteformula <- function(formula, family) {
  if (!is.formula(formula)) {
    stop_("Argument 'formula' is not a formula object")
  }
  family <- substitute(family)
  if (is.character(family)) {
    if (is_supported(family[1])) {
      family <- do.call(paste0(family, "_"), args = list())
    } else {
      stop_("Family '", family[1], "' is not supported")
    }
  } else {
    family_call <- is_valid_family_call(family)
    if (!family_call$supported) {
      stop_("Unsupported family call '", family_call$call_str, "()'")
    } else {
      family <- eval(family_call$call)
    }
  }
  if (has_as_is(deparse(formula))) {
    stop_("The use of I(.) is not supported by dynamiteformula")
  }
  x <- dynamiteformula_(formula, family)
  structure(
    list(
      dynamitechannel(
        formula = x$formula,
        family = x$family,
        response = x$response,
        fixed = x$fixed,
        varying = x$varying,
        specials = x$specials
      )
    ),
    class = "dynamiteformula"
  )
}

#' @describeIn dynamiteformula Internal version of dynamiteformula
#'
#' @noRd
dynamiteformula_ <- function(formula, family) {
  if (is_deterministic(family)) {
    out <- formula_past(formula)
  } else {
    out <- formula_specials(formula)
  }
  out$family <- family
  out$response <- as.character(formula_lhs(formula))
  out
}

#' Create a channel for a dynamiteformula directly
#'
#' @param formula See [`dynamiteformula()`]
#' @param family See [`dynamiteformula()`]
#' @param response \[`character(1)`] Name of the response
#' @param fixed \[`integer()`] Time-invariant covariate indices
#' @param varying \[`integer()`] Time-varying covariate indices
#' @param specials See [`dynamiteformula()`]
#'
#'@noRd
dynamitechannel <- function(formula, family, response,
                            fixed = integer(0), varying = integer(0),
                            specials = list()) {
  list(
    formula = formula,
    family = family,
    response = response,
    fixed = fixed,
    varying = varying,
    specials = specials
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

#' @describeIn dynamiteformula Prepare a deterministic auxiliary channel
#' @export
aux <- function(formula) {
  dynamiteformula(formula, family = "deterministic")
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
    stop_("Method '+.dynamiteformula' is not supported for '",
          class(e1), "' objects")
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
get_responses <- function(x) {
  sapply(x, "[[", "response")
}

#' Get the RHS of all formulas of a dynamiteformula object as a character vector
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_predictors <- function(x) {
  sapply(x, function(y) deparse(formula_rhs(y$formula)))
}

#' Get all formulas of a dynamiteformula object
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_formulas <- function(x) {
  lapply(x, "[[", "formula")
}

#' Get indices of deterministic channels in dynamiteformula
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
which_deterministic <- function(x) {
  which(sapply(x, function(y) is_deterministic(y$family)))
}

#' Get indices of stochastic channels in dynamiteformula
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
which_stochastic <- function(x) {
  which(sapply(x, function(y) !is_deterministic(y$family)))
}

#' Get channels with past value definitions
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
has_past <- function(x) {
  sapply(x, function(y) length(y$specials$past) > 0)
}

#' Get ranks of channels for evaluation order of precedence
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_ranks <- function(x) {
  sapply(x, function(y) y$specials$rank)
}

# TODO can delete? this is not used
# #' Check whether a dynamiteformula contains an intercept
# #'
# #' @param x A `dynamiteformula` object
# #'
# #' @noRd
# has_intercept <- function(x) {
#   attr(terms(x$formula), "intercept") == 1
# }

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
  } else {
    stop_(
      "Unable to add an object of class '", class(e2),
      "' to an object of class 'dynamiteformula'"
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
  resp_list <- list()
  resp_list[[1]] <- get_responses(e1)
  resp_list[[2]] <- get_responses(e2)
  resp_all <- unlist(resp_list)
  resp_duped <- duplicated(resp_all)
  if (any(resp_duped)) {
    stop_("Multiple definitions for response variable(s): ",
          cs(resp_all[resp_duped]))
  }
  if (!is.null(attr(e1, "lags")) && !is.null(attr(e2, "lags"))) {
    stop_("Multiple definitions for lags")
  }
  if (!is.null(attr(e1, "splines")) && !is.null(attr(e2, "splines"))) {
    stop_("Multiple definitions for splines")
  }
  rhs_list <- list()
  rhs_list[[1]] <- get_predictors(e1)
  rhs_list[[2]] <- get_predictors(e2)
  for (i in 1:2) {
    resp_a <- resp_list[[i]]
    resp_b <- resp_list[[3-i]]
    rhs <- rhs_list[[3-i]]
    simul_resp <- which(resp_a %in% rhs)
    if (length(simul_resp) > 0) {
      simul_rhs <- which(rhs %in% resp_a)
      stop_("Simultaneous regression is not supported, response variables '",
            cs(resp_a[simul_resp]), "' appear in the formulas of '",
            cs(resp_b[simul_rhs]), "'")
    }
  }
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

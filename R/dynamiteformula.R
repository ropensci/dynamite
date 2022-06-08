#' Model formula for \pkg{dynamite}
#'
#' TODO description and details of all features
#'
#' @param formula \[`formula`]\cr An R formula describing the model.
#' @param family \[`call`, `character(1)`]\cr
#'   A call to a family function, e.g., `gaussian()` or the family
#'   name as character, e.g., `"gaussian"`.
#' @param random_intercept \[`logical(1)`]\cr If `TRUE`, adds
#'   individual-level intercepts to the channel.
#' @export
dynamiteformula <- function(formula, family, random_intercept = FALSE) {
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
  if (has_as_is(deparse1(formula))) {
    stop_("The use of I(.) is not supported by dynamiteformula")
  }
  x <- dynamiteformula_(formula, family, random_intercept)
  structure(
    list(
      dynamitechannel(
        formula = x$formula,
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

#' @describeIn dynamiteformula Internal version of dynamiteformula
#'
#' @noRd
dynamiteformula_ <- function(formula, family, random_intercept = FALSE) {
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
  out$has_random_intercept <- random_intercept
  if (random_intercept && is_categorical(family)) {
    stop_("Random intercepts are not yet supported for the categorical family")
  }
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
#' @param has_fixed_intercept \[`logical(1)`] Does the channel contain fixed
#'   intercept?
#' @param has_varying_intercept \[`logical(1)`] Does the channel contain
#'   varying intercept?
#' @param has_rand_intercept \[`logical(1)`] Does the channel contain random
#'   individual-level intercept term?
#'@noRd
dynamitechannel <- function(formula, family, response,
                            fixed = integer(0), varying = integer(0),
                            specials = list(),
                            has_fixed_intercept = FALSE,
                            has_varying_intercept = FALSE,
                            has_random_intercept = FALSE) {
  list(
    formula = formula,
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
  unlist(sapply(x, "[[", "response"))
}

#' Get the RHS of all formulas of a dynamiteformula object as a character vector
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_predictors <- function(x) {
  sapply(x, function(y) deparse1(formula_rhs(y$formula)))

}

#' Get terms of all formulas of a dynamiteformula
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_terms <- function(x) {
  lapply(x, function(y) {
    if (is_deterministic(y$family)) {
      character(0)
    } else {
      attr(terms(y$formula), "term.labels")
    }
  })
}

#' Get all formulas of a dynamiteformula object
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_formulas <- function(x) {
  lapply(x, "[[", "formula")
}

#' Get all family objects of a dynamiteformula object
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_families <- function(x) {
  lapply(x, "[[", "family")
}

#' Get a quoted expression of deterministic channel definitions
#'
#' @param x A `dynamiteformula` object
#'
#' @noRd
get_quoted <- function(x) {
  resp <- get_responses(x)
  if (length(resp) > 0) {
    expr <- lapply(x, function(x) deparse1(formula_rhs(x$formula)))
    quote_str <- paste0("`:=`(", paste0(resp, " = ", expr, collapse = ","), ")")
    str2lang(quote_str)
  } else {
    NULL
  }
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
#' @param e1 A `dynamiteformula` object
#' @param e2 A `dynamiteformula` object
#'
#' @noRd
join_dynamiteformulas <- function(e1, e2) {
  out <- c(e1, e2)
  resp_list <- list(
    get_responses(e1),
    get_responses(e2)
  )
  resp_all <- unlist(resp_list)
  resp_duped <- duplicated(resp_all)
  if (any(resp_duped)) {
    stop_("Multiple definitions for response variable(s): ",
          cs(resp_all[resp_duped]))
  }
  if (!is.null(attr(e1, "lags")) && !is.null(attr(e2, "lags"))) {
    stop_("Both dynamiteformulas contain a lags definition")
  }
  if (!is.null(attr(e1, "splines")) && !is.null(attr(e2, "splines"))) {
    stop_("Both dynamiteformulas contain a splines definition")
  }
  rhs_list <- list(
    lapply(get_terms(e1), extract_nonlags),
    lapply(get_terms(e2), extract_nonlags)
  )
  stoch_list <- list(
    which_stochastic(e1),
    which_stochastic(e2)
  )
  for (i in 1:2) {
    resp_a <- resp_list[[i]]
    resp_b <- resp_list[[3-i]][stoch_list[[3-i]]]
    rhs <- rhs_list[[3-i]][stoch_list[[3-i]]]
    if (length(rhs) > 0) {
      for (j in seq_along(resp_a)) {
        simul_lhs <- resp_a[j]
        simul <- sapply(rhs, function(x) simul_lhs %in% x)
        if (any(simul)) {
          simul_rhs <- resp_b[which(simul)[1]]
          stop_("Simultaneous regression is not supported, response variable '",
                simul_lhs, "' appears in the formula of '",
                simul_rhs, "'")
        }
      }
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

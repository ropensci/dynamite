#' Get and Separate All Specials of a Formula Object
#'
#' @param x A `formula` object.
#' @noRd
formula_specials <- function(x) {
  out <- list(formula = NULL, specials = NULL, coefs = NULL)
  xt <- terms(x, specials = formula_special_funs)
  xt_specials <- attr(xt, "specials")[formula_special_funs]
  xt_variables <- attr(xt, "variables")
  xt_terms <- attr(xt, "term.labels")
  for (y in formula_special_funs) {
    out$specials[[y]] <- onlyif(
      !is.null(xt_specials[[y]]),
      xt_variables[[xt_specials[[y]] + 1]][[2]]
    )
  }
  special_vars <- unlist(xt_specials)
  special_vars <- ifelse_(
    !is.null(xt_specials$offset),
    special_vars[!names(special_vars) %in% "offset"],
    special_vars
  )
  x <- drop_terms(
    termobj = xt,
    dropx = get_special_term_indices(
      special_vars,
      xt_variables,
      xt_terms
    )
  )
  xt <- terms(x, specials = c("fixed", "varying", "random"))
  xt_specials <- attr(xt, "specials")[c("fixed", "varying", "random")]
  xt_variables <- attr(xt, "variables")
  xt_terms <- attr(xt, "term.labels")
  fixed_icpt <- 0L
  special_vars <- unlist(xt_specials)
  fixed_terms <- character(0L)
  if (!is.null(xt_specials[["fixed"]])) {
    # eval to ensure fixed_form is a formula
    fixed_form <- eval(xt_variables[[xt_specials[["fixed"]] + 1]][[2]])
    fixed_terms <- formula_terms(fixed_form)
    fixed_icpt <- attr(terms(fixed_form), "intercept")
  }
  varying_terms <- character(0L)
  varying_icpt <- 0L
  if (!is.null(xt_specials[["varying"]])) {
    # eval to ensure varying_form is a formula
    varying_form <- eval(xt_variables[[xt_specials[["varying"]] + 1]][[2]])
    varying_terms <- formula_terms(varying_form)
    varying_icpt <- attr(terms(varying_form), "intercept")
  }
  random_terms <- character(0L)
  random_icpt <- 0L
  if (!is.null(xt_specials[["random"]])) {
    # eval to ensure random_form is a formula
    random_form <- eval(xt_variables[[xt_specials[["random"]] + 1]][[2]])
    random_terms <- formula_terms(random_form)
    random_icpt <- attr(terms(random_form), "intercept")
  }
  x <- drop_terms(
    termobj = xt,
    dropx = get_special_term_indices(
      special_vars,
      xt_variables,
      xt_terms
    )
  )
  form_terms <- formula_terms(x)
  fixed_terms <- union(form_terms, fixed_terms)
  fixed_icpt <- attr(xt, "intercept") || fixed_icpt
  full_terms <- union(fixed_terms, union(varying_terms, random_terms))
  any_icpt <- fixed_icpt || varying_icpt || random_icpt
  if (fixed_icpt && varying_icpt) {
    warning_(c(
      "Both time-independent and time-varying intercept specified:",
      `i` = "Defaulting to time-varying intercept."
    ))
    fixed_icpt <- FALSE
  }
  if (length(full_terms) > 0L) {
    x <- reformulate(
      termlabels = full_terms,
      response = xt_variables[[2]],
      intercept = 1 * any_icpt
    )
  } else {
    y <- as.character(xt_variables[[2]])
    stopifnot_(
      any_icpt,
      c(
        "Invalid formula for response variable {.var {y}}:",
        `x` = "There are no predictors nor an intercept term."
      )
    )
    x <- as.formula(paste0(y, "~ 1"))
  }
  xt <- formula_terms(x)
  out$formula <- x
  out$fixed <- which(xt %in% fixed_terms)
  out$has_fixed_intercept <- as.logical(fixed_icpt)
  out$varying <- which(xt %in% varying_terms)
  out$has_varying_intercept <- as.logical(varying_icpt)
  out$random <- which(xt %in% random_terms)
  out$has_random_intercept <- as.logical(random_icpt)
  out
}


#' Process Formulas for Deterministic Channels and Get Past Value Definitions
#'
#' @param x A `formula` object.
#' @noRd
formula_past <- function(formula) {
  formula_str <- deparse1(formula)
  rhs <- formula_rhs(formula)
  past_def <- NULL
  past_type <- NULL
  if (length(rhs) > 1L && identical(deparse1(rhs[[1L]]), "|")) {
    past_type <- deparse1(rhs[[3L]][[1L]])
    stopifnot_(
      past_type %in% c("init", "past"),
      c(
        "Invalid formula of a deterministic channel:",
        `x` = "Unsupported definition {.var {past_type}}
               on the right-hand side of {.fun |}."
      )
    )
    stopifnot_(
      identical(length(rhs[[3L]]), 2L),
      "Invalid number of arguments supplied to {.fun {past_type}}
       in {.var {formula_str}}."
    )
    past_def <- rhs[[3L]][[2L]]
    formula[[3]] <- rhs[[2L]]
  }
  stopifnot_(
    !grepl("fixed\\(.+\\)", formula_str, perl = TRUE),
    c(
      "The use of {.fun fixed} is not meaningful for deterministic channels:",
      `x` = "Time-invariant definition was found in {.var {formula_str}}."
    )
  )
  stopifnot_(
    !grepl("varying\\(.+\\)", formula_str, perl = TRUE),
    c(
      "The use of {.fun varying} is not meaningful for deterministic channels:",
      `x` = "Time-varying definition was found in {.var {formula_str}}."
    )
  )
  stopifnot_(
    !grepl("random\\(.+\\)", formula_str, perl = TRUE),
    c(
      "The use of {.fun random} is not meaningful for deterministic channels:",
      `x` = "Random effect definition was found in {.var {formula_str}}."
    )
  )
  list(
    formula = formula,
    specials = list(
      past_type = past_type,
      past = past_def,
      rank = Inf
    ),
    fixed = integer(0L),
    varying = integer(0L),
    random = integer(0L)
  )
}

#' Process Response Variables for Deterministic Channels
#'
#' @param y \[`character(1)`] The response variable name.
#' @noRd
deterministic_response <- function(y) {
  if (grepl("factor\\(.*\\)", y, perl = TRUE)) {
    list(
      type = "factor",
      resp = gsub("factor\\((.*)\\)", "\\1", y, perl = TRUE)
    )
  } else if (grepl("numeric\\(.*\\)", y, perl = TRUE)) {
    list(
      type = "numeric",
      resp = gsub("numeric\\((.*)\\)", "\\1", y, perl = TRUE)
    )
  } else if (grepl("integer\\(.*\\)", y, perl = TRUE)) {
    list(
      type = "integer",
      resp = gsub("integer\\((.*)\\)", "\\1", y, perl = TRUE)
    )
  } else if (grepl("logical\\(.*\\)", y, perl = TRUE)) {
    list(
      type = "logical",
      resp = gsub("logical\\((.*)\\)", "\\1", y, perl = TRUE)
    )
  } else {
    warning_(c(
      "No type specified for deterministic channel {.var {y}}:",
      `i` = "Assuming type is {.cls numeric}."
    ))
    list(type = "numeric", resp = y)
  }
}

#' Computes All Specials Defined in a Formula in the Context of the Data
#'
#' @param dformula \[`dformula`]\cr The model formula object.
#' @param data \[`data.table`]\cr Data containing the variables used for
#'   the special definitions in the formula.
#' @noRd
evaluate_specials <- function(dformula, data) {
  lapply(seq_along(dformula), function(i) {
    if (length(dformula[[i]]$specials) > 0L) {
      out <- list()
      for (spec in formula_special_funs) {
        spec_formula <- dformula[[i]]$specials[[spec]]
        if (!is.null(spec_formula)) {
          out[[spec]] <- data[, eval(spec_formula)]
        }
      }
      out
    } else {
      NULL
    }
  })
}

#' Retrieve the Corresponding Term of a Special Variable in a Formula
#'
#' @param special \[`integer()`]\cr A vector of special variable indices.
#' @param vars \[`language`]\cr
#'   The `"variables"` attribute of a `terms` object.
#' @param term_labels \[`character()`]\cr
#'   The `"term.labels"` attribute of a `terms` object.
#'
#' @noRd
get_special_term_indices <- function(special, vars, term_labels) {
  out <- integer(length(special))
  for (i in seq_along(special)) {
    v <- deparse1(vars[[special[i] + 1L]])
    out[i] <- which(term_labels == v)
  }
  out
}

# Supported specials
formula_special_funs <- c(
  "offset",
  "trials"
)

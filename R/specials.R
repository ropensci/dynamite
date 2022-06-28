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
  x <- ifelse_(
    length(special_vars) > 0,
    formula(
      drop.terms(
        xt,
        dropx = get_special_term_indices(
          special_vars,
          xt_variables,
          xt_terms
        ),
        keep.response = TRUE
      )
    ),
    x
  )
  xt <- terms(x, specials = c("fixed", "varying"))
  xt_specials <- attr(xt, "specials")[c("fixed", "varying")]
  xt_variables <- attr(xt, "variables")
  xt_terms <- attr(xt, "term.labels")
  fixed_icpt <- 0
  special_vars <- unlist(xt_specials)
  fixed_terms <- character(0)
  if (!is.null(xt_specials[["fixed"]])) {
    # eval to ensure fixed_form is a formula
    fixed_form <- eval(xt_variables[[xt_specials[["fixed"]] + 1]][[2]])
    fixed_terms <- formula_terms(fixed_form)
    fixed_icpt <- attr(terms(fixed_form), "intercept")
  }
  varying_terms <- character(0)
  varying_icpt <- 0
  if (!is.null(xt_specials[["varying"]])) {
    # eval to ensure varying_form is a formula
    varying_form <- eval(xt_variables[[xt_specials[["varying"]] + 1]][[2]])
    varying_terms <- formula_terms(varying_form)
    varying_icpt <- attr(terms(varying_form), "intercept")
  }
  if (!is.null(special_vars)) {
    if (length(xt_terms) > length(special_vars)) {
      x <- formula(
        drop.terms(
          xt,
          dropx = get_special_term_indices(
            special_vars,
            xt_variables,
            xt_terms
          ),
          keep.response = TRUE
        )
      )
      form_terms <- formula_terms(x)
    } else {
      form_terms <- character(0)
    }
  } else {
    form_terms <- formula_terms(x)
  }
  fixed_terms <- union(form_terms, fixed_terms)
  fixed_icpt <- attr(xt, "intercept") || fixed_icpt
  full_terms <- union(fixed_terms, varying_terms)
  any_icpt <- fixed_icpt || varying_icpt
  if (fixed_icpt && varying_icpt) {
    warning_(c(
      "Both time-independent and time-varying intercept specified:",
      `i` = "Defaulting to time-varying intercept."
    ))
    fixed_icpt <- FALSE
  }
  if (length(full_terms) > 0) {
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
  out
}

#' Process Formulas for Deterministic Channels and Get Past Value Definitions
#'
#' @param x A `formula` object.
#' @noRd
formula_past <- function(formula) {
  formula_str <- deparse1(formula)
  form_comp <- regexpr(
    pattern = "^(?<resp>[^~]+) ~ (?<def>.+?)(?: \\+ past\\((?<past>.+)\\)){0,1}$",
    text = formula_str,
    perl = TRUE
  )
  start <- attr(form_comp, "capture.start")
  end <- start + attr(form_comp, "capture.length") - 1
  form_resp <- substr(formula_str, start[1], end[1])
  form_def <- substr(formula_str, start[2], end[2])
  if (grepl("past\\(", form_def, perl = TRUE)) {
    stop_("Past values term must be the last term of the formula.")
  }
  form_past <- substr(formula_str, start[3], end[3])
  form_both <- c(form_def, form_past)
  if (any(grepl("fixed\\(.+\\)", form_both, perl = TRUE))) {
    warning_(
      "fixed() definitions of a determinstic channel
       {.var {deparse1(formula_lhs(formula))}} will be ignored."
    )
  }
  if (any(grepl("varying\\(.+\\)", form_both, perl = TRUE))) {
    warning_(
      "varying() definitions of a determinstic channel
       {.var {deparse1(formula_lhs(formula))}} will be ignored."
    )
  }
  past_str <- strsplit(form_past, ",")[[1]]
  na_str <- grepl("NA", past_str)
  past_str[na_str] <- NA
  list(
    formula = as.formula(paste0(form_resp, "~", form_def)),
    specials = list(
      past = past_str,
      rank = Inf
    ),
    fixed = integer(0),
    varying = integer(0)
  )
}

#' Process Response Variables for Deterministic Channels
#'
#' @param x a `character` vector of length 1.
#' @noRd
formula_response <- function(x) {
  if (grepl("factor\\(.*\\)", x, perl = TRUE)) {
    list(type = "factor",
         resp = gsub("factor\\((.*)\\)", "\\1", x, perl = TRUE))
  } else if (grepl("numeric\\(.*\\)", x, perl = TRUE)) {
    list(type = "numeric",
         resp = gsub("numeric\\((.*)\\)", "\\1", x, perl = TRUE))
  } else if (grepl("integer\\(.*\\)", x, perl = TRUE)) {
    list(type = "integer",
         resp = gsub("integer\\((.*)\\)", "\\1", x, perl = TRUE))
  } else if (grepl("logical\\(.*\\)", x, perl = TRUE)) {
    list(type = "logical",
         resp = gsub("logical\\((.*)\\)", "\\1", x, perl = TRUE))
  } else {
    warning_(c(
      "No type specified for deterministic channel {.var {x}}:",
      `i` = "Assuming type is {.cls numeric}."
    ))
    list(type = "numeric", resp = x)
  }
}

#' Computes All Specials Defined in a Formula in the Context of the Data
#'
#' @param formula A `dynamiteformula` object.
#' @param data A `data.table` containing the variables present in the special
#'   definitions in the formula.
#' @noRd
evaluate_specials <- function(formula, data) {
  lapply(seq_along(formula), function(i) {
    if (length(formula[[i]]$specials) > 0) {
      out <- list()
      for (spec in formula_special_funs) {
        spec_formula <- formula[[i]]$specials[[spec]]
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
#' @param special An `integer` vector of special variable indices.
#' @param vars `"variables"` attribute of a `terms` object.
#' @param term_labels `"term.labels"` attribute of a `terms` object.
#'
#' @noRd
get_special_term_indices <- function(special, vars, term_labels) {
  out <- integer(length(special))
  for (i in seq_along(special)) {
    v <- deparse1(vars[[special[i] + 1]])
    out[i] <- which(term_labels == v)
  }
  out
}

# Supported specials
formula_special_funs <- c(
  "offset",
  "trials"
)

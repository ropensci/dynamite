#' Get and separate all specials of a formula
#'
#' @param x A `formula` object
#'
#' @noRd
formula_specials <- function(x) {
  out <- list(formula = NULL, specials = NULL, coefs = NULL)
  xt <- terms(x, specials = formula_special_funs)
  xt_specials <- attr(xt, "specials")[formula_special_funs]
  xt_variables <- attr(xt, "variables")
  special_vars <- unlist(xt_specials)
  for (y in formula_special_funs) {
    if (!is.null(xt_specials[[y]])) {
      out$specials[[y]] <- xt_variables[[xt_specials[[y]] + 1]][[2]]
    }
  }
  if (length(special_vars) > 0) {
    x <- formula(drop.terms(xt, special_vars - 1, keep.response = TRUE))
  }
  xt <- terms(x, specials = c("fixed", "varying"))
  xt_specials <- attr(xt, "specials")[c("fixed", "varying")]
  xt_variables <- attr(xt, "variables")
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
    form_terms <- formula_terms(x)[-(special_vars - 1)]
  } else {
    form_terms <- formula_terms(x)
  }
  fixed_terms <- union(form_terms, fixed_terms)
  fixed_icpt <- attr(xt, "intercept") || fixed_icpt
  common_terms <- intersect(fixed_terms, varying_terms)
  if (length(common_terms) > 0) {
    stop_("Variables '", cs(common_terms), "' ",
          "specified as both time-constant and time-varying.")
  }
  full_terms <- c(fixed_terms, varying_terms)
  any_icpt <- fixed_icpt || varying_icpt
  if (fixed_icpt && varying_icpt) {
    warning_("Both time-independent and time-varying intercept specified. ",
             "Defaulting to time-varying intercept.")
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
    if (!any_icpt) {
      stop_("Invalid formula for response variable '", y, "', ",
            "there are no predictors nor an intercept")
    }
    x <- as.formula(paste0(y, "~ 1"))
  }
  out$formula <- x
  out$fixed <- which(full_terms %in% fixed_terms)
  out$has_fixed_intercept <- as.logical(fixed_icpt)
  out$varying <- which(full_terms %in% varying_terms)
  out$has_varying_intercept <- as.logical(varying_icpt)
  out
}

#' Process formulas for deterministic channels and get past value definitions
#'
#' @param x A `formula` object
#'
#' @noRd
formula_past <- function(formula) {
  formula_str <- deparse1(formula)
  form_comp <- regexpr(
    pattern = "^(?<resp>[^~]+) ~ (?<def>[^~]+?)(?: \\+ past\\((?<past>.+)\\)){0,1}$",
    text = formula_str,
    perl = TRUE
  )
  start <- attr(form_comp, "capture.start")
  end <- start + attr(form_comp, "capture.length") - 1
  form_resp <- substr(formula_str, start[1], end[1])
  form_def <- substr(formula_str, start[2], end[2])
  if (grepl("past\\(", form_def, perl = TRUE)) {
    stop_("Past values term must be the last term of the formula")
  }
  form_past <- substr(formula_str, start[3], end[3])
  form_both <- c(form_def, form_past)
  if (any(grepl("fixed\\(.+\\)", form_both))) {
    warning_("fixed() definitions of a determinstic channel for ",
             as.character(formula_lhs(formula)), " will be ignored")
  }
  if (any(grepl("varying\\(.+\\)", form_both))) {
    warning_("varying() definitions of a determinstic channel for ",
             as.character(formula_lhs(formula)), " will be ignored")
  }
  past_str <- strsplit(form_past, ",")[[1]]
  na_str <- grepl("NA", past_str)
  past_str[na_str] <- NA
  list(
    formula = as.formula(paste0(form_resp, "~", form_def)),
    specials = list(
      past = try_(past_str, type = "numeric"),
      rank = Inf
    ),
    fixed = integer(0),
    varying = integer(0)
  )
}

#' Computes all specials defined in a formula in the context of the data
#'
#' @param formula A `dynamiteformula` object
#' @param data A `data.table` containing the variables present in the special
#'   definitions in the formula
#'
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

formula_special_funs <- c(
  "offset",
  "trials"
)

# Get and separate all specials of a formula
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
    if (length(special_vars)) {
        x <- formula(drop.terms(xt, special_vars - 1, keep.response = TRUE))
    }
    xt <- terms(x, specials = c("fixed", "varying"))
    xt_specials <- attr(xt, "specials")[c("fixed", "varying")]
    xt_variables <- attr(xt, "variables")
    special_vars <- unlist(xt_specials)
    fixed_rhs <- character(0)
    varying_rhs <- character(0)
    if (!is.null(xt_specials[["fixed"]])) {
        fixed_form <- eval(xt_variables[[xt_specials[["fixed"]] + 1]][[2]])
        fixed_rhs <- formula_rhs(fixed_form)
    }
    if (!is.null(xt_specials[["varying"]])) {
        varying_form <- eval(xt_variables[[xt_specials[["varying"]] + 1]][[2]])
        varying_rhs <- formula_rhs(varying_form)
    }
    if (!is.null(special_vars)) {
        form_rhs <- formula_rhs(x)[-(special_vars - 1)]
    } else {
        form_rhs <- formula_rhs(x)
    }
    fixed_rhs <- union(form_rhs, fixed_rhs)
    common_rhs <- intersect(fixed_rhs, varying_rhs)
    if (length(common_rhs)) {
        stop_("Variables ", common_rhs, " specified as both time-constant and time-varying.")
    }
    full_rhs <- c(fixed_rhs, varying_rhs)
    if (!identical(full_rhs, form_rhs)) {
        x <- reformulate(termlabels = full_rhs,
                         response = xt_variables[[2]],
                         intercept = attr(xt, "intercept"))
    }
    out$fixed <- which(full_rhs %in% fixed_rhs)
    out$varying <- which(full_rhs %in% varying_rhs)
    out$formula <- x
    out
}

# Computes all specials defined in a formula in the context of the data
evaluate_specials <- function(formula, data) {
    lapply(seq_along(formula), function(i) {
        if (length(formula[[i]]$specials)) {
            out <- list()
            for (spec in formula_special_funs)
                if (!is.null(spec_formula <- formula[[i]]$specials[[spec]])) {
                    out[[spec]] <- with(data, eval(spec_formula))
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

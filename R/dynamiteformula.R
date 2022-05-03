#' Model formula for \pkg{dynamite}
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

#' Checks if argument is a \code{dynamiteformula} object
#'
#' @param x An \R object
#'
#' @export
is.dynamiteformula <- function(x) {
    inherits(x, "dynamiteformula")
}

# Checks if argument is a formula
is.formula <- function(x) {
    inherits(x, "formula")
}

#' @export
`+.dynamiteformula` <- function(e1, e2) {
    if (is.dynamiteformula(e1)) {
        out <- add_dynamiteformula(e1, e2)
    } else {
        stop_("Method '+.dynamiteformula is not supported for ", class(e1), " objects.")
    }
    out
}

# Get all response variables of a dynamiteformula object
get_resp <- function(x) {
    sapply(x, "[[", "response")
}

# Get all predictor variables of a dynamiteformula object
get_pred <- function(x) {
    lapply(x, "[[", "predictors")
}

# Get all formulas of a dynamiteformula object
get_form <- function(x) {
    lapply(x, "[[", "formula")
}

# Check whether the formula contains an intercept
has_intercept <- function(x) {
    attr(terms(x$formula), "intercept") == 1
}

# Internal `+.dynamiteformula` for model constructions
add_dynamiteformula <- function(e1, e2) {
    if (is.dynamiteformula(e2)) {
        out <- join_dynamiteformulas(e1, e2)
    } else if (is.lags(e2)) {
        out <- set_lags(e1, e2)
    } else if (is.splines(e2)) {
        out <- set_splines(e1, e2)
    } else if (is.hiddenstates(e2)) {
        out <- set_hiddenstates(e1, e2)
    } else if (is.modeldata(e2)) {
        out <- set_modeldata(e1, e2)
    } else {
        stop_("Unable to add an object of class ", class(e2),
              " to an object of class 'dynamiteformula'")
    }
    out
}

# Join two model definitions and verify compatibility
join_dynamiteformulas <- function(e1, e2) {
    out <- c(e1, e2)
    resp_all <- get_resp(out)
    duped <- duplicated(resp_all)
    if (any(duped)) {
        stop_("Multiple definitions for response variables: ", resp_all[duped])
    }
    if (!is.null(attr(e1, "hidden")) && !is.null(attr(e2, "hidden"))) {
        stop_("Multiple definitions for hidden states")
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

# Set lag definitions for all channels
set_lags <- function(e1, e2) {
    if (!is.null(attr(e1, "lags"))) {
        stop_("Multiple definitions for lags")
    }
    attr(e1, "lags") <- e2
    e1
}

# Set the regression coefficient splines of the model
set_splines <- function(e1, e2) {
    if (!is.null(attr(e1, "splines")) && !attr(e2, "override")) {
        stop_("Multiple definitions for splines")
    }
    attr(e1, "splines") <- e2
    e1
}

# Set the hidden state process of the model
# set_hiddenstates <- function(e1, e2) {
#     if (!is.null(e1$hidden) && !attr(e2, "replace")) {
#         stop_("Multiple definitions for hidden states")
#     }
#     e1$hidden <- e2
#     e1
# }

# Set the data to be used by the model
# set_modeldata <- function(e1, e2) {
#     if (!is.null(e1$data) && !attr(e2, "replace")) {
#         stop_("Multiple data definitions")
#     }
#     e1$data <- e2
#     e1
# }

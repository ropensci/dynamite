#' Model formula for \pkg{btvcm}
#'
#' @export
#' @import formula.tools
btvcmformula <- function(formula, family, ...) {
    #all_terms <- terms(bf, specials = "t")
    #sp_terms <- attr(all_terms, "specials")
    # terms(bf)
    # term_labels <- labels(all_terms)
    # vars <- attr(all_terms, "variables")

    # Time-dependent coefficients
    # time_coefs <- integer(0)
    # time_terms <- sp_terms$t
    # if (!is.null(time_terms)) {
    #     time_forms <- lapply(vars[1 + time_terms], function(x) {
    #         y <- as.list(match.call(definition = t_, call = x))
    #         list(add = y$f, remove = x, type = if (is.null(y$type)) "rw" else y$type)
    #     })
    # }
    # out <- list(formula = bf, family = family, time_forms = time_forms)
    if (!is.btvcmfamily(family)) {
        stop_("Unsupported family object")
    }
    # TODO can we drop formula.tools dependency and implement these? Only lhs.vars and rhs.vars are used at the moment.
    lhs <- formula.tools::lhs.vars(formula)
    rhs <- formula.tools::rhs.vars(formula)
    out <- list(
        list(
            formula = formula,
            family = family,
            response = lhs,
            predictors = rhs
        )
    )
    class(out) <- c("btvcmformula")
    return(out)
}

#' @rdname btvcmformula
#' @export
obs <- btvcmformula

#' Checks if argument is a \code{btvcmformula} object
#'
#' @param x An \R object
#'
#' @export
is.btvcmformula <- function(x) {
    inherits(x, "btvcmformula")
}

# Checks if argument is a formula
is.formula <- function(x) {
    inherits(x, "formula")
}

#' @export
`+.btvcmformula` <- function(e1, e2) {
    if (is.btvcmformula(e1)) {
        out <- add_btvcmformula(e1, e2)
    } else {
        stop_("Method '+.btvcmformula is not supported for ", class(e1), " objects.")
    }
    out
}

# Get all response variables of a becauseformula object
get_resp <- function(x) {
    sapply(x, "[[", "response")
}

# Get all predictor variables of a becauseformula object
get_pred <- function(x) {
    lapply(x, "[[", "predictors")
}

# Internal `+.btvcmformula` for model constructions
add_btvcmformula <- function(e1, e2) {
    if (is.btvcmformula(e2)) {
        out <- join_btvcmformulas(e1, e2)
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
              " to an object of class 'btvcmformula'")
    }
    out
}

# Join two model definitions and verify compatibility
join_btvcmformulas <- function(e1, e2) {
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
    class(out) <- "btvcmformula"
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

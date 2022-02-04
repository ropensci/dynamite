#' Model formula for \pkg{because}
#'
#' @export
#' @import formula.tools
becauseformula <- function(formula, family = c("gaussian", "binomial", "poisson", "categorical"), ...) {
    # TODO may need more precise formulation for families
    bf <- formula
    family <- match.arg(family)
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
    lhs <- formula.tools::lhs.vars(bf)
    rhs <- formula.tools::rhs.vars(bf)
    out <- list(formulas = list(bf),
                families = list(family),
                resp = list(lhs),
                pred = list(rhs))
    class(out) <- c("becauseformula")
    return(out)
}

#' @rdname becauseformula
#' @export
obs <- becauseformula

#' Checks if argument is a \code{becauseformula} object
#'
#' @param x An \R object
#'
#' @export
is.becauseformula <- function(x) {
    inherits(x, "becauseformula")
}

# Checks if argument is a formula
is.formula <- function(x) {
    inherits(x, "formula")
}

#' @export
`+.becauseformula` <- function(e1, e2) {
    if (is.becauseformula(e1)) {
        out <- add_becauseformula(e1, e2)
    } else {
        stop_("Method '+.becauseformula is not supported for ", class(e1), " objects.")
    }
    out
}

# Internal `+.becauseformula` for model constructions
add_becauseformula <- function(e1, e2) {
    if (is.becauseformula(e2)) {
        out <- join_becauseformulas(e1, e2)
    } else if (is.splines(e2)) {
        out <- set_splines(e1, e2)
    } else if (is.hiddenstates(e2)) {
        out <- set_hiddenstates(e1, e2)
    } else if (is.modeldata(e2)) {
        out <- set_modeldata(e1, e2)
    } else {
        stop_("Unable to add an object of class ", class(e2), " to an object of class 'becauseformula'")
    }
    out
}

# Join two model definitions and verify compatibility
join_becauseformulas <- function(e1, e2) {
    out <- list(
        formulas = c(e1$formulas, e2$formulas),
        families = c(e1$families, e2$families),
        resp = c(e1$resp, e2$resp),
        pred = c(e1$pred, e2$pred)
    )
    uresp <- unlist(out$resp)
    duped <- duplicated(uresp)
    if (any(duped)) {
        stop_("Multiple definitions for response variables: ", uresp[duped])
    }
    if (!is.null(e1$splines) && !is.null(e1$splines)) {
        stop_("Multiple definitions for splines")
    }
    if (!is.null(e1$hidden) && !is.null(e2$hidden)) {
        stop_("Multiple definitions for hidden states")
    }
    if (!is.null(e1$data) && !is.null(e2$data)) {
        stop_("Multiple definitions for data")
    }
    class(out) <- "becauseformula"
    out
}

# Set the regression coefficient splines of the model
set_splines <- function(e1, e2) {
    if (!is.null(e1$splines) && !attr(e2, "replace")) {
        stop_("Multiple definitions for splines")
    }
    e1$splines <- e2
    e1
}

# Set the hidden state process of the model
set_hiddenstates <- function(e1, e2) {
    if (!is.null(e1$hidden) && !attr(e2, "replace")) {
        stop_("Multiple definitions for hidden states")
    }
    e1$hidden <- e2
    e1
}

# Set the data to be used by the model
set_modeldata <- function(e1, e2) {
    if (!is.null(e1$data) && !attr(e2, "replace")) {
        stop_("Multiple data definitions")
    }
    e1$data <- e2
    e1
}

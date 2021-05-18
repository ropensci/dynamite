#' Model formula for \pkg{because}
#'
#' @export
#' @import formula.tools
becauseformula <- function(formula, family = c("gaussian", "binomial", "poisson", "categorical"), data, ...) {
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
                data = list(data),
                resp = list(lhs),
                pred = list(rhs))
    class(out) <- c("becauseformula")
    return(out)
}

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

# Internal `+.becauseformula` for joining formulas and defining hidden states
add_becauseformula <- function(e1, e2) {
    if (is.becauseformula(e2)) {
        out <- join_becauseformulas(e1, e2)
    } else if (is.function(e2)) {
        res <- try(e2, silent = TRUE)
        if (is.hiddenstates(res)) {
            out <- set_hiddenstates(e1, res)
        } else {
            #TODO something informative
            stop_("Unable to add object to an object of class 'becauseformula'")
        }
    }
    out
}

join_becauseformulas <- function(e1, e2) {
    out <- list(
        formulas = c(e1$formulas, e2$formulas),
        families = c(e1$families, e2$families),
        data = c(e1$data, e2$data),
        resp = c(e1$resp, e2$resp),
        pred = c(e1$pred, e2$pred)
    )
    class(out) <- "becauseformula"
    out
}

set_hiddenstates <- function(e1, e2) {
    states <- try(e2, silent = TRUE)
    # TODO react to possible errors!
    e1$hidden <- states
    e1
}

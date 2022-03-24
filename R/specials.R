# Remove all specials from a formula
remove_specials <- function(x) {
    xt <- terms(x, specials = formula_special_funs)
    spec_attrs <- attr(xt, "specials")[formula_special_funs]
    spec_attrs[sapply(spec_attrs, is.null)] <- NULL
    specials <- unlist(spec_attrs)
    if (length(specials)) {
        formula(drop.terms(xt, specials - 1, keep.response = TRUE))
    } else {
        x
    }
}

# Get all specials of a formula
formula_specials <- function(x) {
    out <- list()
    xt <- terms(x, specials = formula_special_funs)
    xt_specials <- attr(xt, "specials")
    xt_variables <- attr(xt, "variables")
    for (spec in formula_special_funs) {
        spec_eval <- xt_specials[[spec]]
        if (!is.null(spec_eval)) {
            out[[spec]] <- xt_variables[[spec_eval + 1]][[2]]
        }
    }
    out
}

formula_special_funs <- c(
    "offset",
    "trials"
)

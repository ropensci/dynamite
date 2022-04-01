# Some stuff that might be useful

# Shortcut to update lists and vectors
`c<-` <- function(x, value) {
    c(x, value)
}

# Internal placeholder for match.call
t_ <- function(f, type = "rw") {
    invisible(NULL)
}

# Remove all whitespace characters from a character vector
rm_ws <- function(x) {
    y <- gsub("[ \t\r\n]", "", x, perl = TRUE)
    dim(y) <- dim(x)
    y
}

# Create a list with names defined by the arguments
named_list <- function(...) {
    y <- list(...)
    names(y) <- as.character(substitute(list(...)))[-1]
    y
}

# Check which elements of x have a prefix in y
# which_prefix <- function(x, y) {
#     ys <- paste0("^(", paste0(y, collapse = "|"), ").*$")
#     grep(pattern = ys, x, perl = TRUE)
# }

# Get left-hand side of a formula
formula_lhs <- function(x) {
    as.character(x[[2]])
}

# Get right-hand side of formula
formula_rhs <- function(x) {
    attr(terms(x), "term.labels")
}

# Collapse argument vector with a newline ignoring zero-length entries
collapse_rows <- function(x) {
    paste0(x[nzchar(x) > 0], collapse = "\n")
}

# Paste argument vectors with a newline
paste_rows <- function(...) {
    pasted <- sapply(list(...), function(x) {
        xlen <- length(x)
        if (xlen == 0) {
            ""
        } else if (xlen == 1) {
            x
        } else {
            paste0(x, collapse = "")
        }
    })
    pasted <- pasted[nzchar(pasted)]
    if (length(pasted)) {
        collapse_rows(pasted)
    } else {
        ""
    }
}

# Create an indenter
indenter_ <- function(m) {
    x <- rep(" ", m)
    idts <- sapply(1:10, function(y) {
        paste0(rep(x, y), collapse = "")
    })
    force(idts)
    function(n) {
        idts[[n]]
    }
}

# Check if x is a string of the form I(.)
is_as_is <- function(x) {
    grepl("I\\(.*\\)", x, perl = TRUE)
}

stop_ <- function(...) {
    stop(..., call. = FALSE)
}

warning_ <- function(...) {
    warning(..., call. = FALSE)
}

# Try to coerce one argument to specific type
try_ <- function(..., type) {
    if (missing(type)) {
        stop_("Argument 'type' must be given")
    }
    dots <- list(...)
    arg_val <- dots[[1]]
    if (is.null(arg_val)) {
        return(NULL)
    }
    arg_name <- names(dots)[1]
    as_type <- get(paste0("as.", type))
    out <- try(as_type(arg_val), silent = TRUE)
    if (identical(class(out), "try-error")) {
        stop_("Unable to coerce argument ", arg_name, " to '", type, "'")
    }
    out
}

# Shorthand for if (test) yes else no)
# NOTE: This is not the ifelse-function which enforces output length
ifelse_ <- function(test, yes, no) {
    if (test) {
        return(yes)
    }
    return(no)
}

# Return yes if test is TRUE, otherwise an empty vector of the same type
onlyif <- function(test, yes) {
    if (test) {
        return(yes)
    }
    do.call(paste0(typeof(yes)), args = list(length = 0))
}

# Combine model.matrix objects of all formulas of a btvcmformula into one
full_model.matrix <- function(formula, data) {
    model_matrices <- lapply(get_form(formula), model.matrix, data)
    model_matrix <- do.call(cbind, model_matrices)
    u_names <- unique(colnames(model_matrix))
    model_matrix <- model_matrix[, u_names, drop = FALSE]
    n_models <- length(model_matrices)
    attr(model_matrix, "assign") <- vector(mode = "list", length = n_models)
    attr(model_matrix, "fixed") <- vector(mode = "list", length = n_models)
    attr(model_matrix, "varying") <- vector(mode = "list", length = n_models)
    for (i in seq_along(model_matrices)) {
        attr(model_matrix, "assign")[[i]] <- which(u_names %in% colnames(model_matrices[[i]]))
        attr(model_matrix, "fixed")[[i]] <- which(attr(model_matrices[[i]], "assign") %in% formula[[i]]$fixed)
        attr(model_matrix, "varying")[[i]] <- which(attr(model_matrices[[i]], "assign") %in% formula[[i]]$varying)
    }
    model_matrix
}

# For prediction
full_model.matrix_fast <- function(formula, data, u_names) {
    model_matrices <- lapply(get_form(formula), model.matrix, data)
    model_matrix <- do.call(cbind, model_matrices)
    model_matrix[, u_names, drop = FALSE]
}

# Startup message for the package
.onAttach <- function(libname, pkgname) {
    # TODO
    invisible(NULL)
}

# Code to execute when loading the package
.onLoad <- function(libname, pkgname) {
    # TODO
    invisible(NULL)
}

# TODO there is ndraws method in posterior package, should probably define ndraws.btvcmfit
ndraws <- function(x) {
    (x$stanfit@sim$n_save[1] - x$stanfit@sim$warmup2[1]) * x$stanfit@sim$chains
}

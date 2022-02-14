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

# Get left-hand side of a formula as character
formula_lhs <- function(x) {
    as.character(x[[2]])
}

# Get right-hand side of fromula as a list with two components:
# @param vars the predictor terms as a character vector.
# @param cond terms after the vertical bar as a character vector.
formula_rhs <- function(x) {
    rhs_form <- as.character(x[[3]])
    out <- list(vars = trimws(strsplit(rhs_form[[2]], "+", fixed = TRUE)[[1]]), cond = NULL)
    if (length(rhs_form) == 3) {
        out$cond <- trimws(strsplit(rhs_form[[3]], "+", fixed = TRUE)[[1]])
    }
    out
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
    arg_name <- names(dots)[1]
    as_type <- get(paste0("as.", type))
    out <- try(as_type(arg_val), silent = TRUE)
    if (identical(class(out), "try-error")) {
        stop_("Unable to coerce argument ", arg_name, " to '", type, "'")
    }
    out
}

# Starup message for the package
.onAttach <- function(libname, pkgname) {
    # TODO
    invisible(NULL)
}

# Code to execute when loading the package
.onLoad <- function(libname, pkgname) {
    # TODO
    invisible(NULL)
}

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

# Get all specials of a formula
formula_specials <- function(x) {
    out <- list()
    offset <- get_offset(x)
    if (!is.null(offset)) {
        out$offset <- offset
    }
    out
}

# Get offset term from formula
get_offset <- function(x) {
    xt <- terms(x)
    if (!is.null(offset <- attr(xt, "offset"))) {
        return(attr(xt, "variables")[[offset + 1]][[2]])
    }
    return(NULL)
}

# Collapse argument vector with a newline
collapse_rows <- function(x) {
    paste0(x, collapse = "\n")
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
        character(0)
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

# TODO there is ndraws method in posterior package, should probably define ndraws.btvcmfit
ndraws <- function(x) {
    (x$stanfit@sim$n_save[1] - x$stanfit@sim$warmup2[1]) * x$stanfit@sim$chains
}

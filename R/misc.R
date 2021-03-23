# Some stuff that might be useful

# Shortcut to update lists and vectors
`c<-` <- function(x, value) {
    c(x, value)
}

stop_ <- function(...) {
    stop(..., call. = FALSE)
}

warning_ <- function(...) {
    warning(..., call. = FALSE)
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
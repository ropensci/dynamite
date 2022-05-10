# Some stuff that might be useful

# This is needed to use tidyselect's where which is not exported.
utils::globalVariables("where")


#' A Shortcut function to update lists and vectors
#'
#' @param x A list or a vector
#' @param value An element to add to the list or vector
#'
#' @noRd
`c<-` <- function(x, value) {
  c(x, value)
}

# Internal placeholder for match.call
# t_ <- function(f, type = "rw") {
#   invisible(NULL)
# }

# TODO check if used somewhere
#' Remove all whitespace characters from a character vector
#'
#' @param x A character vector
#'
#' @noRd
rm_ws <- function(x) {
  y <- gsub("[ \t\r\n]", "", x, perl = TRUE)
  dim(y) <- dim(x)
  y
}

#' Create a list with names defined by the arguments
#'
#' @param ... Any number of `name = value` pairs
#'
#' @noRd
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

#' Get the left-hand side of a formula
#'
#' @param x A `formula` object
#'
#' @noRd
formula_lhs <- function(x) {
  if (length(x) == 3) {
    x[[2]]
  } else {
    NULL
  }
}

#' Get the right-hand side of formula
#'
#' @param x A `formula` object
#'
#' @noRd
formula_rhs <- function(x) {
  if (length(x) == 3) {
    x[[3]]
  } else {
    x[[2]]
  }
}

#' Get the right-hand side terms of a formula
#'
#' @param x A `formula` object
#'
#' @noRd
formula_terms <- function(x) {
  attr(terms(x), "term.labels")
}

#' Replace terms in a formula based on a regular expression
gsub_formula <- function(pattern, replacement, formula, ...) {
  formula_str <- deparse(formula)
  as.formula(gsub(pattern, replacement, formula_str, ...))
}

# Collapse argument vector with a newline ignoring zero-length entries
# collapse_rows <- function(x) {
#     paste0(x[nzchar(x)], collapse = "\n")
# }

#' Create a comma-separated character string to represent a Stan integer array
#'
#' @param x A character vector
#'
#' @noRd
cs <- function(x) {
  paste0(x, collapse = ", ")
}

#' Paste and optionally parse character strings containing glue syntax
#'
#' @param ... Any number of character vectors of arbitrary length
#' @param .indent A character string to prefix each row with
#' @param .parse A logical value indicating whether glue syntax should be
#'   parsed by [glue::glue()].
#'
#' @noRd
paste_rows <- function(..., .indent = "", .parse = TRUE) {
  dots <- list(...)
  ndots <- length(dots)
  if (ndots) {
    idt_vec <- character(ndots)
    idt_vec[1:ndots] <- .indent
  }
  pasted <- rep("", ndots)
  for (i in seq_len(ndots)) {
    x <- dots[[i]]
    xlen <- length(x)
    if (xlen == 0) {
      pasted[i] <- ""
    } else if (xlen == 1) {
      if (nzchar(x)) {
        if (.parse) {
          xglue <- glue::glue(x, .envir = parent.frame(), .trim = FALSE)
          pasted[i] <- paste0(idt_vec[i], xglue, collapse = "\n")
        } else {
          pasted[i] <- paste0(idt_vec[i], x)
        }
      }
    } else {
      x <- x[nzchar(x)]
      if (length(x) > 0) {
        if (.parse) {
          xglue <- sapply(x, function(y) {
            paste0(glue::glue(y, .envir = parent.frame(), .trim = FALSE),
                   collapse = "\n")
          })
          pasted[i] <- paste0(idt_vec[i], xglue, collapse = "\n")
        } else {
          pasted[i] <- paste0(idt_vec[i], x, collapse = "\n")
        }
      }
    }
  }
  pasted <- pasted[nzchar(pasted)]
  if (length(pasted) > 0) {
    paste0(pasted, collapse = "\n")
  } else {
    ""
  }
}

#' Create an indenter
#'
#' @param m An integer denoting how many spaces does one unit of indentation
#'   correspond to (default = 2)
#'
#' @noRd
indenter_ <- function(m = 2) {
  x <- rep(" ", m)
  idts <- sapply(0:10, function(y) {
    paste0(rep(x, y), collapse = "")
  })
  force(idts)
  function(v) {
    unlist(idts[v + 1])
  }
}

#' Check if x is a string of the form I(.)
#'
#' @param x A character string
#'
#' @noRd
is_as_is <- function(x) {
  grepl("I\\(.*\\)", x, perl = TRUE)
}

#' Stop function execution without displaying the call
#'
#' @param ... See [stop()]
#'
#' @noRd
stop_ <- function(...) {
  stop(..., call. = FALSE)
}

#' Generate a warning message without displaying the call
#'
#' @param ... See [warning()]
#'
#' @noRd
warning_ <- function(...) {
  warning(..., call. = FALSE)
}

#' Try to coerce one argument to specific type
#'
#' @param ... A `name = value` pair,
#'   only the first one is used if multiple pairs are supplied.
#'
#' @noRd
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

#' Shorthand for if (test) yes else no
#'
#' @param test An object which can be coerced into logical mode
#' @param yes An \R object to return when `test` evaluates to `TRUE`
#' @param no An \R object to return when `test` evaluates to `FALSE`
#'
#' @noRd
ifelse_ <- function(test, yes, no) {
  if (test) {
    yes
  } else {
    no
  }
}

#' Return yes if test is TRUE, otherwise an empty vector of the same type
#'
#' @param tes An object which can be coerced into logical mode
#' @param yes An \R object to return when `test` evaluates to `TRUE`
#'
#' @noRd
onlyif <- function(test, yes) {
  if (test) {
    yes
  } else {
    do.call(paste0(typeof(yes)), args = list(length = 0))
  }
}

#' Convert a formula into an expression of mathematical equality
#'
#' @param formula A `formula` object without a response
#'
#' @noRd

formula2expression <- function(formula) {
  str2expression(deparse(formula))
}

#' Evaluate a formula as an expression
#'
#' @param formula A `formula` object
#'
#' @noRd
eval_formula <- function(formula, envir) {
  eval(expr = formula2expression(formula_rhs(formula)), envir = envir)
}

# TODO is needed?
# Startup message for the package
.onAttach <- function(libname, pkgname) {
  # TODO
  invisible(NULL)
}

# TODO is needed?
# Code to execute when loading the package
.onLoad <- function(libname, pkgname) {
  # TODO
  invisible(NULL)
}

# TODO there is ndraws method in posterior package, should probably define ndraws.dynamitefit
ndraws <- function(x) {
  (x$stanfit@sim$n_save[1] - x$stanfit@sim$warmup2[1]) * x$stanfit@sim$chains
}

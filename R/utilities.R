# This is needed to use tidyselect's where which is not exported.
utils::globalVariables(c(".", ".I", ".N", ".SD", "where"))

# Data table awareness
.datatable.aware = TRUE

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
#'
#' @param pattern See [gsub()]
#' @param replacement See [gsub()]
#' @param formula A `formula` object
#' @param ... See [gsub()]
#'
#' @noRd
gsub_formula <- function(pattern, replacement, formula, ...) {
  formula_str <- deparse1(formula)
  as.formula(gsub(pattern, replacement, formula_str, ...))
}

#' Add fixed or varying terms to a formula
#'
#' @param formula A `formula` object
#' @param x A `character` vector of terms to add
#' @param type Either `"fixed"` or `"varying"` indicating type of terms to add
#' @param varying_idx Indices of left-hand side terms that have
#'   time-varying coefficients
#' @param varying_icpt Does the formula have a varying intercept?
#'
#' @srrstats {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @noRd
increment_formula <- function(formula, x, type = c("fixed", "varying"),
                              varying_idx, varying_icpt, fixed_icpt) {
  type <- match.arg(type)
  v_icpt <- ifelse_(varying_icpt, "1", "-1")
  n_varying <- length(varying_idx)
  x_plus <- paste0(x, collapse = " + ")
  ft <- terms(formula)
  tr <- attr(ft, "term.labels")
  v <- paste0(tr[varying_idx], collapse = " + ")
  formula_str <- ""
  if (n_varying > 0) {
    if (n_varying < length(tr)) {
      formula <- drop.terms(ft, dropx = varying_idx, keep.response = TRUE)
      ft <- terms(formula)
      tr <- attr(ft, "term.labels")
    } else {
      tr <- character(0)
    }
  }
  if (length(tr) > 0) {
    formula <- reformulate(tr,
                           response = formula_lhs(formula),
                           intercept = fixed_icpt)
    formula_str <- deparse1(formula)
  } else {
    formula_str <- paste0(
      deparse1(formula_lhs(formula)),
      " ~ ",
      ifelse_(fixed_icpt, "1", "-1")
    )
  }
  if (identical(type, "varying")) {
    v <- paste0(" + ", v, ifelse_(nzchar(v), " + ", ""), x_plus)
    out_str <- glue::glue("{formula_str} + varying(~{v_icpt}{v})")
  } else {
    if (nzchar(v)) {
      v <- paste0(" + ", v)
    }
    if (n_varying > 0 || varying_icpt) {
      out_str <- glue::glue("{formula_str} + {x_plus} + varying(~{v_icpt}{v})")
    } else {
      out_str <- glue::glue("{formula_str} + {x_plus}")
    }
  }
  as.formula(out_str)
}

increment_formula_determnistic <- function(formula, x) {
  as.formula(paste0(deparse1(formula), " + ", x))
}

#' Create a comma-separated character string to represent a Stan integer array
#'
#' @param x A `character` vector
#'
#' @noRd
cs <- function(x) {
  paste0(x, collapse = ", ")
}

#' Paste and optionally parse character strings containing glue syntax
#'
#' @param ... Any number of `character` vectors of arbitrary length
#' @param .indent A `character` string to prefix each row with
#' @param .parse A `logical` value indicating whether glue syntax should be
#'   parsed by [glue::glue()].
#'
#' @noRd
paste_rows <- function(..., .indent = "", .parse = TRUE) {
  dots <- list(...)
  ndots <- length(dots)
  if (ndots) {
    idt_vec <- character(ndots)
    idt_vec[1L:ndots] <- .indent
  }
  pasted <- rep("", ndots)
  for (i in seq_len(ndots)) {
    x <- dots[[i]]
    xlen <- length(x)
    if (identical(xlen, 0L)) {
      pasted[i] <- ""
    } else if (identical(xlen, 1L)) {
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
      if (length(x) > 0L) {
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
  if (length(pasted) > 0L) {
    paste0(pasted, collapse = "\n")
  } else {
    ""
  }
}

#' Create an indenter
#'
#' @param m An integer denoting how many spaces does one unit of indentation
#'   correspond to (default = 2L)
#'
#' @noRd
indenter_ <- function(m = 2L) {
  x <- rep(" ", m)
  idts <- sapply(0L:10L, function(y) {
    paste0(rep(x, y), collapse = "")
  })
  force(idts)
  function(v) {
    unlist(idts[v + 1L])
  }
}

#' Check if x contains a string of the form I(.)
#'
#' @param x A character string
#'
#' @noRd
has_as_is <- function(x) {
  grepl("I\\(.+\\)", x, perl = TRUE)
}

#' Fill gaps (NAs) in a vector with the last non-NA observation
#'
#' @param x A vector possibly containing NA values
#'
#' @noRd
locf <- function(x) {
  non_na <- ifelse_(is.na(x[1L]), c(1L, which(!is.na(x))), which(!is.na(x)))
  fill <- diff(c(non_na, length(x) + 1L))
  rep(x[non_na], fill)
}

#' Stop function execution without displaying the call
#'
#' @param message See [cli::cli_abort()]
#' @param ... See [cli::cli_abort()]
#' @param call See See [cli::cli_abort()]
#'
#' @noRd
stop_ <- function(message, ..., call = rlang::caller_env()) {
  cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
}

#' Generate a warning message
#'
#' @param message See [cli::cli_abort()]
#' @param ... See [cli::cli_abort()]
#'
#' @noRd
warning_ <- function(message, ...) {
  cli::cli_warn(message, ..., .envir = parent.frame())
}

#' Generate an informative message
#'
#' @param message See [cli::cli_abort()]
#' @param ... See [cli::cli_abort()]
#'
#' @noRd
message_ <- function(message, ...) {
  cli::cli_inform(message, ..., .envir = parent.frame())
}

#' Try to coerce one argument to specific type
#'
#' @param ... A `name = value` pair,
#'   only the first one is used if multiple pairs are supplied.
#'
#' @noRd
try_ <- function(..., type) {
  if (missing(type)) {
    stop_("Argument {.var type} must be given.")
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
    stop_("Unable to coerce argument {.var {arg_name}} to {.cls {type}}.")
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

# Placeholder for future
# Startup message for the package
# .onAttach <- function(libname, pkgname) {
#   # TODO
#   #invisible(NULL)
# }

# Placeholder for future
# Code to execute when loading the package
# .onLoad <- function(libname, pkgname) {
#   # TODO
#   #invisible(NULL)
# }

# TODO there is ndraws method in posterior package, should probably define ndraws.dynamitefit
ndraws <- function(x) {
  (x$stanfit@sim$n_save[1] - x$stanfit@sim$warmup2[1]) * x$stanfit@sim$chains
}

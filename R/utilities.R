# This is needed to use tidyselect's where which is not exported.
utils::globalVariables(c(".", ".I", ".N", ".SD", "where"))

# Data table awareness
.datatable.aware <- TRUE

#' Get the left-hand side of a formula
#'
#' @param x A `formula` object.
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
#' @param x A `formula` object.
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
#' @param x A `formula` object.
#' @noRd
formula_terms <- function(x) {
  attr(terms(x), "term.labels")
}

#' Replace terms in a formula based on a regular expression
#'
#' @param pattern See [base::gsub()].
#' @param replacement See [base::gsub()].
#' @param formula A `formula` object
#' @param ... Additional arguments passed to [base::gsub()]
#' @noRd
gsub_formula <- function(pattern, replacement, formula, ...) {
  formula_str <- deparse1(formula)
  as.formula(gsub(pattern, replacement, formula_str, ...))
}

#' Add fixed or varying terms to a formula
#'
#' @param formula A `formula` object.
#' @param x A `character` vector of terms to add.
#' @param type Either `"fixed"` or `"varying"` indicating type of terms to add.
#' @param varying_idx Indices of left-hand side terms that have
#'   time-varying coefficients
#' @param varying_icpt Does the formula have a varying intercept?
#' @param fixed_icpt Does the formula have a fixed intercept?
#' @noRd
increment_formula <- function(formula, x, type = c("fixed", "varying"),
                              varying_idx, varying_icpt, fixed_icpt) {
  v_icpt <- ifelse_(varying_icpt, "1", "-1")
  n_varying <- length(varying_idx)
  x_plus <- paste0(x, collapse = " + ")
  ft <- terms(formula)
  tr <- attr(ft, "term.labels")
  v <- paste0(tr[varying_idx], collapse = " + ")
  formula_str <- ""
  if (n_varying > 0L) {
    if (n_varying < length(tr)) {
      formula <- drop.terms(ft, dropx = varying_idx, keep.response = TRUE)
      ft <- terms(formula)
      tr <- attr(ft, "term.labels")
    } else {
      tr <- character(0)
    }
  }
  if (length(tr) > 0) {
    formula <- reformulate(
      termlabels = tr,
      response = formula_lhs(formula),
      intercept = fixed_icpt
    )
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
    if (n_varying > 0L || varying_icpt) {
      out_str <- glue::glue("{formula_str} + {x_plus} + varying(~{v_icpt}{v})")
    } else {
      out_str <- glue::glue("{formula_str} + {x_plus}")
    }
  }
  as.formula(out_str)
}

#' Create a comma-separated character string to represent a Stan integer array
#'
#' @param x A `character` vector.
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
#' @noRd
paste_rows <- function(..., .indent = "", .parse = TRUE) {
  dots <- list(...)
  ndots <- length(dots)
  if (ndots) {
    idt_vec <- character(ndots)
    idt_vec[seq.int(1L, ndots)] <- .indent
  }
  pasted <- rep("", ndots)
  for (i in seq_len(ndots)) {
    x <- dots[[i]]
    x <- x[nzchar(x)]
    xlen <- length(x)
    if (identical(xlen, 0L)) {
      pasted[i] <- ""
    } else if (identical(xlen, 1L)) {
      if (.parse) {
        xglue <- glue::glue(x, .envir = parent.frame(), .trim = FALSE)
        pasted[i] <- paste0(idt_vec[i], xglue, collapse = "\n")
      } else {
        pasted[i] <- paste0(idt_vec[i], x)
      }
    } else {
      if (.parse) {
        xglue <- vapply(
          x,
          function(y) {
            paste0(
              glue::glue(y, .envir = parent.frame(), .trim = FALSE),
              collapse = "\n"
            )
          },
          character(1L)
        )
        pasted[i] <- paste0(idt_vec[i], xglue, collapse = "\n")
      } else {
        pasted[i] <- paste0(idt_vec[i], x, collapse = "\n")
      }
    }
  }
  pasted <- pasted[nzchar(pasted)]
  ifelse_(
    length(pasted) > 0L,
    paste0(pasted, collapse = "\n"),
    ""
  )
}

#' Create an indenter
#'
#' @param m An integer denoting how many spaces does one unit of indentation
#'   correspond to (default = 2L)
#' @noRd
indenter_ <- function(m = 2L) {
  x <- rep(" ", m)
  idts <- vapply(0L:10L, function(y) {
    paste0(rep(x, y), collapse = "")
  }, character(1L))
  force(idts)
  function(v) {
    unlist(idts[v + 1L])
  }
}

#' Check if x contains a string of the form I(.)
#'
#' @param x A `character` vector.
#' @noRd
has_as_is <- function(x) {
  grepl("I\\(.+\\)", x, perl = TRUE)
}

#' Fill gaps (NAs) in a vector with the last non-NA observation
#'
#' @param x A vector possibly containing NA values
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
#' @noRd
stop_ <- function(message, ..., call = rlang::caller_env()) {
  cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
}

stopifnot_ <- function(cond, message, ..., call = rlang::caller_env()) {
  if (!cond) {
    cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
  }
}

#' Generate a warning message
#'
#' @param message See [cli::cli_abort()]
#' @param ... See [cli::cli_abort()]
#' @noRd
warning_ <- function(message, ...) {
  cli::cli_warn(message, ..., .envir = parent.frame())
}

#' Generate an informative message
#'
#' @param message See [cli::cli_abort()]
#' @param ... See [cli::cli_abort()]
#' @noRd
message_ <- function(message, ...) {
  cli::cli_inform(message, ..., .envir = parent.frame())
}

#' Shorthand for `if (test) yes else no`
#'
#' @param test An object which can be coerced into `logical`.
#' @param yes An \R object to return when `test` evaluates to `TRUE`.
#' @param no An \R object to return when `test` evaluates to `FALSE`.
#' @noRd
ifelse_ <- function(test, yes, no) {
  if (test) {
    yes
  } else {
    no
  }
}

#' Return yes if test is TRUE, otherwise return NULL
#'
#' @param tes An object which can be coerced into `logical`.
#' @param yes An \R object to return when `test` evaluates to `TRUE`.
#' @noRd
onlyif <- function(test, yes) {
  if (test) {
    yes
  } else {
    NULL
  }
}

#' Adds NA gaps to fill in missing time points in a data frame
#'
#' @inheritParams dynamite
#' @noRd
fill_time <- function(data, group_var, time_var) {
  if (is.factor(data[[time_var]])) {
    warning_("Time indexing variable {.arg {time_var}} is a {.cls factor},
      converting the variable to {.cls integer} based on its levels.")
    data[[time_var]] <- as.integer(data[[time_var]])
  }
  time <- sort(unique(data[[time_var]]))
  stopifnot_(
    length(time) > 1L,
    "There must be at least two time points in the data."
  )
  time_ivals <- diff(time)
  time_scale <- min(diff(time))
  if (any(time_ivals[!is.na(time_ivals)] %% time_scale > 0)) {
    stop_("Observations must occur at regular time intervals.")
  } else {
    full_time <- seq(time[1L], time[length(time)], by = time_scale)
    groups <- !is.null(group_var)
    if (groups) {
      time_groups <- data |>
        dplyr::group_by(.data[[group_var]]) |>
        dplyr::summarise(has_missing = !identical(.data[[time_var]], full_time))
      if (any(time_groups$has_missing)) {
        full_data_template <- expand.grid(
          time = full_time,
          group = unique(data[[group_var]])
        )
        names(full_data_template) <- c(time_var, group_var)
        data <- full_data_template |>
          dplyr::left_join(data, by = c(group_var, time_var))
      }
    } else {
      if (!identical(data[[time_var]], full_time)) {
        full_data_template <- data.frame(time = full_time)
        names(full_data_template) <- time_var
        data <- full_data_template |>
          dplyr::left_join(data, by = time_var)
      }
    }
  }
  data
}

# Placeholder for future
# Startup message for the package
# .onAttach <- function(libname, pkgname) {
#   #invisible(NULL)
# }

# Placeholder for future
# Code to execute when loading the package
# .onLoad <- function(libname, pkgname) {
#   #invisible(NULL)
# }
#' Return the Number of Posterior Draws of the dynamitefit Object
#' @export
#' @export ndraws
#' @rdname ndraws-dynamitefit
#' @method ndraws dynamitefit
#' @examples
#' ndraws(gaussian_examplefit)
ndraws.dynamitefit <- function(x) {
  as.integer(
    (x$stanfit@sim$n_save[1] - x$stanfit@sim$warmup2[1]) * x$stanfit@sim$chains
  )
}

# For data.table and tidyselect
utils::globalVariables(c(".", ".I", ".N", ".SD", "where"))

# Data table awareness
.datatable.aware <- TRUE

#' Get the Left-hand Side of a Formula
#'
#' @param x A `formula` object.
#' @noRd
formula_lhs <- function(x) {
  if (identical(length(x), 3L)) {
    x[[2]]
  } else {
    NULL
  }
}

#' Get the Right-hand Side of a Formula
#'
#' @param x A `formula` object.
#' @noRd
formula_rhs <- function(x) {
  if (identical(length(x), 3L)) {
    x[[3]]
  } else {
    x[[2]]
  }
}

#' Get the Right-hand Side Terms of a Formula
#'
#' @param x A `formula` object.
#' @noRd
formula_terms <- function(x) {
  attr(terms(x), "term.labels")
}

#' Wrapper of `drop.terms` with Formula Output And Empty RHS Support
#'
#' @inheritParams drop.terms
#' @noRd
drop_terms <- function(termobj, dropx) {
  dropx_len <- length(dropx)
  if (identical(dropx_len, 0L)) {
    formula(termobj)
  } else {
    labs <- attr(termobj, "term.labels")
    if (length(labs) > dropx_len) {
      formula(drop.terms(termobj, dropx, keep.response = TRUE))
    } else {
      icpt <- attr(termobj, "intercept")
      resp <- attr(termobj, "variables")[[attr(termobj, "response") + 1L]]
      as.formula(paste0(as.character(resp), "~", ifelse_(icpt, "1", "-1")))
    }
  }
}

#' Drops Variables from the Data That Are Not Used by Any Model Formula
#'
#' @inheritParams dynamite
#' @noRd
drop_unused <- function(dformula, data, group_var, time_var) {
  used <- c(group_var, time_var)
  for (i in seq_along(dformula)) {
    used <- c(used, all.vars(dformula[[i]]$original))
  }
  unused <- setdiff(names(data), used)
  for (u in unused) {
    data[, (u) := NULL]
  }
}

#' Add Fixed or Varying Terms to a Formula
#'
#' @param formula \[`formula`]\cr A formula object.
#' @param x \[`character()`]\cr A vector of terms to add.
#' @param type \[`character(1)`]\cr Either `"fixed"` or `"varying"`
#'   indicatingthe  type of terms to add.
#' @param varying_idx \[`integer()`] Indices of left-hand side terms that have
#'   time-varying coefficients
#' @param varying_icpt \[`logical(1)`] Does the formula have a varying
#'   intercept?
#' @param fixed_icpt \[`logical(1)`] Does the formula have a fixed intercept?
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
      tr <- character(0L)
    }
  }
  if (length(tr) > 0L) {
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
    v <- ifelse_(nzchar(v), paste0(" + ", v), v)
    out_str <- ifelse_(
      n_varying > 0L || varying_icpt,
      glue::glue("{formula_str} + {x_plus} + varying(~{v_icpt}{v})"),
      glue::glue("{formula_str} + {x_plus}")
    )
  }
  as.formula(out_str)
}

#' Create a Comma-separated Character String
#'
#' @param x A `character` vector.
#' @noRd
cs <- function(x) {
  paste0(x, collapse = ", ")
}

#' Paste And Optionally Parse Character Strings Containing `glue` Syntax
#'
#' @param ... Any number of `character` vectors of arbitrary length.
#' @param .indent \[`character(1)`]\cr A string to prefix each row with
#' @param .parse \[`logical(1)`]\cr Should `glue` syntaxbe parsed
#'   by [glue::glue()].
#' @noRd
paste_rows <- function(..., .indent = "", .parse = TRUE) {
  dots <- list(...)
  ndots <- length(dots)
  if (ndots > 0L) {
    idt_vec <- character(ndots)
    idt_vec[seq_len(ndots)] <- .indent
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

#' Create an Indentation Function
#'
#' @param m \[`integer(1)`]\cr An integer denoting how many spaces does one
#'   unit of indentation correspond to (default is 2).
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

#' Fill Gaps (NAs) in a Vector with the Last Non-NA Observation
#'
#' @param x \[`vector()`]\cr A vector possibly containing NA values.
#' @noRd
locf <- function(x) {
  non_na <- ifelse_(is.na(x[1L]), c(1L, which(!is.na(x))), which(!is.na(x)))
  fill <- diff(c(non_na, length(x) + 1L))
  rep(x[non_na], fill)
}

#' Stop Function Execution Without Displaying the Call
#'
#' @param message See [cli::cli_abort()].
#' @param ... See [cli::cli_abort()].
#' @param call See See [cli::cli_abort()].
#' @noRd
stop_ <- function(message, ..., call = rlang::caller_env()) {
  cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
}

#' Stop Function Execution Unless Condition Is True
#' @inheritParams stop_
#' @param cond \[`logical(1)`] Condition to evaluate.
#' @noRd
stopifnot_ <- function(cond, message, ..., call = rlang::caller_env()) {
  if (!cond) {
    cli::cli_abort(message, ..., .envir = parent.frame(), call = call)
  }
}

#' Generate a Warning Message
#'
#' @param message See [cli::cli_warn()].
#' @param ... See [cli::cli_warn()].
#' @noRd
warning_ <- function(message, ...) {
  cli::cli_warn(message, ..., .envir = parent.frame())
}

#' Generate an Informative Message
#'
#' @param message See [cli::cli_inform()]
#' @param ... See [cli::cli_inform()]
#' @noRd
message_ <- function(message, ...) {
  cli::cli_inform(message, ..., .envir = parent.frame())
}

#' Shorthand for `if (test) yes else no`
#'
#' @param test \[`logical(1)`] Condition to evaluate.
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

#' Return `yes` If `test` Is `TRUE`, Otherwise Return `NULL`
#'
#' @param test \[`logical(1)`] Condition to evaluate.
#' @param yes An \R object to return when `test` evaluates to `TRUE`.
#' @noRd
onlyif <- function(test, yes) {
  if (test) {
    yes
  } else {
    NULL
  }
}

#' Number of Unique Values
#'
#' @inheritParams data.table::uniqueN
#' @noRd
n_unique <- data.table::uniqueN

#' Is Categorical Logit GLM Supported By Current Stan Version
#'
#' @noRd
stan_supports_categorical_logit_glm <- function(backend) {
  backend_version <- ifelse_(
    backend == "rstan",
    as.character(rstan::stan_version()),
    as.character(cmdstanr::cmdstan_version())
  )
  utils::compareVersion(backend_version, "2.23") >= 0
}
#' Row-wise log-sum-exp
#'
#' @noRd
log_sum_exp <- function(x) {
  maxs <- apply(x, 1, max)
  maxs + log(rowSums(exp(x - maxs)))
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

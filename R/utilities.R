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
    x[[2L]]
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
    x[[3L]]
  } else {
    x[[2L]]
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
      formula(stats::drop.terms(termobj, dropx, keep.response = TRUE))
    } else {
      icpt <- attr(termobj, "intercept")
      resp <- attr(termobj, "variables")[[attr(termobj, "response") + 1L]]
      as.formula(paste0(as.character(resp), "~", ifelse_(icpt, "1", "-1")))
    }
  }
}

#' Faster Column Bind for `data.table` Objects
#'
#' @param ... `data.table` objects.
#' @noRd
cbind_datatable <- function(...) {
  data.table::setattr(
    do.call(c, list(...)), "class", c("data.table", "data.frame")
  )
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
#'   indicating the type of terms to add.
#' @param varying_idx \[`integer()`]\cr Indices of left-hand side terms that
#'   have time-varying coefficients.
#' @param random_idx \[`integer()`]\cr Indices of the left-hand side term that
#'   have random coefficient terms.
#' @param varying_icpt \[`logical(1)`]\cr Does the formula have a varying
#'   intercept?
#' @param fixed_icpt \[`logical(1)`]\cr
#'   Does the formula have a fixed intercept?
#' @param random_icpt \[`logical(1)`]\cr
#'   Does the formula have a random intercept
#' @noRd
increment_formula <- function(formula, x,
                              type = c("fixed", "varying", "random"),
                              varying_idx, fixed_idx, random_idx,
                              varying_icpt, fixed_icpt, random_icpt) {
  v_icpt <- ifelse_(varying_icpt, "1", "-1")
  r_icpt <- ifelse_(random_icpt, "1", "-1")
  n_varying <- length(varying_idx)
  n_random <- length(random_idx)
  x_plus <- paste0(x, collapse = " + ")
  ft <- terms(formula)
  tr <- attr(ft, "term.labels")
  f <- paste0(tr[fixed_idx], collapse = " + ")
  v <- paste0(tr[varying_idx], collapse = " + ")
  r <- paste0(tr[random_idx], collapse = " + ")
  resp <- attr(ft, "variables")[[2L]]
  formula_str <- ifelse_(
    fixed_icpt,
    paste0(resp, "~1"),
    paste0(resp, "~-1")
  )
  v_out <- ifelse_(
    identical(type, "varying"),
    paste0(" + ", v, ifelse_(nzchar(v), " + ", ""), x_plus),
    ifelse_(nzchar(v), paste0(" + ", v), v)
  )
  v_out <- ifelse_(
    identical(type, "varying") || n_varying > 0L || varying_icpt,
    glue::glue(" + varying(~{v_icpt}{v_out})"),
    ""
  )
  r_out <- ifelse_(
    identical(type, "random"),
    paste0(" + ", r, ifelse_(nzchar(r), " + ", ""), x_plus),
    ifelse_(nzchar(r), paste0(" + ", r), r)
  )
  r_out <- ifelse_(
    identical(type, "random") || n_random > 0L || random_icpt,
    glue::glue(" + random(~{r_icpt}{r_out})"),
    ""
  )
  f_out <- ifelse_(
    identical(type, "fixed"),
    paste0(" + ", f, ifelse_(nzchar(f), " + ", ""), x_plus),
    ifelse_(nzchar(f), paste0(" + ", f), f)
  )
  out_str <- glue::glue("{formula_str}{f_out}{v_out}{r_out}")
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
#' @param .parse \[`logical(1)`]\cr Should `glue` syntax be parsed
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
  idts <- vapply(
    seq.int(0L, 10L),
    function(y) paste0(rep(x, y), collapse = ""),
    character(1L)
  )
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

#' Computes a Topological Ordering for the Vertices of a DAG.
#'
#' @param A \[`matrix`]\cr An adjacency matrix.
#' @return An `integer` vector giving a topological order of the vertices.
#' @noRd
topological_order <- function(A) {
  n <- ncol(A)
  v <- seq.int(n)
  ord <- integer(n)
  roots <- which(!colSums(A))
  n_roots <- length(roots)
  j <- 1L
  while (n_roots > 0) {
    ord[seq.int(j, j + n_roots - 1L)] <- v[roots]
    v <- v[-roots]
    A <- A[-roots, -roots, drop = FALSE]
    j <- j + n_roots
    roots <- which(!colSums(A))
    n_roots <- length(roots)
  }
  if (nrow(A) > 0) {
    integer(0L)
  } else {
    ord
  }
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
#'
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

#' Unlist lapply
#'
#' @inheritParams base::lapply
#' @noRd
ulapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...))
}

#' Actual rank for an increasing sequence
#'
#' @param x An increasing `numeric` sequence
#' @noRd
rank_ <- function(x) {
  1L + c(0L, cumsum(diff(x) > 0))
}

#' Intersect Matrix Columns
#'
#' @param x List of matrices of identical dimensions
#' @noRd
matrix_intersect <- function(x) {
  nc <- ncol(x[[1L]])
  nr <- nrow(x[[1L]])
  out <- matrix(0L, nrow = nr, ncol = nc)
  for (i in seq_len(nc)) {
    tmp <-Reduce(intersect, lapply(x, function(y) y[,i]))
    out[seq_along(tmp), i] <- tmp
  }
  out
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
log_sum_exp_rows <- function(x) {
  maxs <- apply(x, 1L, max)
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

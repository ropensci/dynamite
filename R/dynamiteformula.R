#' Model formula for \pkg{dynamite}
#'
#' Defines a new observational or a new auxiliary channel for the model using
#' standard \R formula syntax. Formulas of individual response variables can be
#' joined together via `+`. See 'Details' and the package vignette for more
#' information. The function `obs` is a shorthand alias for `dynamiteformula`,
#' and `aux` is a shorthand alias for
#' `dynamiteformula(formula, family = "deterministic")`.
#'
#' Currently the \pkg{dynamite} package supports the following
#' distributions for the observations:
#'
#' * Categorical: `categorical` (with a softmax link using the first category
#'   as reference). See the documentation of the `categorical_logit_glm` in the
#'   Stan function reference manual (https://mc-stan.org/users/documentation/).
#' * Gaussian: `gaussian` (identity link, parameterized using mean and standard
#'   deviation).
#' * Poisson: `poisson` (log-link, with an optional known offset variable).
#' * Negative-binomial: `negbin` (log-link, using mean and dispersion
#'   parameterization, with an optional known offset variable). See the
#'   documentation on `NegBinomial2` in the Stan function reference manual.
#' * Bernoulli: `bernoulli` (logit-link).
#' * Binomial: `binomial` (logit-link).
#' * Exponential: `exponential` (log-link).
#' * Gamma: `gamma` (log-link, using mean and shape parameterization).
#' * Beta: `beta` (logit-link, using mean and precision parameterization).
#'
#' The models in the \pkg{dynamite} package are defined by combining the
#' channel-specific formulas defined via \R formula syntax.
#' Each channel is defined via the `obs` function, and the channels are
#' combined with `+`. For example a formula
#' `obs(y ~ lag(x), family = "gaussian") + obs(x ~ z, family = "poisson")`
#' defines a model with two channels;
#' first we declare that `y` is a gaussian variable depending on a previous
#' value of `x` (`lag(x)`), and then we add a second channel declaring `x` as
#' Poisson distributed depending on some exogenous variable `z`
#' (for which we do not define any distribution).
#'
#' In addition to declaring response variables via `obs`, we can also use
#' the function `aux` to define auxiliary channels which are deterministic
#' functions of other variables. The values of auxiliary variables are computed
#' dynamically during prediction, making the use of lagged values and other
#' transformations possible. Note that the auxiliary channel can also depend
#' on other variables without lags. The function `aux` also does not use the
#' `family` argument, which is automatically set to `deterministic` and is a
#' special channel type of `obs`. Note that lagged values of deterministic
#' `aux` channels do not imply fixed time points. Instead they must be given
#' starting values using a special function `init` that directly initializes
#' the lags to specified values, or by `past` which computes the initial values
#' based on an R expression. Both `init` and `past` should appear on the
#' right hand side of the model formula, separated from the primary defining
#' expression via `|`.
#'
#' The formula within `obs` can also contain an additional special
#' function `varying`, which defines the time-varying part of the model
#' equation, in which case we could write for example
#' `obs(x ~ z + varying(~ -1 + w), family = "poisson")`, which defines a model
#' equation with a constant intercept and time-invariant effect of `z`, and a
#' time-varying effect of `w`. We also remove the duplicate intercept with `-1`
#' in order to avoid identifiability issues in the model estimation
#' (we could also define a time varying intercept, in which case we would write
#' `obs(x ~ -1 + z + varying(~ w), family = "poisson)`). The part of the formula
#' not wrapped with `varying` is assumed to correspond to the fixed part of the
#' model, so `obs(x ~ z + varying(~ -1 + w), family = "poisson")` is actually
#' identical to
#' `obs(x ~ -1 + fixed(~ z) + varying(~ -1 + w), family = "poisson")` and
#' `obs(x ~ fixed(~ z) + varying(~ -1 + w), family = "poisson")`.
#'
#' When defining varying effects, we also need to define how the these
#' time-varying regression coefficient behave. For this, a `splines` component
#' should be added to the model, e.g.,
#' `obs(x ~ varying(~ -1 + w), family = "poisson) + splines(df = 10)` defines a
#' cubic B-spline with 10 degrees of freedom for the time-varying coefficient
#' corresponding to the `w`. If the model contains multiple time-varying
#' coefficients, same spline basis is used for all coefficients, with unique
#' spline coefficients and their standard deviation.
#'
#' If the desired model contains lagged predictors of each response in each
#' channel, these can be quickly added to the model as either time-invariant
#' or time-varying predictors via [lags()] instead of writing them manually
#' for each channel.
#'
#' It is also possible to define a random intercept term for each group by
#' using the component [random()] where the first argument defines for which
#' channels the intercept should be added, and second argument defines whether
#' or not these intercepts should be correlated between channels. This leads
#' to a model where in addition to the common intercept, each individual/group
#' has their own intercept with zero-mean normal prior and unknown standard
#' deviation (or multivariate gaussian in case `correlated = TRUE`),
#' analogously with the typical mixed models. Note however that if the channel
#' already contains the lagged response variable, the "intercept" is actually a
#' slope of (linear) trend as dynamite does not do any centering of variables.
#'
#' @export
#' @param formula \[`formula`]\cr An \R formula describing the model.
#' @param family \[`character(1)`]\cr The family name. See 'Details' for the
#'   supported families.
#' @return A `dynamiteformula` object.
#' @srrstats {G2.3b} Uses tolower.
#' @srrstats {RE1.0} Uses a formula interface.
#' @examples
#' # A single gaussian response channel with a time-varying effect of 'x',
#' # and a time-varying effect of the lag of 'y' using B-splines with
#' # 20 degrees of freedom for the coefficients of the time-varying terms.
#' obs(y ~ -1 + varying(~x), family = "gaussian") +
#'   lags(type = "varying") +
#'   splines(df = 20)
#'
#' # A two-channel categorical model with time-invariant predictors
#' # here, lag terms are specified manually
#' obs(x ~ z + lag(x) + lag(y), family = "categorical") +
#'   obs(y ~ z + lag(x) + lag(y), family = "categorical")
#'
#' # The same categorical model as above, but with the lag terms
#' # added using 'lags'
#' obs(x ~ z, family = "categorical") +
#'   obs(y ~ z, family = "categorical") +
#'   lags(type = "fixed")
#'
#' # A multichannel model with a gaussian, Poisson and a Bernoulli response and
#' # an auxiliary channel for the logarithm of 'p' plus one
#' obs(g ~ lag(g) + lag(logp), family = "gaussian") +
#'   obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
#'   obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
#'   aux(numeric(logp) ~ log(p + 1))
#'
dynamiteformula <- function(formula, family) {
  stopifnot_(
    !missing(formula),
    "Argument {.arg formula} is missing."
  )
  stopifnot_(
    !missing(family),
    "Argument {.arg family} is missing."
  )
  stopifnot_(
    inherits(formula, "formula"),
    "Argument {.arg formula} must be a {.cls formula} object."
  )
  stopifnot_(
    length(formula) == 3L,
    "Argument {.arg formula} must include
     a response and the model specification."
  )
  stopifnot_(
    checkmate::test_string(x = family, na.ok = FALSE),
    "Argument {.arg family} must be a single {.cls character} string."
  )
  family <- tolower(family)
  stopifnot_(
    is_supported(family),
    "Family {.val {family}} is not supported."
  )
  family <- do.call(paste0(family, "_"), args = list())
  stopifnot_(
    !"I" %in% all.names(formula),
    "{.code I(.)} is not supported by {.fun dynamiteformula}."
  )
  x <- dynamiteformula_(formula, formula, family)
  structure(
    list(
      dynamitechannel(
        formula = x$formula,
        original = formula,
        family = x$family,
        response = x$response,
        fixed = x$fixed,
        varying = x$varying,
        specials = x$specials,
        has_fixed_intercept = x$has_fixed_intercept,
        has_varying_intercept = x$has_varying_intercept
      )
    ),
    class = "dynamiteformula"
  )
}

#' Internal Version of `dynamiteformula`
#'
#' @inheritParams dynamiteformula
#' @param original The original `formula` definition.
#' @noRd
dynamiteformula_ <- function(formula, original, family) {
  if (is_deterministic(family)) {
    out <- formula_past(formula)
    resp_parsed <- deterministic_response(deparse1(formula_lhs(formula)))
    out$specials$resp_type <- resp_parsed$type
    out$response <- resp_parsed$resp
  } else {
    out <- formula_specials(formula)
    if (is_binomial(family)) {
      stopifnot_(
        "trials" %in% names(out$specials),
        "Formula for a binomial channel must include a trials term."
      )
    }
    out$response <- deparse1(formula_lhs(formula))
  }
  out$family <- family
  out$original <- original
  out
}

#' Create a Channel For a `dynamiteformula` Object Directly
#'
#' @inheritParams dynamiteformula
#' @param response \[`character(1)`]\cr Name of the response.
#' @param fixed \[`integer()`]\cr Time-invariant covariate indices.
#' @param varying \[`integer()`]\cr Time-varying covariate indices.
#' @param has_fixed_intercept \[`logical(1)`]\cr Does the channel contain fixed
#'   intercept?
#' @param has_varying_intercept \[`logical(1)`]\cr Does the channel contain
#'   varying intercept?
#' @param has_rand_intercept \[`logical(1)`]\cr Does the channel contain random
#'   individual-level intercept term?
#' @noRd
dynamitechannel <- function(formula, original = NULL, family, response,
                            fixed = integer(0L), varying = integer(0L),
                            specials = list(),
                            has_fixed_intercept = FALSE,
                            has_varying_intercept = FALSE) {
  list(
    formula = formula,
    original = original,
    family = family,
    response = response,
    fixed = fixed,
    varying = varying,
    specials = specials,
    has_fixed_intercept = has_fixed_intercept,
    has_varying_intercept = has_varying_intercept
  )
}

#' @rdname dynamiteformula
#' @export
obs <- dynamiteformula

#' Is the Argument a `dynamiteformula` Object?
#'
#' @param x An \R object.
#' @noRd
is.dynamiteformula <- function(x) {
  inherits(x, "dynamiteformula")
}

#' @rdname dynamiteformula
#' @export
aux <- function(formula) {
  dynamiteformula(formula, family = "deterministic")
}

#' @rdname dynamiteformula
#' @param e1 \[`dynamiteformula`]\cr A model formula specification.
#' @param e2 \[`dynamiteformula`]\cr A model formula specification.
#' @export
#' @examples
#' obs(y ~ x, family = "gaussian") + obs(z ~ w, family = "exponential")
#'
`+.dynamiteformula` <- function(e1, e2) {
  stopifnot_(
    is.dynamiteformula(e1),
    "Method {.fun +.dynamiteformula} is not supported for
     {.cls {class(e1)}} objects."
  )
  add_dynamiteformula(e1, e2)
}

#' @rdname dynamiteformula
#' @param x \[`dynamiteformula`]\cr The model formula.
#' @param ... Ignored.
#' @export
#' @examples
#' x <- obs(y ~ x, family = "gaussian") +
#'   obs(z ~ w, family = "exponential") +
#'   aux(d ~ log(y) | init(c(0, 1))) +
#'   lags(k = 2) +
#'   splines(df = 5) +
#'   random(responses = c("y", "z"), correlated = TRUE)
#' print(x)
#'
print.dynamiteformula <- function(x, ...) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamiteformula(x),
    "Argument {.arg x} must be a {.cls dynamiteformula} object."
  )
  out <- data.frame(
    Family = get_family_names(x),
    Formula = vapply(get_originals(x), function(y) deparse1(y), character(1L))
  )
  rownames(out) <- get_responses(x)
  print.data.frame(out, right = FALSE)
  if (!is.null(attr(x, "lags"))) {
    k <- attr(x, "lags")$k
    cat("\nLagged responses added as predictors with: k = ", cs(k), sep = "")
  }
  if (!is.null(attr(x, "random"))) {
    resp <- attr(x, "random")$responses
    co <- attr(x, "random")$correlated
    cat(
      ifelse_(co, "\nCorrelated random ", "\nRandom "),
      "intercepts added for response(s): ",
      cs(resp),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' Get All Response Variables of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_responses <- function(x) {
  vapply(x, function(y) y$response, character(1L))
}

#' Get The RHS of All Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_predictors <- function(x) {
  vapply(x, function(y) deparse1(formula_rhs(y$formula)), character(1L))
}

#' Get Terms of All Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_terms <- function(x) {
  lapply(x, function(y) {
    if (is_deterministic(y$family)) {
      character(0L)
    } else {
      attr(terms(y$formula), "term.labels")
    }
  })
}

#' Get All Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_formulas <- function(x) {
  lapply(x, "[[", "formula")
}

#' Get All Original Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_originals <- function(x) {
  lapply(x, "[[", "original")
}

#' Get All Family Objects of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_families <- function(x) {
  lapply(x, "[[", "family")
}

#' Get All Family Names of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_family_names <- function(x) {
  vapply(x, function(x) x$family$name, character(1L))
}

#' Get a Quoted Expression of Deterministic Channel Definitions
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_quoted <- function(x) {
  resp <- get_responses(x)
  if (length(resp) > 0L) {
    expr <- lapply(x, function(x) deparse1(formula_rhs(x$formula)))
    quote_str <- paste0(
      "`:=`(",
      paste0(resp, " = ", expr, collapse = ","),
      ")"
    )
    str2lang(quote_str)
  } else {
    NULL
  }
}

#' Get Indices of Deterministic Channels in a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
which_deterministic <- function(x) {
  which(vapply(x, function(y) is_deterministic(y$family), logical(1L)))
}

#' Get Indices of Stochastic Channels in a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object
#' @noRd
which_stochastic <- function(x) {
  which(vapply(x, function(y) !is_deterministic(y$family), logical(1L)))
}

#' Get Channels with Past Value Definitions of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object
#' @noRd
has_past <- function(x) {
  vapply(x, function(y) !is.null(y$specials$past), logical(1L))
}

#' Internal `+.dynamiteformula` For Model Construction
#'
#' @param e1 An \R object.
#' @param e2 An \R object.
#' @noRd
add_dynamiteformula <- function(e1, e2) {
  if (is.dynamiteformula(e2)) {
    out <- join_dynamiteformulas(e1, e2)
  } else if (inherits(e2, "lags")) {
    out <- set_lags(e1, e2)
  } else if (inherits(e2, "splines")) {
    out <- set_splines(e1, e2)
  } else if (inherits(e2, "random")) {
    out <- set_random(e1, e2)
  } else {
    stop_(
      "Unable to add an object of class {.cls {class(e2)}}
      to an object of class {.cls dynamiteformula}."
    )
  }
  out
}

#' Join Two Model Definitions and Verify Compatibility
#'
#' @param e1 A `dynamiteformula` object.
#' @param e2 A `dynamiteformula` object.
#' @noRd
join_dynamiteformulas <- function(e1, e2) {
  out <- c(e1, e2)
  resp_list <- list(
    get_responses(e1),
    get_responses(e2)
  )
  resp_all <- unlist(resp_list)
  resp_duped <- duplicated(resp_all)
  stopifnot_(
    !any(resp_duped),
    "Multiple definitions for response variable{?s}
     {.var {resp_all[resp_duped]}}."
  )
  stopifnot_(
    is.null(attr(e1, "lags")) || is.null(attr(e2, "lags")),
    "Both dynamiteformulas contain a lags definition."
  )
  stopifnot_(
    is.null(attr(e1, "splines")) || is.null(attr(e2, "splines")),
    "Both dynamiteformulas contain a splines definition."
  )
  stopifnot_(
    is.null(attr(e1, "random")) || is.null(attr(e2, "random")),
    "Both dynamiteformulas contain a random intercepts definition."
  )
  rhs_list <- list(
    lapply(get_terms(e1), extract_nonlags),
    lapply(get_terms(e2), extract_nonlags)
  )
  stoch_list <- list(
    which_stochastic(e1),
    which_stochastic(e2)
  )
  for (i in 1L:2L) {
    resp_a <- resp_list[[i]]
    resp_b <- resp_list[[3L - i]][stoch_list[[3L - i]]]
    rhs <- rhs_list[[3L - i]][stoch_list[[3L - i]]]
    if (length(rhs) > 0L) {
      for (j in seq_along(resp_a)) {
        simul_lhs <- resp_a[j]
        simul <- vapply(rhs, function(x) simul_lhs %in% x, logical(1L))
        stopifnot_(
          !any(simul),
          c(
            "Simultaneous regression is not supported:",
            `x` = "Response variable {.var {simul_lhs}} appears in
                  the formula of {.var {resp_b[which(simul)[1L]]}}."
          )
        )
      }
    }
  }
  attributes(out) <- c(attributes(e1), attributes(e2))
  class(out) <- "dynamiteformula"
  out
}

#' Set Lags Definitions for All Channels in a `dynamiteformula` Object
#'
#' @param e1 A `dynamiteformula` object.
#' @param e2 A `lags` object.
#' @noRd
set_lags <- function(e1, e2) {
  stopifnot_(
    is.null(attr(e1, "lags")),
    "Multiple definitions for lags."
  )
  attr(e1, "lags") <- e2
  e1
}

#' Set the Regression Coefficient Splines of the Model
#'
#' @param e1 A `dynamiteformula` object.
#' @param e2 A `splines` object.
#' @noRd
set_splines <- function(e1, e2) {
  stopifnot_(
    is.null(attr(e1, "splines")) || attr(e2, "override"),
    "Multiple definitions for splines."
  )
  attr(e1, "splines") <- e2
  e1
}

#' Set the Random Intercepts of the Model
#'
#' @param e1 A `dynamiteformula` object.
#' @param e2 A `random` object.
#' @noRd
set_random <- function(e1, e2) {
  stopifnot_(
    is.null(attr(e1, "random")) || attr(e2, "random"),
    "Multiple definitions for random intercepts."
  )
  attr(e1, "random") <- e2
  e1
}

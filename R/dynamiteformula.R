#' Model formula for \pkg{dynamite}
#'
#' Defines a new observational or a new auxiliary channel for the model.
#' See 'Details' or the package vignette for more information.
#'
#' @details Currently the `dynamite` package supports the following
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
#' * Gamma: `gamma` (log-link, using mean and shape parameterisation).
#' * Beta: `beta` (logit-link, using mean and precision parameterisation).
#'
#' The models in the \pkg{dynamite} package are defined by combining the
#' channel-specific formulas defined via  \R formula syntax.
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
#' special channel type of `obs`.
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
#' It is also possible to define a random intercept term for each group by
#' using component `random` where the first argument defines for which channels
#' the intercept should be added, and second argument defines whether or not
#' these intercepts should be correlated between channels. This leads to a
#' model where the in addition to the common intercept each individual/group
#' has their own intercept with zero-mean normal prior and unknown standard
#' deviation (or multivariate gaussian in case `correlated = TRUE`),
#' analogously with the typical mixed models. Note however that with
#' a large number of time points these intercepts can become challenging
#' sample with default priors. This is because with large group sizes the
#' group-level intercepts tend to be behave similarly to fixed group-factor
#' variable so the model becomes overparameterized given these and the common
#' intercept term. In these cases, a better option might be to use fixed group
#' effects, i.e. to include the grouping variable to the model instead of
#' random intercept (this way the one group is included in the intercept).
#'
#' @param formula \[`formula`]\cr An \R formula describing the model.
#' @param family \[`character(1)`]\cr The family name. See 'Details' for the
#' supported families.
#' @return An object of class `dynamiteformula`.
#' @export
#' @examples
#' obs(y ~ -1 + varying(~x), family = "gaussian") +
#'   lags(type = "varying") + splines(df = 20)
#'
#' @srrstats {G2.3b} Uses tolower.
#' @srrstats {RE1.0} Uses a formula interface.
dynamiteformula <- function(formula, family) {
  stopifnot_(
    is.formula(formula),
    "Argument {.arg formula} must be a {.cls formula} object."
  )
  stopifnot_(
    checkmate::test_string(x = family, na.ok = FALSE),
    "Argument {.arg family} must be a single {.cls character} string."
  )
  family <- tolower(family)
  if (is_supported(family)) {
    family <- do.call(paste0(family, "_"), args = list())
  } else {
    stop_("Family {.val {family}} is not supported.")
  }
  stopifnot_(
    !has_as_is(deparse1(formula)),
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

#' @describeIn dynamiteformula Internal Version of `dynamiteformula`
#' @noRd
dynamiteformula_ <- function(formula, original, family) {
  if (is_deterministic(family)) {
    out <- formula_past(formula)
    resp_parsed <- formula_response(deparse1(formula_lhs(formula)))
    out$specials$resp_type <- resp_parsed$type
    if (!is.null(out$specials$past)) {
      out$specials$past <- do.call(
        what = paste0("as.", resp_parsed$type),
        args = list(out$specials$past)
      )
    }
    out$response <- resp_parsed$resp
  } else {
    out <- formula_specials(formula)
    out$response <- deparse1(formula_lhs(formula))
  }
  out$family <- family
  out$original <- original
  out
}

#' Create a Channel for a `dynamiteformula` Directly
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
#'@noRd
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

#' Is the Argument a `dynamiteformula` Object
#'
#' @param x An \R object.
#' @noRd
is.dynamiteformula <- function(x) {
  inherits(x, "dynamiteformula")
}

#' @describeIn dynamiteformula Prepare a Deterministic Auxiliary Channel
#' @export
aux <- function(formula) {
  dynamiteformula(formula, family = "deterministic")
}

#' Join Two `dynamiteformula` Objects
#'
#' @param e1 An \R object.
#' @param e2 An \R object.
#' @export
#' @examples
#' obs(y ~ x, family = "gaussian") + obs(z ~ w, family = "exponential")
`+.dynamiteformula` <- function(e1, e2) {
  if (is.dynamiteformula(e1)) {
    out <- add_dynamiteformula(e1, e2)
  } else {
    stop_("Method {.fun +.dynamiteformula} is not supported for
          {.cls {class(e1)}} objects.")
  }
  out
}

#' Is the Argument a `formula` Object.
#'
#' @param x An \R object
#' @noRd
is.formula <- function(x) {
  inherits(x, "formula")
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
  vapply(x, function(y) length(y$specials$past) > 0L, logical(1L))
}

#' Internal `+.dynamiteformula` For Model Construction
#'
#' @param e1 An \R object.
#' @param e2 An \R object.
#' @noRd
add_dynamiteformula <- function(e1, e2) {
  if (is.dynamiteformula(e2)) {
    out <- join_dynamiteformulas(e1, e2)
  } else if (is.lags(e2)) {
    out <- set_lags(e1, e2)
  } else if (is.splines(e2)) {
    out <- set_splines(e1, e2)
  } else if (is.random(e2)) {
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
#' @param e1 A `dynamiteformula` object,
#' @param e2 A `lags` object,
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

#' Model Formula for \pkg{dynamite}
#'
#' Defines a new observational or a new auxiliary channel for the model using
#' standard \R formula syntax. Formulas of individual response variables can be
#' joined together via `+`. See 'Details' and the package vignettes for more
#' information. The function `obs` is a shorthand alias for `dynamiteformula`,
#' and `aux` is a shorthand alias for
#' `dynamiteformula(formula, family = "deterministic")`.
#'
#' Currently the \pkg{dynamite} package supports the following
#' distributions for the observations:
#'
#' * Categorical: `categorical` (with a softmax link using the first category
#'   as reference). See the documentation of the `categorical_logit_glm` in the
#'   Stan function reference manual <https://mc-stan.org/users/documentation/>.
#' * Multinomial: `multinomial` (softmax link, first category is reference).
#' * Gaussian: `gaussian` (identity link, parameterized using mean and standard
#'   deviation).
#' * Multivariate Gaussian: `mvgaussian` (identity link, parameterized using
#'   mean vector, standard deviation vector and the Cholesky decomposition of
#'   the correlation matrix).
#' * Poisson: `poisson` (log-link, with an optional known offset variable).
#' * Negative-binomial: `negbin` (log-link, using mean and dispersion
#'   parameterization, with an optional known offset variable). See the
#'   documentation on `NegBinomial2` in the Stan function reference manual.
#' * Bernoulli: `bernoulli` (logit-link).
#' * Binomial: `binomial` (logit-link).
#' * Exponential: `exponential` (log-link).
#' * Gamma: `gamma` (log-link, using mean and shape parameterization).
#' * Beta: `beta` (logit-link, using mean and precision parameterization).
#' * Student t: `student` (identity link, parameterized using degrees of
#'   freedom, location and scale)
#'
#' The models in the \pkg{dynamite} package are defined by combining the
#' channel-specific formulas defined via \R formula syntax.
#' Each channel is defined via the `obs` function, and the channels are
#' combined with `+`. For example a formula
#' `obs(y ~ lag(x), family = "gaussian") + obs(x ~ z, family = "poisson")`
#' defines a model with two channels;
#' first we declare that `y` is a Gaussian variable depending on a previous
#' value of `x` (`lag(x)`), and then we add a second channel declaring `x` as
#' Poisson distributed depending on some exogenous variable `z`
#' (for which we do not define any distribution).
#'
#' Number of trials for binomial channels should be defined via a `trials`
#' model component, e.g., `obs(y ~ x + trials(n), family = "binomial")`,
#' where `n` is a data variable defining the number of trials. For multinomial
#' channels, the number of trials is automatically defined to be the sum
#' of the observations over the categories, but can also be defined using
#' the `trials` component, for example for prediction.
#'
#' Multivariate channels are defined by providing a single formula for all
#' components or by providing component-specific formulas separated by a `|`.
#' The response variables that correspond to the components should be joined by
#' `c()`. For instance, the following would define `c(y1, y2)` as multivariate
#' gaussian with `x` as a predictor for the mean of the first component and
#' `x` and `z` as predictors for the mean of the second component:
#' `obs(c(y1, y2) ~ x | x + z, family = "mvgaussian")`. A multinomial channel
#' should only have a single formula.
#'
#' In addition to declaring response variables via `obs`, we can also use
#' the function `aux` to define auxiliary channels which are deterministic
#' functions of other variables. The values of auxiliary variables are computed
#' dynamically during prediction, making the use of lagged values and other
#' transformations possible. The function `aux` also does not use the
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
#' It is also possible to define group-specific (random) effects term
#' using the special syntax `random()` similarly as `varying()`. For example,
#' `random(~1)` leads to a model where in addition to the common intercept,
#' each individual/group has their own intercept with zero-mean normal prior and
#' unknown standard deviation analogously with the typical mixed models. An
#' additional model component [random_spec()] can be used to define
#' whether the random effects are allowed to correlate within and across
#' channels and whether to use centered or noncentered parameterization for
#' the random effects.
#'
#' @export
#' @family formulas
#' @param formula \[`formula`]\cr An \R formula describing the model.
#' @param family \[`character(1)`]\cr The family name. See 'Details' for the
#'   supported families.
#' @param link \[`character(1)`]\cr The name of the link function to use or
#'   `NULL`. See details for supported link functions and default values of
#'   specific families.
#' @return A `dynamiteformula` object.
#' @srrstats {G2.3b} Uses tolower.
#' @srrstats {RE1.0} Uses a formula interface.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
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
dynamiteformula <- function(formula, family, link = NULL) {
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
    is_supported_family(family),
    "Family {.val {family}} is not supported."
  )
  family <- do.call(paste0(family, "_"), args = list(link = link))
  stopifnot_(
    !"I" %in% all.names(formula),
    "{.code I(.)} is not supported by {.fun dynamiteformula}."
  )
  if (is_deterministic(family)) {
    dims <- list(formula_past(formula))
    resp_parsed <- deterministic_response(deparse1(formula_lhs(formula)))
    dims[[1L]]$specials$resp_type <- resp_parsed$type
    dims[[1L]]$response <- resp_parsed$resp
    dims[[1L]]$original <- formula
    dims[[1L]]$name <- stan_name(resp_parsed$resp)
  } else {
    dims <- parse_formula(formula, family)
    if (is_binomial(family) || is_multinomial(family)) {
      stopifnot_(
        "trials" %in% names(dims[[1L]]$specials),
        "Formula for a {family$name} channel must include a trials term."
      )
    }
  }
  structure(
    lapply(dims, function(x) {
      dynamitechannel(
        formula = x$formula,
        original = x$original,
        family = family,
        response = x$response,
        name = x$name,
        fixed = x$fixed,
        varying = x$varying,
        random = x$random,
        specials = x$specials,
        has_fixed_intercept = x$has_fixed_intercept,
        has_varying_intercept = x$has_varying_intercept,
        has_random_intercept = x$has_random_intercept
      )
    }),
    class = "dynamiteformula",
    channel_groups = rep(1L, length(dims)),
    model_topology = 1L
  )
}

#' Create a Channel For a `dynamiteformula` Object Directly
#'
#' @inheritParams dynamiteformula
#' @param response \[`character(1)`]\cr Name of the response.
#' @param fixed \[`integer()`]\cr Time-invariant covariate indices.
#' @param varying \[`integer()`]\cr Time-varying covariate indices.
#' @param random \[`integer()`]\cr Random effect covariate indices.
#' @param has_fixed_intercept \[`logical(1)`]\cr Does the channel contain
#'   a fixed intercept?
#' @param has_varying_intercept \[`logical(1)`]\cr Does the channel contain
#'   a varying intercept?
#' @param has_random_intercept \[`logical(1)`]\cr Does the channel contain
#'   a random group-level intercept term?
#' @noRd
dynamitechannel <- function(formula, original = NULL,
                            family, response, name = NULL,
                            fixed = integer(0L), varying = integer(0L),
                            random = integer(0L), specials = list(),
                            has_fixed_intercept = FALSE,
                            has_varying_intercept = FALSE,
                            has_random_intercept = FALSE) {
  list(
    formula = formula,
    original = original,
    family = family,
    response = response,
    name = name,
    fixed = fixed,
    varying = varying,
    random = random,
    specials = specials,
    has_fixed_intercept = has_fixed_intercept,
    has_varying_intercept = has_varying_intercept,
    has_random_intercept = has_random_intercept
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

#' Parse Channel Formulas for `dynamiteformula`
#'
#' @param x A `formula` object.
#' @noRd
parse_formula <- function(x, family) {
  responses <- all.vars(formula_lhs(x))
  formula_str <- deparse1(formula_rhs(x))
  formula_parts <- strsplit(formula_str, "|", fixed = TRUE)[[1L]]
  n_formulas <- length(formula_parts)
  n_responses <- length(responses)
  if (is_multivariate(family)) {
    stopifnot_(
      n_responses > 1L,
      "A multivariate channel must have more than one response variable."
    )
    if (is_multinomial(family)) {
      stopifnot_(
        n_formulas == 1L,
        "A multinomial channel must have only one formula component."
      )
    } else {
      stopifnot_(
        n_formulas == n_responses || n_formulas == 1L,
        c(
          "Number of component formulas must be 1 or
           the number of dimensions: {n_responses}",
          `x` = "{n_formulas} formulas were provided."
        )
      )
    }
  } else {
    stopifnot_(
      n_responses == 1L,
      "A univariate channel must have only one response variable."
    )
    stopifnot_(
      n_formulas == 1L,
      "A univariate channel must have only one formula component."
    )
  }
  formula_parts <- ifelse_(
    n_formulas == 1L,
    rep(formula_parts, n_responses),
    formula_parts
  )
  responses <- str_quote(responses)
  formulas <- lapply(paste0(responses, "~", formula_parts), as.formula)
  predictors <- lapply(
    formulas,
    function(y) {
      find_nonlags(formula_rhs(y))
    }
  )
  resp_pred <- vapply(
    seq_along(responses),
    function(i) {
      responses[i] %in% predictors[[i]]
    },
    logical(1L)
  )
  p <- sum(resp_pred)
  stopifnot_(
    !any(resp_pred),
    c(
      "Contemporaneous self-dependency found in model formula:",
      `x` = "{cli::qty(p)} Variable{?s} {.arg {cs(responses[resp_pred])}}
             appear{?s/} on both sides of the formula for ({cs(responses)})."
    )
  )
  lapply(
    formulas,
    formula_specials,
    original = x,
    family = family
  )
}

#' @rdname dynamiteformula
#' @param e1 \[`dynamiteformula`]\cr A model formula specification.
#' @param e2 \[`dynamiteformula`]\cr A model formula specification.
#' @export
#' @examples
#' data.table::setDTthreads(1) # For CRAN
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
#' data.table::setDTthreads(1) # For CRAN
#' x <- obs(y ~ x + random(~ 1 + lag(d)), family = "gaussian") +
#'   obs(z ~ varying(~w), family = "exponential") +
#'   aux(numeric(d) ~ log(y) | init(c(0, 1))) +
#'   lags(k = 2) +
#'   splines(df = 5) +
#'   random_spec(correlated = FALSE)
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
  cg <- attr(x, "channel_groups")
  n_cg <- n_unique(cg)
  rn <- character(n_cg)
  out <- data.frame(
    Family = rep(NA_character_, n_cg),
    Formula = rep(NA_character_, n_cg)
  )
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    j <- cg_idx[1L]
    rn[i] <- ifelse_(
      is_multivariate(x[[j]]$family),
      paste(get_names(x[cg_idx]), collapse = "_"),
      x[[j]]$name
    )
    out[i, "Family"] <- x[[j]]$family$name
    out[i, "Formula"] <- deparse1(x[[j]]$original)
  }
  rownames(out) <- rn
  print.data.frame(out, right = FALSE)
  if (!is.null(attr(x, "lags"))) {
    k <- attr(x, "lags")$k
    type <- attr(x, "lags")$type
    cat(
      "\nLagged responses added as ", type, " predictors with: k = ", cs(k),
      sep = ""
    )
  }
  if (!is.null(attr(x, "random_spec"))) {
    rand <- which_random(x)
    co <- attr(x, "random_spec")$correlated
    cat(
      ifelse_(co, "\nCorrelated random ", "\nRandom "),
      "effects added for response(s): ",
      cs(get_names(x)[rand]),
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
get_rhs <- function(x) {
  vapply(x, function(y) deparse1(formula_rhs(y$formula)), character(1L))
}

#' Get the Names of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_names <- function(x) {
  vapply(x, function(y) y$name, character(1L))
}

#' Get Lagged Terms of All Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_lag_terms <- function(x) {
  lapply(x, function(y) {
    unique(find_lags(formula_rhs(y$formula)))
  })
}

#' Get Non-Lagged Terms of All Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_nonlag_terms <- function(x) {
  lapply(x, function(y) {
    unique(find_nonlags(formula_rhs(y$original)))
  })
}

#' Get the Order of All Lag Terms of All Formulas of a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object.
#' @noRd
get_lag_orders <- function(x) {
  tmp <- lapply(x, function(y) {
    tmp_ <- find_lag_orders(formula_rhs(y$original))
    if (nrow(tmp_) > 0L) {
      tmp_$resp <- y$response
    }
    tmp_
  })
  unique(rbindlist_(tmp[vapply(tmp, nrow, integer(1L)) > 0]))
}

#' Get Special Type Formula of a Dimension in a `dynamiteformula`
#'
#' @param x A channel of a `dynamiteformula`
#' @noRd
get_type_formula <- function(x, type = c("fixed", "varying", "random")) {
  has_icpt <- ifelse_(
    type %in% c("fixed", "varying"),
    x$has_fixed_intercept || x$has_varying_intercept,
    x$has_random_intercept
  )
  icpt <- ifelse_(has_icpt, "1", "-1")
  idx <- x[[type]]
  ft <- terms(x$formula)
  resp <- attr(ft, "variables")[[2L]]
  tr <- attr(ft, "term.labels")
  rhs <- paste0(tr[idx], collapse = " + ")
  rhs_out <- ifelse_(nzchar(rhs), paste0(" + ", rhs), "")
  out_str <- paste0(str_quote(deparse1(resp)), " ~ ", icpt, rhs_out)
  ifelse_(
    has_icpt || nzchar(rhs_out),
    as.formula(out_str),
    NULL
  )
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

#' Get a Directed Acyclic Graph (DAG) of a `dynamiteformula` Object.
#'
#' @param x A `dynamiteformula` object.
#' @param project A `logical` value. If `TRUE`, deterministic responses are
#'   projected out of the DAG.
#' @param covariates A `logical` value. If `TRUE`, fixed covariates are also
#'   included in the DAG.
#' @param format A `character` string describing how to label variable names.
#' @noRd
get_dag <- function(x, project = FALSE, covariates = FALSE,
                    format = c("default", "expression", "lag")) {
  f_ <- dag_formatter(format)
  resp <- get_responses(x)
  contemp_dep <- ifelse_(
    covariates,
    get_nonlag_terms(x),
    lapply(get_nonlag_terms(x), function(y) y[y %in% resp])
  )
  lag_dep <- get_lag_orders(x)
  cg <- attr(x, "channel_groups")
  lag_dep <- ifelse_(
    covariates,
    lag_dep,
    lag_dep[lag_dep$var %in% resp, ]
  )
  if (project) {
    resp_det <- which_deterministic(x)
    for (i in resp_det) {
      contemp_pa <- contemp_dep[[i]]
      k <- length(contemp_pa)
      if (k > 0L) {
        contemp_dep[-i] <- lapply(
          contemp_dep[-i],
          function(y) {
            ifelse_(
              contemp_pa %in% y,
              union(y, contemp_pa),
              y
            )
          }
        )
      }
      lag_dep_pa <- lag_dep[lag_dep$resp == resp[i], ]
      lag_dep_ch <- lag_dep[lag_dep$var == resp[i], ]
      lag_dep_new <- vector(mode = "list", length = nrow(lag_dep_ch))
      if (nrow(lag_dep_pa) > 0L || k > 0L) {
        for (j in seq_len(nrow(lag_dep_ch))) {
          lag_dep_new[[j]] <- data.frame(
            var = c(contemp_pa, lag_dep_pa$var),
            order = c(rep(0L, k), lag_dep_pa$order) + lag_dep_ch$order[j],
            resp = lag_dep_ch$resp[j]
          )
        }
      }
      lag_dep <- rbind(
        lag_dep[lag_dep$resp != resp[i] & lag_dep$var != resp[i], ],
        rbindlist_(lag_dep_new)
      )
    }
    resp_stoch <- which_stochastic(x)
    contemp_dep <- contemp_dep[resp_stoch]
    resp <- resp[resp_stoch]
    cg <- cg[resp_stoch]
  }
  max_lag <- max(1L, max(lag_dep$order))
  all_vars <- c(
    resp,
    unique(setdiff(union(unlist(contemp_dep), lag_dep$var), resp))
  )
  resp_lag <- expand.grid(var = all_vars, order = seq_len(max_lag))
  v <- c(
    f_(resp_lag$var, -1L * resp_lag$order),
    f_(resp_lag$var, resp_lag$order),
    f_(all_vars, 0L)
  )
  n <- length(v)
  m <- nrow(lag_dep)
  A <- matrix(
    0L,
    nrow = n,
    ncol = n,
    dimnames = replicate(2L, v, simplify = FALSE)
  )
  if (covariates) {
    p <- length(resp)
    resp <- all_vars
    layout_y <- c(
      rev(seq_along(resp[order(cg)])),
      p + seq_len(length(all_vars) - p)
    )
  } else {
    layout_y <- rev(seq_along(resp[order(cg)]))
  }
  resp_past <- f_(resp, -1L)
  resp_future <- f_(resp, 1L)
  resp_t <- f_(resp, 0L)
  layout <- data.frame(var = v, x = NA, y = NA)
  layout[layout$var %in% resp_t, "x"] <- 0.0
  layout[layout$var %in% resp_t, "y"] <- layout_y
  layout[layout$var %in% resp_past, "x"] <- -1.0
  layout[layout$var %in% resp_past, "y"] <- layout_y
  layout[layout$var %in% resp_future, "x"] <- 1.0
  layout[layout$var %in% resp_future, "y"] <- layout_y
  var_past <- f_(lag_dep$var, -1 * lag_dep$order)
  resp_future <- f_(lag_dep$resp, lag_dep$order)
  resp_t <- f_(lag_dep$resp, 0L)
  var_t <- f_(lag_dep$var, 0L)
  A[cbind(var_past, resp_t)] <- 1L
  A[cbind(var_t, resp_future)] <- 1L
  edgelist <- list()
  edgelist[[1L]] <- data.frame(from = var_past, to = resp_t)
  edgelist[[2L]] <- data.frame(from = var_t, to = resp_future)
  e_idx <- 3L
  for (i in seq_len(max_lag - 1L)) {
    include <- (lag_dep$order + i) <= max_lag
    var_past <- f_(lag_dep$var, -1 * (lag_dep$order + i))
    var_past <- var_past[include]
    resp_past <- f_(lag_dep$resp, -1 * i)
    resp_past <- resp_past[include]
    var_future <- f_(lag_dep$var, i)
    var_future <- var_future[include]
    resp_future <- f_(lag_dep$resp, (lag_dep$order + i))
    resp_future <- resp_future[include]
    A[cbind(var_past, resp_past)] <- 1L
    A[cbind(var_future, resp_future)] <- 1L
    edgelist[[e_idx]] <- data.frame(from = var_past, to = resp_past)
    edgelist[[e_idx + 1L]] <- data.frame(from = var_future, to = resp_future)
    resp_past <- f_(resp, -1L * (i + 1L))
    resp_future <- f_(resp, i + 1L)
    layout[layout$var %in% resp_past, "x"] <- (-1.0) * (i + 1L)
    layout[layout$var %in% resp_past, "y"] <- layout_y
    layout[layout$var %in% resp_future, "x"] <- (1.0) * (i + 1L)
    layout[layout$var %in% resp_future, "y"] <- layout_y
    e_idx <- e_idx + 2L
  }
  for (i in seq_along(contemp_dep)) {
    k <- length(contemp_dep[[i]])
    if (k > 0L) {
      resp_ti <- f_(resp[i], 0L)
      contemp_t <- f_(contemp_dep[[i]], 0L)
      A[contemp_t, resp_ti] <- 1L
      edgelist[[e_idx]] <- data.frame(from = contemp_t, to = resp_ti)
      e_idx <- e_idx + 1L
      for (j in seq_len(max_lag)) {
        contemp_past <- f_(contemp_dep[[i]], -1L * j)
        contemp_future <- f_(contemp_dep[[i]], j)
        resp_past <- f_(resp[i], -1L * j)
        resp_future <- f_(resp[i], j)
        A[contemp_past, resp_past] <- 1L
        A[contemp_future, resp_future] <- 1L
        edgelist[[e_idx]] <-
          data.frame(from = contemp_past, to = resp_past)
        edgelist[[e_idx + 1L]] <-
          data.frame(from = contemp_future, to = resp_future)
        e_idx <- e_idx + 2L
      }
    }
  }
  list(A = A, edgelist = rbindlist_(edgelist), layout = layout)
}

#' Get the Markov Blanket of a Response Variable
#'
#' @param g Output of `get_dag()`.
#' @param y A `character` string naming the response variable.
#' @noRd
get_markov_blanket <- function(g, y) {
  A <- g$A
  v <- colnames(A)
  pa <- v[which(A[, y] == 1L)]
  ch <- v[which(A[y, ] == 1L)]
  pa_ch <- onlyif(
    length(ch) > 0L,
    unique(v[which(A[, ch, drop = FALSE] == 1L, arr.ind = TRUE)[, "row"]])
  )
  setdiff(union(c(ch, pa), pa_ch), y)
}

#' Create a Function to Format DAG Vertex Names
#'
#' @inheritParams get_dag
#' @noRd
dag_formatter <- function(format) {
  switch(
    format,
    default = function(v, k) {
      suffix <- c(" - ", "", " + ")
      k_str <- as.character(abs(k))
      k_str[k == 0L] <- ""
      s <- sign(k) + 2L
      paste0(v, "_{t", suffix[s], k_str, "}")
    },
    expression = function(v, k) {
      suffix <- c(" - ", "", " + ")
      k_str <- as.character(abs(k))
      k_str[k == 0L] <- ""
      s <- sign(k) + 2L
      paste0(v, "[t", suffix[s], k_str, "]")
    },
    lag = function(v, k) {
      suffix <- c("_lag", "", "_lead")
      k_str <- as.character(abs(k))
      k_str[k == 0L] <- ""
      s <- sign(k) + 2L
      paste0(v, suffix[s], k_str)
    }
  )
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

#' Get Responses with Random Effects in a `dynamiteformula` Object
#'
#' @param x A `dynamiteformula` object
#' @noRd
which_random <- function(x) {
  which(
    vapply(
      x,
      function(y) {
        !is_deterministic(y$family) &&
          (y$has_random_intercept || length(y$random) > 0L)
      },
      logical(1L)
    )
  )
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
  } else if (inherits(e2, "random_spec")) {
    out <- set_random_spec(e1, e2)
  } else if (inherits(e2, "latent_factor")) {
    out <- set_lfactor(e1, e2)
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
    is.null(attr(e1, "random_spec")) || is.null(attr(e2, "random_spec")),
    "Both dynamiteformulas contain a random_spec definition."
  )
  pred <- c(get_nonlag_terms(e1), get_nonlag_terms(e2))
  cg1 <- attr(e1, "channel_groups")
  cg2 <- attr(e2, "channel_groups")
  cg <-  c(cg1, cg2 + max(cg1))
  n_cg <- n_unique(cg)
  dep <- matrix(0L, nrow = n_cg, ncol = n_cg)
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    for (j in cg_idx) {
      dep[unique(cg[which(resp_all %in% pred[[j]])]), i] <- 1L
    }
  }
  topo <- topological_order(dep)
  stopifnot_(
    length(topo) > 0L,
    "The model must be acyclic."
  )
  attributes(out) <- c(attributes(e1), attributes(e2))
  attr(out, "channel_groups") <- cg
  attr(out, "model_topology") <- topo
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

#' Set the Random Effect Correlations of the Model
#'
#' @param e1 A `dynamiteformula` object.
#' @param e2 A `random_spec` object.
#' @noRd
set_random_spec <- function(e1, e2) {
  stopifnot_(
    is.null(attr(e1, "random_spec")) || attr(e2, "random_spec"),
    "Multiple definitions for random effect specifications."
  )
  attr(e1, "random_spec") <- e2
  e1
}

#' Set the Latent Factors of the Model
#'
#' @param e1 A `dynamiteformula` object.
#' @param e2 A `latent_factor` object.
#' @noRd
set_lfactor <- function(e1, e2) {
  stopifnot_(
    is.null(attr(e1, "lfactor")) || attr(e2, "lfactor"),
    "Multiple definitions for latent factors."
  )
  attr(e1, "lfactor") <- e2
  e1
}

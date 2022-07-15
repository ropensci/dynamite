#' Estimate a Bayesian Dynamic Multivariate Panel Model
#'
#' Any univariate unbounded continuous distributions supported by Stan can be
#' used as a prior for model parameters (the distribution is automatically
#' truncated to positive side for constrained parameters). In addition, any
#' univariate distribution bounded to positive real line can be used as a prior
#' for parameters constrained to be positive. See Stan function reference at
#' \url{https://mc-stan.org/users/documentation/} for details. For custom
#' priors, you should first get the default priors with `get_priors` function,
#' and then modify the `priors` column of the obtained data frame before
#' supplying it to the `dynamite`.
#'
#' The default priors for regression coefficients are based on the standard
#' deviation of the covariates at the first non-fixed time point. In case this
#' is 0 or NA, it is transformed to (arbitrary) 0.5. The final prior is then
#' normal distribution with zero mean and two times this standard deviation.
#'
#' The prior for the correlation structure of the random intercepts is defined
#' via the Cholesky decomposition of the correlation matrix, as
#' `lkj_corr_cholesky(1)`. See
#' \url{https://mc-stan.org/docs/functions-reference/cholesky-lkj-correlation-distribution.html}
#' for details.
#'
#' See more details in the package vignette on how to define a dynamite model.
#'
#' The best-case scalability of the dynamite in terms of data size should be
#' approximately linear in terms of number of time points and and number of
#' groups, but as wall-clock time of the MCMC algorithms provided by Stan can
#' depend on the discrepancy of the data and the model (and the subsequent
#' shape of the posterior), this can vary greatly.
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula. See 'Details'.
#'   Input `data` is converted to channel specific matrix representations
#'   via [stats::model.matrix.lm].
#' @param data
#'   \[`data.frame`, `tibble::tibble`, or `data.table::data.table`]\cr
#'   The data frame, tibble or a data.table containing the variables in the
#'   model. Supported column types are `integer`, `logical`, `double`,
#'   `factor`. `character` columns will be converted to factors.
#'   Unused factor levels will be dropped. The `data` can contain missing
#'   values which will simply be ignored in the estimation in a case-wise
#'   fashion (per time-point and per channel).
#' @param group \[`character(1)`]\cr A column name of `data` that denotes the
#'   unique groups, or `NULL` corresponding to a scenario without grouping.
#' @param time \[`character(1)`]\cr A column name of `data` that denotes the
#'   time index of observations. If this variable is a factor, the integer
#'   representation of its levels are used internally for defining the time
#'   indexing.
#' @param priors \[`data.frame`]\cr An optional data frame with prior
#'   definitions. See 'Details'.
#' @param verbose \[`logical(1)`]\cr All warnings and messages are suppressed
#'   if set to `FALSE`. Defaults to `TRUE`.
#' @param debug \[`list()`]\cr A named list of form `name = TRUE` indicating
#'   additional objects in the environment of the `dynamite` function which are
#'   added to the return object. Additionally, values `no_compile = TRUE` and
#'   `no_sampling = TRUE` can be used to skip the compilation of the Stan code
#'   and sampling steps respectively. This can be useful for debugging when
#'   combined with `model_code = TRUE`.
#' @param ... Additional arguments to [rstan::sampling()].
#' @export
#' @examples
#' \dontrun{
#' fit <- dynamite(obs(y ~ -1 + varying(~x), family = "gaussian") +
#'   lags(type = "varying") + splines(df = 20), gaussian_example, "id", "time",
#'   chains = 1, refresh = 0)
#'
#' library(dplyr)
#' library(ggplot2)
#' cf <- coef(fit) %>%
#'   group_by(time, variable) %>%
#'   summarise(
#'     mean = mean(value),
#'     lwr = quantile(value, 0.025),
#'     upr = quantile(value, 0.975))
#'
#' cf %>%
#'   ggplot(aes(time, mean)) + theme_bw() +
#'     geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.7) +
#'     geom_line() +
#'     facet_wrap(~ variable, scales = "free_y")
#' }
#'
#' @srrstats {G2.9} Potential loss of information is reported by `dynamite`.
#' @srrstats {RE1.1} Documented in `dformula` parameter.
#' @srrstats {RE1.4} Documented in `data` parameter.
#' @srrstats {RE2.0} Transformations are documented (in parameter `data`) and
#'   warnings are issued when appropriate.
#' @srrstats {RE2.1} Non-finite values are not allowed, NAs are appropriately
#'   considered.
#' @srrstats {RE3.0} Convergence diagnostics are delegated to Stan
#' @srrstats {RE3.1} Stan messages can be suppressed
#' @srrstats {RE4.0} `dynamite` returns a `dynamitefit` object.
#' @srrstats {RE4.1} `dynamitefit` object can be generated without sampling.
#' @srrstats {RE4.4} The model is specified via formula objects.
#' @srrstats {RE4.8} `dynamitefit` object contains the response variables
#' @srrstats {RE4.13} `dynamitefit` object contains the predictor variables
#' @srrstats {BS1.1} Data input is documented in the `data` parameter.,
#' @srrstats {BS1.3, BS1.3b} Computational parameters are delegated to Stan.
#' @srrstats {BS2.6} Checks for computational parameters are performed by Stan.
#' @srrstats {BS2.7} Starting values can be controlled via `...`.
#' @srrstats {BS2.9} Chains have different starting values by default.
#' @srrstats {BS2.12} Verbosity of output can be controlled via `...`.
#' @srrstats {BS2.13} Progress indicators can be suppressed without suppressing
#'   other output, via `...`.
#' @srrstats {BS2.14} Warnings can be suppressed by setting `verbose = FALSE`-
#' @srrstats {BS3.0} Documented in `data` parameter.
#' @srrstats {BS4.0} Stan is referenced.
#' @srrstats {BS5.0, BS5.1, BS5.2, BS5.3, BS5.5}
#'   Available from the resulting `dynamitefit` object
#' @srrstats {RE5.0} The scalability of the algorithms are studied to some
#' extent in the tests and note also here. As the computational algorithms are
#' based on Stan, the  scalability of the package depends directly on the
#' scalability of Stan.
dynamite <- function(dformula, data, group = NULL, time,
                     priors = NULL, verbose = TRUE, debug = NULL, ...) {
  stopifnot_(
    is.dynamiteformula(dformula),
    "Argument {.arg dformula} must be a {.cls dynamiteformula} object."
  )
  stopifnot_(
    is.data.frame(data),
    "Argument {.arg data} must be a {.cls data.frame} object."
  )
  stopifnot_(
    checkmate::test_string(
      x = group,
      null.ok = TRUE
    ),
    "Argument {.arg group} must be a single character string or {.code NULL}."
  )
  if (!is.null(group)) {
    stopifnot_(
      !is.null(data[[group]]),
      "Can't find grouping variable {.var {group}} in {.arg data}."
    )
  }
  stopifnot_(
    !missing(time),
    "Argument {.var time} is missing."
  )
  stopifnot_(
    checkmate::test_string(x = time),
    "Argument {.arg time} must be a single character string."
  )
  stopifnot_(
    !is.null(data[[time]]),
    "Can't find time index variable {.var {time}} in {.arg data}."
  )
  stopifnot_(
    checkmate::test_flag(x = verbose),
    "Argument {.arg verbose} must be a single {.cls logical} value."
  )
  data <- parse_data(data, dformula, group, time, verbose)
  dformulas <- parse_lags(data, dformula, group, time)
  evaluate_deterministic(data, dformulas, group, time)
  stan <- prepare_stan_data(
    data,
    dformulas$stoch,
    group,
    time,
    priors,
    fixed = attr(dformulas$all, "max_lag")
  )
  model_code <- create_blocks(
    dformula = dformulas$stoch,
    indent = 2L,
    vars = stan$model_vars
  )
  model <- onlyif(
    is.null(debug) || !isTRUE(debug$no_compile),
    {
      onlyif(verbose, rlang::inform("Compiling Stan model."))
      rstan::stan_model(model_code = model_code)
    }
  )
  # don't save redundant parameters by default
  # could also remove omega_raw
  dots <- list(...)
  if (is.null(dots$pars) && stan$sampling_vars$M > 0) {
    dots$pars <- c("nu_raw", "nu", "L")
    dots$include <- FALSE
  }
  stanfit <- onlyif(
    !isTRUE(debug$no_compile) && !isTRUE(debug$no_sampling),
    do.call(rstan::sampling,
      c(list(object = model, data = stan$sampling_vars), dots))
  )
  out <- structure(
    list(
      stanfit = stanfit,
      dformulas = dformulas,
      data = data,
      stan = stan,
      group_var = group,
      time_var = time,
      priors = dplyr::bind_rows(stan$priors),
      verbose
    ),
    class = "dynamitefit"
  )
  # Adds any object in the environment of this function to the return object
  # if its name is included in the debug argument in the form `name = TRUE`
  for (opt in names(debug)) {
    if (debug[[opt]]) {
      got <- try(get(x = opt), silent = TRUE)
      out[[opt]] <- onlyif(!"try-error" %in% class(got), got)
    }
  }
  out
}

#' Access The Model Formula of a Dynamite Model
#'
#' Returns the model definition as a quoted expression.
#'
#' @param x \[`dynamitefit`\]\cr The model fit object.
#' @param ... Not used.
#' @export
#' @examples
#' formula(gaussian_example_fit)
formula.dynamitefit <- function(x, ...) {
  formula_str <- vapply(
    get_originals(x$dformulas$all),
    deparse1,
    character(1L)
  )
  ch_stoch <- which_stochastic(x$dformulas$all)
  ch_det <- which_deterministic(x$dformulas$all)
  family_str <- vapply(
    get_families(x$dformulas$all),
    function(y) y$name,
    character(1L)
  )
  lag_defs <- attr(x$dformulas$all, "lags")
  spline_defs <- attr(x$dformulas$stoch, "splines")
  obs_str <- ""
  aux_str <- ""
  lags_str <- ""
  spline_str <- ""
  n_stoch <- length(ch_stoch)
  if (n_stoch > 0L) {
    obs_str <- paste0(
      glue::glue(
        "obs({formula_str[ch_stoch]}, family = {family_str[ch_stoch]}())"
      ),
      collapse = " +\n"
    )
  }
  n_det <- length(ch_det)
  if (n_det > 0L) {
    aux_str <- paste0(
      glue::glue("aux({formula_str[ch_det]})"),
      collapse = " +\n"
    )
  }
  if (!is.null(lag_defs)) {
    lags_str <- glue::glue("lags(k = {lag_defs$k}, type = {lag_defs$type})")
  }
  if (!is.null(spline_defs)) {
    spline_str <- paste0(
      "splines(",
      "shrinkage = ", spline_defs$shrinkage, ", ",
      "override = FALSE, ",
      "df = ", spline_defs$bs_opts$df, ", ",
      "degree = ", spline_defs$bs_opts$degree, ", ",
      "lb_tau = ", spline_defs$lb_tau, ", ",
      "noncentered = ",  spline_defs$noncentered, ")"
    )
  }
  str2lang(
    paste0(
      "{\n",
      paste(obs_str, aux_str, lags_str, spline_str, sep = " +\n"),
      "\n}"
    )
  )
}

#' Is The Argument a `dynamitefit` Object
#'
#' @param x An \R object.
#' @noRd
is.dynamitefit <- function(x) {
  inherits(x, "dynamitefit")
}

#' Parse Data for Model Fitting
#'
#' @inheritParams dynamite
#' @param group_var \[`character(1)`] Grouping variable name.
#' @param time_var \[`character(1)`] Time index variable name.
#' @srrstats {G2.4d, G2.5} Factors and ordered factors are considered.
#' @srrstats {G2.6} Columns are preprocessed.
#' @srrstats {G2.11} Data checks rely on type, not class, except for factors.
#' @srrstats {G2.12} List columns are not supported.
#' @srrstats {G2.16} Non-finite values are not supported.
#' @noRd
parse_data <- function(data, dformula, group_var, time_var, verbose) {
  data <- droplevels(data)
  data_names <- names(data)
  stopifnot_(!is.character(data[[time_var]]),
    "Time indexing variable {.arg {time_var}} is of type {.cls character}, but
    it should be of type {.cls numeric} or {.cls factor}.")
  data <- data |> dplyr::mutate(dplyr::across(where(is.character), as.factor))
  valid_types <- c("integer", "logical", "double")
  col_types <- vapply(data, typeof, character(1L))
  factor_cols <- vapply(data, is.factor, logical(1L))
  valid_cols <- (col_types %in% valid_types) | factor_cols
  if (any(!valid_cols)) {
    invalid_cols <- data_names[!valid_cols]
    invalid_types <- col_types[!valid_cols]
    stop_(c(
      "Column{?s} {.var {invalid_cols}} of {.arg data} {?is/are} invalid:",
      `x` = "Column type{?s} {.cls {invalid_types}} {?is/are} not supported."
    ))
  }
  coerce_cols <- valid_cols & !factor_cols
  if (any(coerce_cols)) {
    for (i in which(coerce_cols)) {
      data[,i] <- do.call(
        paste0("as.", typeof(data[[i]])),
        args = list(data[[i]])
      )
    }
  }
  resp <- get_responses(dformula)
  ordered_factor_resp <- vapply(
    seq_along(resp),
    function(i) {
      is_categorical(dformula[[i]]$family) && is.ordered(data[, resp[i]])
    },
    logical(1L)
  )
  if (any(ordered_factor_resp)) {
    rof <- resp[ordered_factor_resp]
    if (verbose) {
      warning_(c(
        "Response variable{?s} {.var {rof}} {?is/are} of class
         {.cls ordered factor} whose channel{?s} {?is/are} categorical:",
        `i` = "{.var {rof}} will be converted to {?an/} unordered factor{?s}."
      ))
    }
    for (i in seq_along(rof)) {
      class(data[, rof[i]]) <- "factor"
    }
  }
  finite_cols <- vapply(
    data,
    function(x) all(is.finite(x) | is.na(x)),
    logical(1L)
  )
  stopifnot_(
    all(finite_cols),
    "Non-finite values in variable{?s} {.var {data_names[!finite_cols]}} of
     {.arg data}."
  )
  data <- fill_time(data, group_var, time_var)
  data <- data.table::as.data.table(data)
  data.table::setkeyv(data, c(group_var, time_var))
  data
}

#' Parse Lag and Lags Definitions of a `dynamiteformula` Object
#'
#' Also processes the checks on random component and adds it to dformulas.
#' @inheritParams parse_data
#' @noRd
parse_lags <- function(data, dformula, group_var, time_var) {
  channels_det <- which_deterministic(dformula)
  channels_stoch <- which_stochastic(dformula)
  resp_all <- get_responses(dformula)
  resp_stoch <- resp_all[channels_stoch]
  n_rows <- data[, .N]
  n_channels <- length(resp_all)
  max_lag <- 0L
  for (i in seq_len(n_channels)) {
    fix_rhs <- complete_lags(formula_rhs(dformula[[i]]$formula))
    dformula[[i]]$formula <- as.formula(
      paste0(resp_all[i], "~", deparse1(fix_rhs))
    )
  }
  data_names <- names(data)
  non_lags <- extract_nonlags(
    unique(unlist(get_terms(dformula[channels_stoch])))
  )
  valid_resp <- c(resp_stoch, data_names)
  mis_vars <- which(!non_lags %in% valid_resp)
  stopifnot_(
    identical(length(mis_vars), 0L),
    "Can't find variable{?s} {.var {non_lags[mis_vars]}} in {.arg data}."
  )
  predictors <- get_predictors(dformula)
  lag_map <- extract_lags(predictors)
  gl <- parse_global_lags(dformula, lag_map, resp_stoch, channels_stoch)
  dformula <- gl$dformula
  lag_map <- gl$lag_map
  max_lag <- gl$max_lag
  mis_lags <- which(!lag_map$var %in% c(resp_all, data_names))
  stopifnot_(
    identical(length(mis_lags), 0L),
    c(
      "Unable to construct lagged values of
       {.var {cs(lag_map$var[mis_lags])}}:",
      `x` = "Can't find such variables in {.var data}."
    )
  )
  stoch_k <- lag_map |>
    dplyr::filter(.data$var %in% resp_stoch) |>
    dplyr::pull(.data$k)
  if (length(stoch_k) > 0) {
    max_lag <- max(max_lag, lag_map$k[stoch_k])
  }
  sl <- parse_singleton_lags(dformula, lag_map, max_lag, valid_resp)
  dformula <- sl$dformula
  dformula_det <- dformula[channels_det]
  dformula_lag_pred <- sl$channels[sl$pred]
  dformula_lag_stoch <- c(
    gl$channels[gl$stoch],
    sl$channels[sl$stoch & !sl$pred]
  )
  dformula_lag_det <- c(
    gl$channels[!gl$stoch],
    sl$channels[!sl$stoch & !sl$pred]
  )
  attr(dformula_lag_pred, "rank_order") <- order(sl$rank[sl$pred])
  attr(dformula_lag_det, "rank_order") <-
    order(c(gl$rank[!gl$stoch], sl$rank[!sl$stoch & !sl$pred]))
  attr(dformula_lag_pred, "original_response") <- sl$resp[sl$pred]
  attr(dformula, "max_lag") <- max_lag
  random_defs <- attr(dformula, "random")
  if (!is.null(random_defs)) {
    families <- unlist(get_families(dformula[channels_stoch]))
    valid_channels <- resp_stoch[!(families %in% "categorical")]
    # default, use all channels except categorical
    if (is.null(random_defs$channels)) {
      random_defs$channels <- valid_channels
      stopifnot_(length(valid_channels) > 0,
        c("No valid channels for random intercept component:",
          `x` = "Random intercepts are not supported for the categorical
          family."
        )
      )
    } else {
      nu_channels <- random_defs$channels %in% resp_stoch
      stopifnot_(all(nu_channels),
        c("Argument {.arg channel} of {.fun random} contains variables
       {.var {cs(resp_stoch[nu_channels])}}:",
          `x` = "No such response variables in the model."
        )
      )
      nu_channels <- random_defs$channels %in% valid_channels
      stopifnot_(all(nu_channels),
        c("Random intercepts are not supported for the categorical family:",
          `x` = "Found random intercept declaration for categorical variable{?s}
          {.var {cs(valid_channels[nu_channels])}}."
        )
      )
    }
    if(length(random_defs$channels) < 2) {
      random_defs$correlated <- FALSE
    }
  }
  list(
    all = dformula,
    det = dformula_det,
    stoch = structure(
      dformula[channels_stoch],
      splines = attr(dformula, "splines"),
      random = random_defs
    ),
    lag_pred = dformula_lag_pred,
    lag_det = dformula_lag_det,
    lag_stoch = dformula_lag_stoch
  )
}

#' Parse a lags definition in a dynamiteformula
#'
#' @inheritParams parse_data
#' @param lag_map Output of `extract_lags`.
#' @param resp_stoch A `character` of stochastic response variable names.
#' @param channels_stoch
#'   A `logical` vector indicating which channels are stochastic.
#' @noRd
parse_global_lags <- function(dformula, lag_map, resp_stoch, channels_stoch) {
  lags_def <- attr(dformula, "lags")
  idx <- 0L
  k <- lags_def$k
  type <- lags_def$type
  resp_all <- get_responses(dformula)
  max_lag <- ifelse_(is.null(lags_def), 0L, max(k))
  n_lag <- max_lag * length(resp_stoch)
  channels <- vector(mode = "list", length = n_lag)
  stoch <- logical(n_lag)
  rank <- integer(n_lag)
  lhs <- character(n_lag)
  increment <- logical(n_lag)
  for (i in seq_len(max_lag)) {
    for (j in seq_along(channels_stoch)) {
      y <- resp_stoch[j]
      idx <- idx + 1L
      lhs[idx] <- paste0(y, "_lag", i)
      if (identical(i, 1L)) {
        rhs <- y
        stoch[idx] <- TRUE
      } else {
        rhs <- paste0(y, "_lag", i - 1L)
      }
      rank[idx] <- i
      increment[idx] <- i %in% k
      channels[[idx]] <- dynamitechannel(
        formula = as.formula(paste0(lhs[idx], " ~ ", rhs)),
        family = deterministic_(),
        response = lhs[idx]
      )
    }
  }
  if (!is.null(lags_def)) {
    for (j in channels_stoch) {
      dformula[[j]] <- dynamiteformula_(
        formula = increment_formula(
          formula = dformula[[j]]$formula,
          x = lhs[increment],
          type = type,
          varying_idx = dformula[[j]]$varying,
          varying_icpt = dformula[[j]]$has_varying_intercept,
          fixed_icpt = dformula[[j]]$has_fixed_intercept
        ),
        original = dformula[[j]]$original,
        family = dformula[[j]]$family
      )
    }
    lag_map <- lag_map |>
      dplyr::filter(!(.data$var %in% resp_all & .data$k %in% k))
  }
  list(
    dformula = dformula,
    channels = channels,
    lag_map = lag_map,
    max_lag = max_lag,
    rank = rank,
    stoch = stoch
  )
}

#' Parse manual lag terms in a dynamiteformula
#'
#' @inheritParams parse_global_lags
#' @param max_lag Largest shift value of the model in any lag.
#' @param valid_resp A `character` vector of valid LHS variables that can
#'   appear in the model formulas.
#' @noRd
parse_singleton_lags <- function(dformula, lag_map, max_lag, valid_resp) {
  n_lag <- nrow(lag_map)
  n_resp <- length(dformula)
  resp_all <- get_responses(dformula)
  channels <- vector(mode = "list", length = n_lag)
  resp_lag <- character(n_lag)
  pred <- logical(n_lag)
  stoch <- logical(n_lag)
  rank <- integer(n_lag)
  lag_var <- unique(lag_map$var)
  idx <- 0L
  for (resp in lag_var) {
    y <- prepare_lagged_response(
      dformula,
      lag_map,
      max_lag,
      resp,
      resp_all,
      valid_resp
    )
    if (y$deterministic) {
      dformula[[y$idx]]$specials$past <- NULL
    }
    for (i in seq_along(y$lag_idx)) {
      idx <- idx + 1L
      rhs <- ifelse_(identical(i, 1L), resp, paste0(resp, "_lag", i - 1L))
      lhs <- paste0(resp, "_lag", i)
      rank[idx] <- i
      stoch[idx] <- !y$deterministic
      pred[idx] <- !y$is_resp
      resp_lag[idx] <- resp
      spec <- NULL
      if (y$is_resp && !is.null(y$past) && i > y$past_offset) {
        y$past_idx <- y$past_idx + 1L
        spec <- list(
          past = y$past[y$past_idx],
          past_offset = max_lag,
          resp_type = y$type
        )
      }
      channels[[idx]] <- dynamitechannel(
        formula = as.formula(paste0(lhs, " ~ ", rhs)),
        family = deterministic_(),
        response = lhs,
        specials = spec
      )
      dformula <- parse_present_lags(dformula, lag_map, y, i, lhs)
    }
  }
  list(
    dformula = dformula,
    channels = channels,
    pred = pred,
    rank = rank,
    resp = resp_lag,
    stoch = stoch
  )
}

#' Parse lags that actually appear in a dynamiteformula
#'
#' @inheritParams parse_singleton_lags
#' @param y Output of `prepare_lagged_response`.
#' @param i Row index of lag_map
#' @param lhs Name of the new lagged response.
#' @noRd
parse_present_lags <- function(dformula, lag_map, y, i, lhs) {
  k <- y$lag_idx[i]
  if (lag_map$present[k]) {
    for (j in seq_along(dformula)) {
      dformula[[j]]$formula <- gsub_formula(
        pattern = lag_map$src[k],
        replacement = lhs,
        formula = dformula[[j]]$formula,
        fixed = TRUE
      )
    }
  }
  dformula
}

#' Prepare a response variable for a new lag channel
#'
#' @inheritParams parse_singleton_lags
#' @param resp Name of the response variable being lagged
#' @param resp_all All responses in the model as a `character` vector.
#' @noRd
prepare_lagged_response <- function(dformula, lag_map, max_lag, resp,
                                    resp_all, valid_resp) {
  n_resp <- length(dformula)
  y <- list()
  y_obs <- NULL
  y_self <- NULL
  y$resp <- resp
  y$lag_idx <- which(lag_map$var == resp)
  y$src <- lag_map$src[y$lag_idx]
  y$idx <- which(resp_all == resp)
  y$is_resp <- length(y$idx) > 0L
  y$deterministic <- FALSE
  if (y$is_resp) {
    y$past <- NULL
    y$past_idx <- NULL
    y$past_offset <- NULL
    y$deterministic <- is_deterministic(dformula[[y$idx]]$family)
    y$type <- dformula[[y$idx]]$specials$resp_type
    if (y$deterministic) {
      y_form <- deparse1(dformula[[y$idx]]$formula)
      y$past <- dformula[[y$idx]]$specials$past
      y_past_len <- length(y$past)
      y_self <- max(extract_self_lags(y_form, resp))
      y_max <- max(lag_map$k[y$lag_idx])
      y_obs <- extract_lags(y_form) |>
        dplyr::filter(.data$var %in% valid_resp) |>
        dplyr::pull(.data$k)
      y_obs <- ifelse_(length(y_obs) > 0L, max(y_obs), 0L)
      y_req_past <- min(y_max - max_lag + max(y_self, y_obs), y_max)
      y$past_offset <- y_max - y_req_past
      stopifnot_(
        y_past_len >= y_req_past,
        c(
          "Deterministic channel {.var {resp}} requires {y_req_past} initial
           value{?s}:",
          `x` = "You've supplied {cli::no(y_past_len)}
                 initial value{?s//s}."
        )
      )
      y$past_idx <- 0L
    }
  }
  y
}

#' Evaluate Definitions of Deterministic Channels
#'
#' @inheritParams parse_data
#' @param dformulas The return object of [dynamite::parse_lags()].
#' @noRd
evaluate_deterministic <- function(data, dformulas, group_var, time_var) {
  fixed <- as.integer(attr(dformulas$all, "max_lag"))
  n_time <- length(unique(data[[time_var]]))
  n_id <- ifelse_(is.null(group_var), 1L, length(unique(data[[group_var]])))
  dd <- dformulas$det
  dlp <- dformulas$lag_pred
  dld <- dformulas$lag_det
  dls <- dformulas$lag_stoch
  cl <- get_quoted(dd)
  initialize_deterministic(data, dd, dlp, dld, dls)
  idx <- seq.int(1L, n_time * n_id, by = n_time) - 1L
  assign_initial_values(data, dd, dlp, dld, dls, idx, fixed, group_var)
  if (n_time > fixed + 1L) {
    ro_stoch <- seq_len(length(dls))
    ro_det <- ifelse_(
      is.null(attr(dld, "rank_order")),
      integer(0L),
      attr(dld, "rank_order")
    )
    lhs_det <- get_responses(dld)
    rhs_det <- get_predictors(dld)
    lhs_stoch <- get_responses(dls)
    rhs_stoch <- get_predictors(dls)
    idx <- idx + fixed + 1L
    for (i in seq.int(fixed + 2L, n_time)) {
      idx <- idx + 1L
      assign_lags(data, ro_det, idx, lhs_det, rhs_det)
      assign_lags(data, ro_stoch, idx, lhs_stoch, rhs_stoch)
      assign_deterministic(data, cl, idx)
    }
  }
}

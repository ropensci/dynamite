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
#' \url{https://mc-stan.org/docs/functions-reference/
#' cholesky-lkj-correlation-distribution.html}
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
#' @param data
#'   \[`data.frame`, `tibble::tibble`, or `data.table::data.table`]\cr
#'   The data frame, tibble or a data.table containing the variables in the
#'   model. Supported column types are `integer`, `logical`, `double`,
#'   `factor`. `character` columns will be converted to factors.
#'   Unused factor levels will be dropped. The `data` can contain missing
#'   values which will simply be ignored in the estimation in a case-wise
#'   fashion (per time-point and per channel). Input `data` is converted to
#'   channel specific matrix representations via [stats::model.matrix.lm].
#' @param group \[`character(1)`]\cr A column name of `data` that denotes the
#'   unique groups, or `NULL` corresponding to a scenario without any groups.
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
#'   combined with `model_code = TRUE`, which adds the Stan model code to the
#'   return object.
#' @param ... Additional arguments to [rstan::sampling()].
#' @export
#' @examples
#' \dontrun{
#' fit <- dynamite(
#'   dformula = obs(y ~ -1 + varying(~x), family = "gaussian") +
#'     lags(type = "varying") +
#'     splines(df = 20), gaussian_example, "id", "time",
#'   chains = 1,
#'   refresh = 0
#' )
#'
#' library(dplyr)
#' library(ggplot2)
#' cf <- coef(fit) %>%
#'   group_by(time, variable) %>%
#'   summarise(
#'     mean = mean(value),
#'     lwr = quantile(value, 0.025),
#'     upr = quantile(value, 0.975)
#'   )
#'
#' cf %>%
#'   ggplot(aes(time, mean)) +
#'   theme_bw() +
#'   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.7) +
#'   geom_line() +
#'   facet_wrap(~variable, scales = "free_y")
#' }
#'
#' @srrstats {G2.9} Potential loss of information is reported by `dynamite`.
#' @srrstats {RE1.1} Documented in `dformula` parameter.
#' @srrstats {RE1.4} Documented in `data` parameter.
#' @srrstats {RE2.0} Transformations are documented (in parameter `data`) and
#'   warnings are issued when appropriate.
#' @srrstats {RE2.1} Non-finite values are not allowed, NAs are appropriately
#'   considered.
#' @srrstats {RE3.0} Convergence diagnostics are delegated to Stan.
#' @srrstats {RE3.1} Stan messages can be suppressed.
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
#'   Available from the resulting `dynamitefit` object.
#' @srrstats {RE5.0} The scalability of the algorithms is studied to some
#'   extent in the tests and noted here. As the computational algorithms are
#'   based on Stan, the  scalability of the package depends directly on the
#'   scalability of Stan.
dynamite <- function(dformula, data, group = NULL, time,
                     priors = NULL, verbose = TRUE, debug = NULL, ...) {
  stopifnot_(
    !missing(dformula),
    "Argument {.arg dformula} is missing."
  )
  stopifnot_(
    !missing(data),
    "Argument {.arg data} is missing."
  )
  stopifnot_(
    !missing(time),
    "Argument {.var time} is missing."
  )
  stopifnot_(
    is.dynamiteformula(dformula),
    "Argument {.arg dformula} must be a {.cls dynamiteformula} object."
  )
  stopifnot_(
    length(which_stochastic(dformula)) > 0L,
    "Argument {.arg dformula} must contain at least one stochastic channel."
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
  stopifnot_(
    is.null(group) || !is.null(data[[group]]),
    "Can't find grouping variable {.var {group}} in {.arg data}."
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
  data_name <- deparse1(substitute(data))
  data <- parse_data(dformula, data, group, time, verbose)
  dformula <- parse_past(dformula, data, group, time)
  dformulas <- parse_lags(dformula, data, group, time, verbose)
  evaluate_deterministic(dformulas, data, group, time)
  stan <- prepare_stan_input(
    dformulas$stoch,
    data,
    group,
    time,
    priors,
    fixed = attr(dformulas$all, "max_lag"),
    verbose
  )
  model_code <- create_blocks(
    dformula = dformulas$stoch,
    indent = 2L,
    vars = stan$model_vars
  )
  model <- onlyif(
    is.null(debug) || !isTRUE(debug$no_compile),
    {
      onlyif(verbose, message_("Compiling Stan model."))
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
    do.call(
      rstan::sampling,
      c(list(object = model, data = stan$sampling_vars), dots)
    )
  )
  out <- structure(
    list(
      stanfit = stanfit,
      dformulas = dformulas,
      data = data,
      data_name = data_name,
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
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
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
  obs_str <- onlyif(
    length(ch_stoch) > 0L,
    obs_str <- paste0(
      glue::glue(
        "obs({formula_str[ch_stoch]}, family = {family_str[ch_stoch]}())"
      ),
      collapse = " +\n"
    )
  )
  aux_str <- onlyif(
    length(ch_det) > 0L,
    aux_str <- paste0(
      glue::glue("aux({formula_str[ch_det]})"),
      collapse = " +\n"
    )
  )
  lags_str <- onlyif(
    !is.null(lag_defs),
    glue::glue("lags(k = {lag_defs$k}, type = {lag_defs$type})")
  )
  spline_str <- onlyif(
    !is.null(spline_defs),
    paste0(
      "splines(",
      "shrinkage = ", spline_defs$shrinkage, ", ",
      "override = FALSE, ",
      "df = ", spline_defs$bs_opts$df, ", ",
      "degree = ", spline_defs$bs_opts$degree, ", ",
      "lb_tau = ", spline_defs$lb_tau, ", ",
      "noncentered = ", spline_defs$noncentered, ")"
    )
  )
  str2lang(
    paste0(
      "{\n",
      paste_rows(obs_str, aux_str, lags_str, spline_str),
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
parse_data <- function(dformula, data, group_var, time_var, verbose) {
  data <- droplevels(data)
  data_names <- names(data)
  stopifnot_(
    !is.character(data[[time_var]]),
    c(
      "Time index variable {.arg {time_var}} is of type {.cls character}:",
      `i` = "It must be of type {.cls numeric} or {.cls factor}."
    )
  )
  if (is.factor(data[[time_var]])) {
    if (verbose) {
      warning_(c(
        "Time index variable {.arg {time_var}} is a {.cls factor}:",
        `i` = "Converting the variable to {.cls integer} based on its levels."
      ))
    }
    data[[time_var]] <- as.integer(data[[time_var]])
  }
  data <- data |> dplyr::mutate(dplyr::across(where(is.character), as.factor))
  valid_types <- c("integer", "logical", "double")
  col_types <- vapply(data, typeof, character(1L))
  factor_cols <- vapply(data, is.factor, logical(1L))
  valid_cols <- (col_types %in% valid_types) | factor_cols
  stopifnot_(
    all(valid_cols),
    c(
      "Column{?s} {.var {data_names[!valid_cols]}} of {.arg data}
      {?is/are} invalid:",
      `x` = "Column type{?s} {.cls {col_types[!valid_cols]}}
             {?is/are} not supported."
    )
  )
  for (i in which(valid_cols & !factor_cols)) {
    data[, i] <- do.call(
      paste0("as.", typeof(data[[i]])),
      args = list(data[[i]])
    )
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
    "Non-finite values were found in variable{?s}
    {.var {data_names[!finite_cols]}} of {.arg data}."
  )
  data <- fill_time(data, group_var, time_var)
  data <- data.table::as.data.table(data)
  drop_unused(dformula, data, group_var, time_var)
  data.table::setkeyv(data, c(group_var, time_var))
  data
}

#' Evaluate Past Value Definitions of Each Deterministic Channel
#'
#' @inheritParams parse_data
#' @noRd
parse_past <- function(dformula, data, group_var, time_var) {
  n_time <- length(unique(data[[time_var]]))
  n_group <- ifelse_(is.null(group_var), 1L, length(unique(data[[group_var]])))
  past <- has_past(dformula)
  for (i in seq_along(dformula)) {
    if (past[i]) {
      y <- dformula[[i]]$response
      past_eval <- try(data[, eval(dformula[[i]]$specials$past)], silent = TRUE)
      stopifnot_(
        !"try-error" %in% class(past_eval),
        c(
          "Unable to evaluate past definition of
          deterministic channel {.var {y}}:",
          `x` = attr(past_eval, "condition")$message
        )
      )
      past_type <- dformula[[i]]$specials$past_type
      past_len <- length(past_eval)
      stopifnot_(
        identical(past_type, "init") || identical(past_len, nrow(data)),
        c(
          "Incompatible past definition of deterministic channel {.var {y}}:",
          `x` = "The definition evaluates to length {past_len} but the data has
                {nrow(data)} rows."
        )
      )
      dformula[[i]]$specials$past <- past_eval
    } else {
      dformula[[i]]$specials$past_type <- "init"
    }
  }
  dformula
}

#' Parse Lag and Lags Definitions of a `dynamiteformula` Object
#'
#' Also processes the checks on random component and adds it to dformulas.
#'
#' @inheritParams parse_data
#' @noRd
parse_lags <- function(dformula, data, group_var, time_var, verbose) {
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
  max_lag <- ifelse_(
    length(stoch_k) > 0L,
    max(max_lag, stoch_k),
    max_lag
  )
  sl <- parse_singleton_lags(
    dformula,
    data,
    group_var,
    lag_map,
    valid_resp,
    verbose
  )
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
    if (is.null(random_defs$responses)) {
      random_defs$responses <- valid_channels
      stopifnot_(
        length(valid_channels) > 0L,
        c(
          "No valid responses for random intercept component:",
          `x` = "Random intercepts are not supported for the categorical
                family."
        )
      )
    } else {
      nu_channels <- random_defs$responses %in% resp_stoch
      stopifnot_(
        all(nu_channels),
        c(
          "Argument {.arg responses} of {.fun random} contains variables
          {.var {cs(resp_stoch[nu_channels])}}:",
          `x` = "No such response variables in the model."
        )
      )
      nu_channels <- random_defs$responses %in% valid_channels
      stopifnot_(
        all(nu_channels),
        c(
          "Random intercepts are not supported for the categorical family:",
          `x` = "Found random intercept declaration for categorical variable{?s}
                {.var {cs(valid_channels[nu_channels])}}."
        )
      )
    }
    if (length(random_defs$responses) < 2L) {
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

#' Parse a `lags` Definition in a `dynamiteformula` Object
#'
#' @inheritParams parse_lags
#' @param lag_map \[`data.frame`]\cr
#'   Output of `extract_lags`.
#' @param resp_stoch \[`character()`]\cr
#'   A vector of stochastic response variable names.
#' @param channels_stoch \[`logical()`]\cr
#'   A vector indicating which channels are stochastic.
#' @noRd
parse_global_lags <- function(dformula, lag_map, resp_stoch, channels_stoch) {
  lags_def <- attr(dformula, "lags")
  idx <- 0L
  k <- lags_def$k
  type <- lags_def$type
  resp_all <- get_responses(dformula)
  max_lag <- ifelse_(is.null(lags_def), 0L, max(k))
  n_stoch <- length(resp_stoch)
  n_lag <- max_lag * n_stoch
  channels <- vector(mode = "list", length = n_lag)
  dterms <- get_terms(dformula[channels_stoch])
  stoch <- logical(n_lag)
  rank <- integer(n_lag)
  lhs <- character(n_lag)
  include <- logical(n_lag)
  increment <- replicate(n_stoch, logical(n_lag), simplify = FALSE)
  for (i in seq_len(max_lag)) {
    for (j in seq_len(n_stoch)) {
      y <- resp_stoch[j]
      idx <- idx + 1L
      stoch[idx] <- identical(i, 1L)
      lhs[idx] <- paste0(y, "_lag", i)
      rank[idx] <- i
      rhs <- ifelse_(stoch[idx], y, paste0(y, "_lag", i - 1L))
      new_term <- logical(n_stoch)
      lag_term <- paste0("lag(", y, ", ", i, ")")
      for (l in seq_len(n_stoch)) {
        new_term[l] <- !lag_term %in% dterms[[l]]
        increment[[l]][idx] <- (i %in% k) && new_term[l]
      }
      include[idx] <- any(new_term)
      channels[[idx]] <- dynamitechannel(
        formula = as.formula(paste0(lhs[idx], " ~ ", rhs)),
        family = deterministic_(),
        response = lhs[idx]
      )
    }
  }
  list(
    dformula = parse_new_lags(dformula, channels_stoch, increment, type, lhs),
    channels = channels[include],
    max_lag = max_lag,
    rank = rank[include],
    stoch = stoch[include]
  )
}

#' Parse Manual Lag Terms in a `dynamiteformula` Object
#'
#' @inheritParams parse_lags
#' @param lag_map \[`data.frame`]\cr
#'   Output of `extract_lags`.
#' @param valid_resp \[`character()`]\cr
#'   A vector of valid LHS variables that can  appear in the model formulas.
#' @noRd
parse_singleton_lags <- function(dformula, data, group_var,
                                 lag_map, valid_resp, verbose) {
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
      resp,
      resp_all,
      verbose
    )
    if (y$deterministic) {
      dformula[[y$idx]]$specials$past <- NULL
      dformula[[y$idx]]$specials$past_type <- NULL
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
      if (y$is_resp && !is.null(y$past_val)) {
        if (identical(y$past_type, "past")) {
          past_out <- lag_(y$past_val, i)
          past_out[data[, .I[seq_len(..i)], by = group_var]$V1] <- NA
          spec <- list(
            past = past_out,
            resp_type = y$type
          )
        } else if (identical(y$past_type, "init")) {
          spec <- list(
            past = y$past_val[i],
            resp_type = y$type
          )
        }
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

#' Prepare a New Channel for a Lagged Response
#'
#' @inheritParams parse_singleton_lags
#' @param resp \[`character(1)`]\cr Name of the response variable being lagged
#' @param resp_all  \[`character()`]\cr Vector of all responses in the model.
#' @noRd
prepare_lagged_response <- function(dformula, lag_map,
                                    resp, resp_all, verbose) {
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
    y$past_val <- dformula[[y$idx]]$specials$past
    y$past_idx <- NULL
    y$past_offset <- NULL
    y$deterministic <- is_deterministic(dformula[[y$idx]]$family)
    y$type <- dformula[[y$idx]]$specials$resp_type
    y$past_type <- dformula[[y$idx]]$specials$past_type
    if (y$deterministic && !identical(y$past_type, "past")) {
      y_max_lag <- max(lag_map$k[y$lag_idx])
      y_past_len <- length(y$past_val)
      if (y_past_len < y_max_lag) {
        if (verbose) {
          warning_(c(
            "Deterministic channel {.var {resp}} has a maximum lag of
            {y_max_lag} but you've supplied {cli::no(y_past_len)}
            initial value{?s//s}:",
            `i` = "This may result in NA values for {.var {resp}}."
          ))
        }
        y$past_val <- c(y$past_val, rep(NA, y_max_lag - y_past_len))
      }
    }
  }
  y
}

#' Parse and Add Lags Defined via `lags` to a `dynamiteformula` Object
#'
#' @inheritParams parse_global_lags
#' @param channels_stoch \[`logical()`]\cr A vector indicating which channels
#'   are stochastic.
#' @param increment \[`logical()`]\cr  A vector indicating whether to add
#'   the new lag term or not (e.g.,, whether it was already present or not).
#' @param type \[`character(1)`]\cr Either `"fixed"` or `"varying"`.
#' @param lhs \[`character()`]\cr A vector of the new lagged variable names.
#' @noRd
parse_new_lags <- function(dformula, channels_stoch, increment, type, lhs) {
  for (i in channels_stoch) {
    if (any(increment[[i]])) {
      dformula[[i]] <- dynamiteformula_(
        formula = increment_formula(
          formula = dformula[[i]]$formula,
          x = lhs[increment[[i]]],
          type = type,
          varying_idx = dformula[[i]]$varying,
          varying_icpt = dformula[[i]]$has_varying_intercept,
          fixed_icpt = dformula[[i]]$has_fixed_intercept
        ),
        original = dformula[[i]]$original,
        family = dformula[[i]]$family
      )
    }
  }
  dformula
}

#' Parse Lags That Actually Appear in a `dynamiteformula` Object
#'
#' @inheritParams parse_singleton_lags
#' @param y \[`list()`]\cr Output of `prepare_lagged_response`.
#' @param i \[`integer(1)`]\cr Row index of lag_map.
#' @param lhs \[`character(1)`]\cr Name of the new lagged response.
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

#' Adds NA Gaps to Fill In Missing Time Points in a Data Frame
#'
#' @inheritParams dynamite
#' @noRd
fill_time <- function(data, group_var, time_var) {
  has_groups <- !is.null(group_var)
  time <- sort(unique(data[[time_var]]))
  stopifnot_(
    length(time) > 1L,
    "There must be at least two time points in the data."
  )
  if (has_groups) {
    time_count <- data |>
      dplyr::group_by(.data[[group_var]]) |>
      dplyr::count(.data[[time_var]]) |>
      dplyr::summarise(unique = all(.data[["n"]] == 1L))
    d <- time_count[[group_var]][!time_count$unique]
    stopifnot_(
      all(time_count$unique),
      c(
        "Each time index must correspond to a single observation per group:",
        `x` = "{cli::qty(d)}Group{?s} {.var {d}} of {.var {group_var}}
               {cli::qty(d)}{?has/have} duplicate observations."
      )
    )
  } else {
    stopifnot_(
      all(!duplicated(data[[time_var]])),
      "Each time index must correspond to a single observation."
    )
  }
  time_ivals <- diff(time)
  time_scale <- min(diff(time))
  if (any(time_ivals[!is.na(time_ivals)] %% time_scale > 0)) {
    stop_("Observations must occur at regular time intervals.")
  } else {
    full_time <- seq(time[1L], time[length(time)], by = time_scale)
    if (has_groups) {
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

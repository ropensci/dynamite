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
#' See more details in the package vignette on how to define a dynamite model.
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula. See 'Details'.
#' @param data
#'   \[`data.frame`, `tibble::tibble`, or `data.table::data.table`]\cr
#'   The data frame, tibble or a data.table containing the variables in the
#'   model. Supported column types are `integer`, `logical`, `double`,
#'   `factor`. `character` columns will be converted to factors.
#'   Unused factor levels will be dropped. The `data` can contain missing
#'   values which will simply be ignored in the estimation in a case-wise
#'   fashion (per time-point and per channel).
#' @param group \[`character(1)`]\cr A column name of `data` that denotes the
#'   unique groups.
#' @param time \[`character(1)`]\cr A column name of `data` that denotes the
#'   time axis.
#' @param priors \[`data.frame`]\cr An optional data frame with prior
#'   definitions. See details.
#' @param debug TODO
#' @param ... Additional arguments to [rstan::sampling()].
#' @export
#' @examples
#' \dontrun{
#' fit <- dynamite(obs(y ~ -1 + varying(~x), family = "gaussian" +
#'   lags(type = "varying") + splines(df = 20), gaussian_example, id, time,
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
#' @srrstats {RE1.2} *Regression Software should document expected format (types or classes) for inputting predictor variables, including descriptions of types or classes which are not accepted.*
#' @srrstats {BS5.2} *Bayesian Software should either return the input function or prior distributional specification in the return object; or enable direct access to such via additional functions which accept the return object as single argument.*
#' @srrstats {G2.8} *Software should provide appropriate conversion or dispatch routines as part of initial pre-processing to ensure that all other sub-functions of a package receive inputs of a single defined class or type.*
#' @srrstats {G2.9} *Software should issue diagnostic messages for type conversion in which information is lost (such as conversion of variables from factor to character; standardisation of variable names; or removal of meta-data such as those associated with [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as insertion of variable or column names where none were provided).*
#' @srrstats {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.*
#' @srrstats {G2.13} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#' @srrstats {G2.14} *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
#' @srrstats {G2.14a} *error on missing data*
#' @srrstats {G2.14b} *ignore missing data with default warnings or messages issued*
#' @srrstats {G2.14c} *replace missing data with appropriately imputed values*
#' @srrstats {G2.15} *Functions should never assume non-missingness, and should never pass data with potential missing values to any base routines with default `na.rm = FALSE`-type parameters (such as [`mean()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/mean.html), [`sd()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/sd.html) or [`cor()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html)).*
#' @srrstats {BS1.2} *Description of how to specify prior distributions, both in textual form describing the general principles of specifying prior distributions, along with more applied descriptions and examples, within:*
#' @srrstats {BS1.1} *Descriptions of how to enter data, both in textual form and via code examples. Both of these should consider the simplest cases of single objects representing independent and dependent data, and potentially more complicated cases of multiple independent data inputs.*
#' @srrstats {BS1.3} *Description of all parameters which control the computational process (typically those determining aspects such as numbers and lengths of sampling processes, seeds used to start them, thinning parameters determining post-hoc sampling from simulated values, and convergence criteria). In particular:*
#' @srrstats {BS1.3b} *Where applicable, Bayesian software should document, both in text and examples, how to use different sampling algorithms for a given model.*
#' @srrstats {BS2.6} *Check that values for computational parameters lie within plausible ranges.*
#' @srrstats {BS2.7} *Enable starting values to be explicitly controlled via one or more input parameters, including multiple values for software which implements or enables multiple computational "chains."*
#' @srrstats {BS2.9} *Ensure each chain is started with a different seed by default.*
#' @srrstats {BS2.12} *Bayesian Software should implement at least one parameter controlling the verbosity of output, defaulting to verbose output of all appropriate messages, warnings, errors, and progress indicators.*
#' @srrstats {BS2.13} *Bayesian Software should enable suppression of messages and progress indicators, while retaining verbosity of warnings and errors. This should be tested.*
#' @srrstats {BS2.14} *Bayesian Software should enable suppression of warnings where appropriate. This should be tested.*
#' @srrstats {BS2.15} *Bayesian Software should explicitly enable errors to be caught, and appropriately processed either through conversion to warnings, or otherwise captured in return values. This should be tested.*
#' @srrstats {BS3.0} *Explicitly document assumptions made in regard to missing values; for example that data is assumed to contain no missing (`NA`, `Inf`) values, and that such values, or entire rows including any such values, will be automatically removed from input data.*
#' @srrstats {BS4.0} *Packages should document sampling algorithms (generally via literary citation, or reference to other software)*
#' @srrstats {BS4.5} *Ensure that appropriate mechanisms are provided for models which do not converge.*
#' @srrstats {BS5.0} *Return values should include starting value(s) or seed(s), including values for each sequence where multiple sequences are included*
#' @srrstats {BS5.1} *Return values should include appropriate metadata on types (or classes) and dimensions of input data*
#' @srrstats {BS5.2} *Bayesian Software should either return the input function or prior distributional specification in the return object; or enable direct access to such via additional functions which accept the return object as single argument.*
#' @srrstats {BS5.3} *Bayesian Software should return convergence statistics or equivalent*
#' @srrstats {BS5.5} *Appropriate diagnostic statistics to indicate absence of convergence should either be returned or immediately able to be accessed.*
#' @srrstats {RE1.3a} *Where otherwise relevant information is not transferred, this should be explicitly documented.*
#' @srrstats {RE2.0} *Regression Software should document any transformations applied to input data, for example conversion of label-values to `factor`, and should provide ways to explicitly avoid any default transformations (with error or warning conditions where appropriate).*
#' @srrstats {RE2.1} *Regression Software should implement explicit parameters controlling the processing of missing values, ideally distinguishing `NA` or `NaN` values from `Inf` values (for example, through use of `na.omit()` and related functions from the `stats` package).*
#' @srrstats {RE3.0} *Issue appropriate warnings or other diagnostic messages for models which fail to converge.*
#' @srrstats {RE3.1} *Enable such messages to be optionally suppressed, yet should ensure that the resultant model object nevertheless includes sufficient data to identify lack of convergence.*
#' @srrstats {RE4.0} *Regression Software should return some form of "model" object, generally through using or modifying existing class structures for model objects (such as `lm`, `glm`, or model objects from other packages), or creating a new class of model objects.*
#' @srrstats {RE4.1} *Regression Software may enable an ability to generate a model object without actually fitting values. This may be useful for controlling batch processing of computationally intensive fitting algorithms.*
#' @srrstats {RE4.4} *The specification of the model, generally as a formula (via `formula()`)*
#' @srrstats {RE4.8} *Response variables, and associated "metadata" where applicable.*
#' @srrstats {RE4.13} *Predictor variables, and associated "metadata" where applicable.*
dynamite <- function(dformula, data, group, time,
                     priors = NULL, debug = NULL, ...) {
  stopifnot_(
    is.data.frame(data),
    "Argument {.arg data} is not a {.cls data.frame} object."
  )
  if (missing(group) || is.null(group)) {
    group <- NULL
    group_var <- NULL
  } else {
    group_var <- try_type(group, "character")
    stopifnot_(
      !is.null(data[[group_var]]),
      "Can't find grouping variable {.var {group_var}} in {.arg data}."
    )
  }
  stopifnot_(!missing(time), "Argument {.var time} is missing.")
  time_var <- try_type(time, "character")
  stopifnot_(
    !is.null(data[[time_var]]),
    "Can't find time index variable {.var {time_var}} in {.arg data}."
  )
  data <- parse_data(data, dformula, group_var, time_var)
  dformulas <- parse_lags(data, dformula, group_var, time_var)
  evaluate_deterministic(data, dformulas, group_var, time_var)
  stan <- prepare_stan_data(
    data,
    dformulas$stoch,
    group_var,
    time_var,
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
    rstan::stan_model(model_code = model_code)
  )
  stanfit <- onlyif(
    !isTRUE(debug$no_compile) && !isTRUE(debug$no_sampling),
    rstan::sampling(model, data = stan$sampling_vars, ...)
  )
  out <- structure(
    list(
      stanfit = stanfit,
      dformulas = dformulas,
      data = data,
      stan = stan,
      group_var = group_var,
      time_var = time_var,
      priors = dplyr::bind_rows(stan$priors)
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
#' @srrstats {G2.4d} *explicit conversion to factor via `as.factor()`*
#' @srrstats {G2.5} *Where inputs are expected to be of `factor` type, secondary documentation should explicitly state whether these should be `ordered` or not, and those inputs should provide appropriate error or other routines to ensure inputs follow these expectations.*
#' @srrstats {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
#' @srrstats {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.*
#' @srrstats {G2.16} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.*
#' @noRd
parse_data <- function(data, dformula, group_var, time_var) {
  data <- droplevels(data)
  data_names <- names(data)
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
      is_categorical(dformula[[i]]$family) &&
      all(c("ordered", "factor") %in% class(data[, resp[i]]))
    },
    logical(1L)
  )
  if (any(ordered_factor_resp)) {
    rof <- resp[ordered_factor_resp]
    warning_(c(
      "Response variable{?s} {.var {rof}} {?is/are} of class
       {.cls ordered factor} whose channel{?s} {?is/are} categorical:",
      `i` = "{.var {rof}} will be converted to {?an/} unordered factor{?s}."
    ))
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
  time <- sort(unique(data[[time_var]]))
  stopifnot_(
    length(time) > 1L,
    "There must be at least two time points in the data."
  )
  data <- fill_time(data, time, group_var, time_var)
  data <- data.table::as.data.table(data)
  data.table::setkeyv(data, c(group_var, time_var))
  data
}

#' Parse Lag and Lags Definitions of a `dynamiteformula` Object
#'
#' @inheritParams parse_data
#' @noRd
parse_lags <- function(data, dformula, group_var, time_var) {
  channels_det <- which_deterministic(dformula)
  channels_stoch <- which_stochastic(dformula)
  resp_all <- get_responses(dformula)
  resp_stoch <- resp_all[channels_stoch]
  n_rows <- data[,.N]
  n_channels <- length(resp_all)
  max_lag <- 0L
  for (i in seq_len(n_channels)) {
    dformula[[i]]$formula <- gsub_formula(
      pattern = "lag\\(([^\\,\\)]+)\\)",
      replacement = "lag\\(\\1, 1\\)",
      formula = dformula[[i]]$formula,
      fixed = FALSE,
      perl = TRUE,
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
  list(
    all = dformula,
    det = dformula_det,
    stoch = structure(
      dformula[channels_stoch],
      splines = attr(dformula, "splines")
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
        family = dformula[[j]]$family,
        random_intercept = dformula[[j]]$has_random_intercept
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
    if (lag_map$increment[k]) {
      for (j in seq_along(dformula)) {
        if (y$origin[i, j]) {
          dformula[[j]] <- increment_dynamiteformula(
            dformula[[j]], y, y$term_type[i, j], lhs
          )
        }
      }
    } else {
      for (j in seq_along(dformula)) {
        if (y$origin[i, j]) {
          dformula[[j]]$formula <- gsub_formula(
            pattern = lag_map$src[k],
            replacement = lhs,
            formula = dformula[[j]]$formula,
            fixed = TRUE
          )
        }
      }
    }
  }
  dformula
}

#' Add a new lag term to a dynamiteformula
#'
#' @inheritParams parse_present_lags
#' @param type Type of term to add, either `"fixed"` or `"varying"`.
#' @noRd
increment_dynamiteformula <- function(dformula, y, type, lhs) {
  if (y$deterministic) {
    dformula$formula <- increment_formula_deterministic(dformula$formula, lhs)
  } else {
    dformula <- dynamiteformula_(
      formula = increment_formula(
        formula = dformula$formula,
        x = lhs,
        type = type,
        varying_idx = dformula$varying,
        varying_icpt = dformula$has_varying_intercept,
        fixed_icpt = dformula$has_fixed_intercept
      ),
      original = dformula$original,
      family = dformula$family,
      random_intercept = dformula$has_random_intercept
    )
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
  get_origin(dformula, y)
}

#' Determine which channel a lagged response originated from
#'
#' @inheritParams parse_present_lags
#' @noRd
get_origin <- function(dformula, y) {
  n_resp <- length(dformula)
  n_src <- length(y$src)
  y$origin <- matrix(FALSE, n_src, n_resp)
  y$term_type <- matrix("", n_src, n_resp)
  if (y$deterministic) {
    form_str <- vapply(get_formulas(dformula), deparse1, character(1L))
    for (i in seq_len(n_src)) {
      y$origin[i, ] <- grepl(y$src[i], form_str, fixed = TRUE)
    }
  } else {
    for (j in seq_len(n_resp)) {
      dform_terms <- get_terms(dformula[j])[[1L]]
      for (i in seq_len(n_src)) {
        term_idx <- which(dform_terms == y$src[i])
        if (length(term_idx) > 0L) {
          y$origin[i, j] <- TRUE
          y$term_type[i, j] <- ifelse_(
            term_idx %in% dformula[[j]]$fixed,
            "fixed",
            "varying"
          )
        }
      }
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

#' Estimate a Bayesian Dynamic Multivariate Panel Model
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula. See 'Details'.
#' @param data \[`data.frame`]\cr The data frame containing the variables in
#'   the model.
#' TODO document here allowed data types etc.
#' @param group \[`character(1)`]\cr A column name of `data` that denotes the
#'   unique groups.
#' @param time \[`character(1)`]\cr A column name of `data` that denotes the
#'   time axis.
#' @param priors TODO
#' @param debug TODO
#' @param ... Additional arguments to [rstan::sampling()].
#' @export
#' @examples
#' \dontrun{
#' fit <- dynamite(obs(y ~ -1 + varying(~x), family = gaussian()) +
#'   lags("varying") + splines(df = 20), gaussian_example, id, time,
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
# TODO all
# TODO document priors
# TODO document what missingness means
# TODO warn ordered factor as response
dynamite <- function(dformula, data, group, time,
                     priors = NULL, debug = NULL, ...) {
  # stored for return object
  if (!is.data.frame(data)) {
    stop_("Argument {.var data} is not a {.cls data.frame} object.")
  }
  if (missing(group) || is.null(group)) {
    group <- group_var <- NULL
  } else {
    group_var <- try_(group = group, type = "character")
    if (is.null(data[[group_var]])) {
      stop_(
        "Can't find grouping variable {.var {group_var}} in {.var data}."
      )
    }
  }
  if (missing(time)) {
    stop_("Argument {.var time} is missing.")
  }
  time_var <- try_(time = time, type = "character")
  if (is.null(data[[time_var]])) {
    stop_(
      "Can't find time index variable {.var {time_var}} in {.var data}."
    )
  }
  data <- parse_data(data, dformula, group_var, time_var)
  dformulas <- parse_lags(data, dformula, group_var, time_var)
  evaluate_deterministic(data, dformulas, group_var, time_var)
  # TODO check for NAs
  stan <- prepare_stan_data(data, dformulas$stoch, group_var, time_var, priors,
                            fixed = attr(dformulas$all, "max_lag"))
  model_code <- create_blocks(dformula = dformulas$stoch, indent = 2L,
                              vars = stan$model_vars)
  # TODO needs to be NULL?
  # model_vars[grep("_prior_distr_", names(model_vars))] <- NULL
  # debug <- dots$debug
  model <- if (!is.null(debug) && isTRUE(debug$no_compile)) {
    NULL
  } else {
    rlang::inform("Compiling Stan model.")
    rstan::stan_model(model_code = model_code)
  }
  stanfit <- if (isTRUE(debug$no_compile) || isTRUE(debug$no_sampling)) {
    NULL
  } else {
    rstan::sampling(model, data = stan$sampling_vars, ...)
  }
  # TODO return the function call for potential update method?
  out <- structure(
    list(
      # TODO what else do we need to return?
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
      if (!"try-error" %in% class(got)) {
        out[[opt]] <- got
      }
    }
  }
  out
}

#' Access model formula of a dynamite model fit
#'
#' Return a list containing the formulas defining each channel of the model.
#'
#' @param x \[`dynamitefit`\]\cr The model fit object
#' @param ... Not used
#'
#' @export
formula.dynamitefit <- function(x, ...) {
  get_originals(x$dformulas$all)
}

#' Parse data for model fitting
#'
#' @param data \[`data.frame`]\cr The data frame containing the variables in
#'   the model.
#' @param dformula \[`dynamiteformula`] The model formula
#' @param group_var \[`character(1)`,`NULL`]\cr Group variable name or `NULL`
#'   if there is only one group
#' @param time_var \[`character(1)`]\cr Time index variable name
#'
#' @srrstats {G2.4d} *explicit conversion to factor via `as.factor()`*
#' @srrstats {G2.5} *Where inputs are expected to be of `factor` type, secondary documentation should explicitly state whether these should be `ordered` or not, and those inputs should provide appropriate error or other routines to ensure inputs follow these expectations.*
#' @srrstats {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
#' @srrstats {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.*
#' @srrstats {G2.16} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.*
parse_data <- function(data, dformula, group_var, time_var) {
  data <- droplevels(data) # TODO document this in return value
  data_names <- names(data)
  data <- data |>
    dplyr::mutate(dplyr::across(where(is.character), as.factor))
  valid_types <- c("integer", "logical", "double")
  col_types <- sapply(data, typeof)
  factor_cols <- sapply(data, is.factor)
  valid_cols <- (col_types %in% valid_types) | factor_cols
  if (any(!valid_cols)) {
    invalid_cols <- data_names[!valid_cols]
    invalid_types <- col_types[!valid_cols]
    stop_(c(
      "Column{?s} {.var {invalid_cols}} of {.var data} {?is/are} invalid:",
      `x` = "Column type{?s} {.cls {invalid_types}} {?is/are} not supported."
    ))
  }
  coerce_cols <- valid_cols & !factor_cols
  if (any(coerce_cols)) {
    for (i in which(coerce_cols)) {
      data[,i] <- do.call(paste0("as.", typeof(data[,i])),
                          args = list(data[,i]))
    }
  }
  resp <- get_responses(dformula)
  ordered_factor_resp <- sapply(seq_along(resp), function(i) {
    is_categorical(dformula[[i]]$family) &&
      all(c("ordered", "factor") %in% class(data[, resp[i]]))
  })
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
  finite_cols <- sapply(data, function(x) all(is.finite(x) | is.na(x)))
  if (any(!finite_cols)) {
    stop_(
      "Non-finite values in variable{?s} {.var {data_names[!finite_cols]}} of
      {.var data}."
    )
  }
  time <- sort(unique(data[[time_var]]))
  if (length(time) == 1) {
    stop_("There must be at least two time points in the data.")
  }
  full_time <- NULL
  # TODO convert Dates etc. to integers before this
  time_ivals <- diff(time)
  time_scale <- min(diff(time))
  if (any(time_ivals[!is.na(time_ivals)] %% time_scale > 0)) {
    stop_("Observations must occur at regular time intervals.")
  } else {
    full_time <- seq(time[1], time[length(time)], by = time_scale)
    groups <- !is.null(group_var)
    if (groups) {
      time_groups <- data |>
        dplyr::group_by(.data[[group_var]]) |>
        dplyr::summarise(has_missing = !identical(.data[[time_var]], full_time))
      if (any(time_groups$has_missing)) {
        full_data_template <- expand.grid(
          time = time,
          group = unique(data[[group_var]])
        )
        names(full_data_template) <- c(time_var, group_var)
        data <- full_data_template |>
          dplyr::left_join(data, by = c(group_var, time_var))
      }
    } else {
      if (!identical(data[[time_var]], full_time)) {
        full_data_template <- expand.grid(
          time = time
        )
        names(full_data_template) <- time_var
        data <- full_data_template |>
          dplyr::left_join(data, by = time_var)
      }
    }
  }
  data <- data.table::as.data.table(data)
  data.table::setkeyv(data, c(group_var, time_var))
  data
}

parse_lags <- function(data, dformula, group_var, time_var) {
  channels_det <- which_deterministic(dformula)
  channels_stoch <- which_stochastic(dformula)
  resp_all <- get_responses(dformula)
  resp_stoch <- resp_all[channels_stoch]
  n_rows <- data[,.N]
  n_channels <- length(resp_all)
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
  predictors <- get_predictors(dformula)
  terms_stoch <- unique(unlist(get_terms(dformula[channels_stoch])))
  non_lags <- extract_nonlags(terms_stoch)
  mis_vars <- which(!non_lags %in% c(resp_all, data_names))
  if (length(mis_vars) > 0) {
    stop_(
      "Can't find variable{?s} {.var {non_lags[mis_vars]}} in {.var data}."
    )
  }
  lag_map <- extract_lags(predictors)
  lag_ext <- attr(dformula, "lags")
  lags_channel <- list()
  lags_rank <- integer(0)
  lags_lhs <- character(0)
  lags_stoch <- logical(0)
  lags_max <- 0
  if (!is.null(lag_ext)) {
    idx <- 0
    lag_k <- lag_ext$k
    lag_type <- lag_ext$type
    lags_max <- max(lag_k)
    n_lag <- lags_max * length(resp_stoch)
    lags_channel <- vector(mode = "list", length = n_lag)
    lags_stoch <- logical(n_lag)
    lags_rank <- integer(n_lag)
    lags_lhs <- character(n_lag)
    lags_rank <- integer(n_lag)
    lags_increment <- logical(n_lag)
    for (i in seq_len(lags_max)) {
      for (j in seq_along(channels_stoch)) {
        y <- resp_stoch[j]
        idx <- idx + 1
        lags_lhs[idx] <- paste0(y, "_lag", i)
        if (i == 1) {
          lags_rhs <- y
          lags_stoch[idx] <- TRUE
        } else {
          lags_rhs <- paste0(y, "_lag", i - 1)
        }
        lags_rank[idx] <- i
        lags_increment[idx] <- i %in% lag_k
        lags_channel[[idx]] <- dynamitechannel(
          formula = as.formula(paste0(lags_lhs[idx], " ~ ", lags_rhs)),
          family = deterministic_(),
          response = lags_lhs[idx]
        )
      }
    }
    for (j in channels_stoch) {
      dformula[[j]] <- dynamiteformula_(
        formula = increment_formula(
          formula = dformula[[j]]$formula,
          x = lags_lhs[lags_increment],
          type = lag_type,
          varying_idx = dformula[[j]]$varying,
          varying_icpt = dformula[[j]]$has_varying_intercept,
          fixed_icpt = dformula[[j]]$has_fixed_intercept
        ),
        family = dformula[[j]]$family,
        random_intercept = dformula[[j]]$has_random_intercept
      )
    }
    lag_map <- lag_map |>
      dplyr::filter(!(.data$var %in% resp_all & .data$k %in% lag_ext$k))
  }
  mis_vars <- which(!lag_map$var %in% c(resp_all, data_names))
  if (length(mis_vars) > 0) {
    stop_(c(
      "Unable to construct lagged values of {.var {cs(lag_map$var[mis_vars])}}:",
      `x` = "Can't find such variables in {.var data}."
    ))
  }
  n_lag <- lag_map |> nrow()
  map_channel <- list()
  map_resp <- character(0)
  map_rank <- integer(0)
  map_pred <- logical(0)
  map_stoch <- logical(0)
  map_stoch_k <- lag_map |>
    dplyr::filter(.data$var %in% resp_stoch) |>
    dplyr::pull(.data$k)
  max_lag <- lags_max
  if (length(map_stoch_k) > 0) {
    max_lag <- max(lags_max, map_stoch_k)
  }
  if (n_lag > 0) {
    if (any(lag_map$k <= 0)) {
      stop_("Shift values must be positive in {.fun lag}.")
    }
    map_channel <- vector(mode = "list", length = n_lag)
    map_resp <- character(n_lag)
    map_pred <- logical(n_lag)
    map_stoch <- logical(n_lag)
    map_rank <- integer(n_lag)
    lag_resp <- unique(lag_map$var)
    idx <- 0
    for (y in lag_resp) {
      lag_idx <- which(lag_map$var == y)
      y_idx <- which(resp_all == y)
      y_resp <- length(y_idx) > 0
      y_deterministic <- TRUE
      if (y_resp) {
        y_past <- NULL
        y_past_idx <- NULL
        y_past_offset <- NULL
        y_obs <- NULL
        y_self <- NULL
        y_stoch <- TRUE
        y_deterministic <- is_deterministic(dformula[[y_idx]]$family)
        y_type <- dformula[[y_idx]]$specials$resp_type
        if (y_deterministic) {
          y_form <- deparse1(dformula[[y_idx]]$formula)
          y_past <- dformula[[y_idx]]$specials$past
          y_past_len <- length(y_past)
          y_self <- max(extract_self_lags(y_form, y))
          y_max <- max(lag_map$k[lag_idx])
          y_obs <- extract_lags(y_form) |>
            dplyr::filter(.data$var %in% c(resp_stoch, data_names)) |>
            dplyr::pull(.data$k)
          if (length(y_obs) > 0) {
            y_obs <- max(y_obs)
          } else {
            y_obs <- 0
          }
          # TODO better variable names
          y_req_past <- min(y_max - max_lag + max(y_self, y_obs), y_max)
          y_past_offset <- y_max - y_req_past
          if (y_past_len < y_req_past) {
            stop_(c(
              "Deterministic channel {.var {y}} requires {y_req_past} initial
               value{?s}:",
              `x` = "You've supplied {cli::no(y_past_len)} initial value{?s//s}."
            ))
          }
          y_stoch <- FALSE
          y_past_idx <- 0
          dformula[[y_idx]]$specials$past <- NULL
        }
      }
      for (i in seq_along(lag_idx)) {
        idx <- idx + 1
        j <- lag_idx[i]
        if (i == 1) {
          map_rhs <- y
        } else {
          map_rhs <- paste0(y, "_lag", i - 1)
        }
        map_lhs <- paste0(y, "_lag", i)
        map_rank[idx] <- i
        map_stoch[idx] <- !y_deterministic
        map_pred[idx] <- !y_resp
        map_resp[idx] <- y
        spec <- NULL
        if (y_resp) {
          if (!is.null(y_past)) {
            if (i > y_past_offset) {
              y_past_idx <- y_past_idx + 1L
              spec <- list(past = y_past[y_past_idx],
                           past_offset = max_lag,
                           resp_type = y_type)
            }
          }
        }
        map_channel[[idx]] <- dynamitechannel(
          formula = as.formula(paste0(map_lhs, " ~ ", map_rhs)),
          family = deterministic_(),
          response = map_lhs,
          specials = spec
        )
        if (lag_map$present[j]) {
          for (k in seq_len(n_channels)) {
            dformula[[k]]$formula <- gsub_formula(
              pattern = lag_map$src[j],
              replacement = map_lhs,
              formula = dformula[[k]]$formula,
              fixed = TRUE
            )
          }
        }
      }
    }
  }
  dformula_det <- dformula[channels_det]
  dformula_lag_pred <- map_channel[map_pred]
  dformula_lag_stoch <- c(lags_channel[lags_stoch],
                          map_channel[map_stoch & !map_pred])
  dformula_lag_det <- c(lags_channel[!lags_stoch],
                        map_channel[!map_stoch & !map_pred])
  attr(dformula_lag_pred, "rank_order") <- order(map_rank[map_pred])
  attr(dformula_lag_det, "rank_order") <-
    order(c(lags_rank[!lags_stoch], map_rank[!map_stoch & !map_pred]))
  attr(dformula_lag_pred, "original_response") <- map_resp[map_pred]
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

evaluate_deterministic <- function(data, dformulas, group_var, time_var) {
  fixed <- as.integer(attr(dformulas$all, "max_lag"))
  n_time <- length(unique(data[[time_var]]))
  n_id <- 1L
  if (!is.null(group_var)) {
    n_id <- length(unique(data[[group_var]]))
  }
  dd <- dformulas$det
  dlp <- dformulas$lag_pred
  dld <- dformulas$lag_det
  dls <- dformulas$lag_stoch
  n_det <- length(dd)
  n_lag_det <- length(dld)
  n_lag_stoch <- length(dls)
  cl <- get_quoted(dd)
  initialize_deterministic(data, dd, dlp, dld, dls)
  idx <- seq.int(1L, n_time * n_id, by = n_time) - 1L
  assign_initial_values(data, dd, dlp, dld, dls, idx, fixed, group_var)
  if (n_time > fixed + 1L) {
    ro_stoch <- 1:n_lag_stoch
    ro_det <- attr(dld, "rank_order")
    lhs_det <- get_responses(dld)
    rhs_det <- get_predictors(dld)
    lhs_stoch <- get_responses(dls)
    rhs_stoch <- get_predictors(dls)
    idx <- idx + fixed + 1L
    for (i in seq.int(fixed + 2L, n_time)) {
      idx <- idx + 1L
      if (n_lag_det > 0) {
        assign_lags(data, ro_det, idx, lhs_det, rhs_det)
      }
      if (n_lag_stoch > 0) {
        assign_lags(data, ro_stoch, idx, lhs_stoch, rhs_stoch)
      }
      if (n_det > 0) {
        assign_deterministic(data, cl, idx)
      }
    }
  }
}

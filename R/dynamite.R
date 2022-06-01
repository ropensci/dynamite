#' Estimate a Bayesian Dynamic Multivariate Panel Model
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula. See 'Details'.
#' @param data \[`data.frame`]\cr The data frame containing the variables in
#'   the model.
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
dynamite <- function(dformula, data, group, time,
                     priors = NULL, debug = NULL, ...) {
  # stored for return object
  if (!is.data.frame(data)) {
    stop_("Argument 'data' is not a data.frame object")
  }
  if (missing(group) || is.null(group)) {
    group <- group_var <- NULL
  } else {
    group_var <- try_(group = group, type = "character")
    if (is.null(data[[group_var]])) {
      stop_("Grouping variable '", group_var, "' is not present in the data")
    }
  }
  if (missing(time)) {
    stop_("Argument 'time' is missing.")
  }
  time_var <- try_(time = time, type = "character")
  if (is.null(data[[time_var]])) {
    stop_("Time index variable '", time_var, "' is not present in the data")
  }
  original_dformula <- dformula
  data <- parse_data(data, group_var, time_var)
  dformulas <- parse_lags(data, dformula, group_var, time_var)
  evaluate_deterministic(data, dformulas$det, dformulas$lag, group_var, time_var)
  # TODO check for NAs
  stan <- prepare_stan_data(data, dformulas$stoch, group_var, time_var, priors)
  model_code <- create_blocks(dformula = dformulas$stoch, indent = 2L,
                              vars = stan$model_vars)
  # TODO needs to be NULL?
  # model_vars[grep("_prior_distr_", names(model_vars))] <- NULL
  # debug <- dots$debug
  model <- if (!is.null(debug) && isTRUE(debug$no_compile)) {
    NULL
  } else {
    message("Compiling Stan model")
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
      stanfit = stanfit,
      # TODO what else do we need to return?
      dformula = original_dformula,
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

parse_data <- function(data, group_var, time_var) {

  data <- droplevels(data) # TODO document this in return value
  data <- data |>
    dplyr::mutate(dplyr::across(where(is.character), as.factor))
  time <- sort(unique(data[[time_var]]))
  if (length(time) == 1) {
    stop_("There must be at least two time points in the data")
  }
  full_time <- NULL
  # TODO convert Dates etc. to integers before this
  time_ivals <- diff(time)
  time_scale <- min(diff(time))
  if (any(time_ivals[!is.na(time_ivals)] %% time_scale > 0)) {
    stop_("Observations must occur at regular time intervals")
  } else {
    full_time <- seq(time[1], time[length(time)], by = time_scale)
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
  }
  data.table::setDT(data, key = c(group_var, time_var))
  data
}

parse_lags <- function(data, dformula, group_var, time_var) {
  channels_det <- which_deterministic(dformula)
  channels_stoch <- which_stochastic(dformula)
  resp_all <- get_responses(dformula)
  resp_stoch <- resp_all[channels_stoch]
  n_rows <- data[,.N]
  n_channels <- length(resp_all)
  data_names <- names(data)
  predictors <- get_predictors(dformula)
  lag_map <- extract_lags(predictors)
  lag_ext <- attr(dformula, "lags")
  lags_channel <- list()
  lags_rank <- integer(0)
  lags_lhs <- NULL
  lags_max <- 0
  if (!is.null(lag_ext)) {
    lag_k <- lag_ext$k
    lag_type <- lag_ext$type
    lags_max <- max(lag_k)
    n_lag <- lags_max * length(resp_stoch)
    lags_channel <- vector(mode = "list", length = n_lag)
    lags_rank <- integer(n_lag)
    lags_lhs <- character(n_lag)
    lags_rank <- integer(n_lag)
    lags_increment <- logical(n_lag)
    for (i in seq_len(lags_max)) {
      for (j in seq_along(channels_stoch)) {
        y <- resp_stoch[j]
        idx <- (i - 1) * n_lag + j
        lags_lhs[idx] <- paste0(y, "_lag_", i)
        if (i == 1) {
          lags_rhs <- y
        } else {
          lags_rhs <- paste0(y, "_lag_", i - 1)
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
    for (j in seq_len(n_channels)) {
      dformula[[j]] <- dynamiteformula_(
        formula = increment_formula(
          formula = dformula[[j]]$formula,
          x = lags_lhs[lags_increment],
          type = lag_type,
          varying_idx = dformula[[j]]$varying
        ),
        family = dformula[[j]]$family
      )
    }
    lag_map <- lag_map |>
      dplyr::filter(!(.data$var %in% resp_all & .data$k %in% lag_ext$k))
  }
  mis_vars <- which(!lag_map$var %in% c(resp_all, data_names))
  if (length(mis_vars)) {
    stop_("Unable to construct lagged values of '",
          cs(lag_map$var[mis_vars]), "', ",
          "no such variables are present in the data")
  }
  n_lags <- lag_map |> dplyr::filter(.data$var %in% resp_all) |> nrow()
  map_lhs <- NULL
  map_channel <- list()
  map_rank <- integer(0)
  max_lag <- max(lags_max, lag_map$k)
  if (n_lags > 0) {
    if (any(lag_map$k <= 0)) {
      stop_("Only positive shift values are allowed in lag()")
    }
    map_channel <- vector(mode = "list", length = n_lags)
    map_rank <- integer(n_lags)
    lag_resp <- unique(lag_map$var)
    idx <- 0
    for (y in lag_resp) {
      lag_idx <- which(lag_map$var == y)
      y_idx <- which(resp_all == y)
      if (length(y_idx) == 1) {
        y_past <- NULL
        y_past_idx <- NULL
        y_obs_lags <- NULL
        y_stoch <- TRUE
        y_deterministic <- is_deterministic(dformula[[y_idx]]$family)
        if (y_deterministic) {
          y_past <- dformula[[y_idx]]$specials$past
          y_past_len <- length(y_past)
          y_self_lags <- max(lag_map$k[lag_idx])
          if (y_past_len < y_self_lags) {
            stop_("Deterministic channel '", y, "' requires ", y_self_lags, " ",
                  "initial values, but only ", y_past_len, " values ",
                  "have been specified")
          }
          y_stoch <- FALSE
          y_past_idx <- 0
          dformula[[y_idx]]$specials$past <- NULL
        }
        for (i in seq_along(lag_idx)) {
          j <- lag_idx[i]
          y_past_idx <- y_past_idx + 1
          if (i == 1) {
            map_rhs <- y
          } else {
            map_rhs <- paste0(y, "_lag", i - 1)
          }
          map_lhs <- paste0(y, "_lag", i)
          idx <- idx + 1
          map_rank[idx] <- i
          if (!is.null(y_past)) {
            spec <- list(past = y_past[y_past_idx],
                         past_offset = max_lag)
          } else {
            spec <- NULL
          }
          if (!y_deterministic) {
            data[, (map_lhs) := lapply(.SD, lag_, k = i),
                 .SDcols = y, by = group_var]
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
      } else {
        data_idx <- which(data_names == y)
        for (i in seq_along(lag_idx)) {
          j <- lag_idx[i]
          if (lag_map$present[j]) {
            lag_k <- lag_map$k[j]
            lag_var <- paste0(data_names[data_idx], "_lag", lag_k)
            data[, (lag_var) := lapply(.SD, lag_, k = lag_k),
                 .SDcols = y, by = group_var]
            for (k in seq_len(n_channels)) {
              dformula[[k]]$formula <- gsub_formula(
                pattern = lag_map$src[j],
                replacement = lag_var,
                formula = dformula[[k]]$formula,
                fixed = TRUE
              )
            }
          }
        }
      }
    }
  }
  dformula_det <- dformula[channels_det]
  dformula_lag <- c(lags_channel, map_channel)
  attr(dformula_lag, "rank_order") <- order(c(lags_rank, map_rank))
  attr(dformula_lag, "max_lag") <- max_lag
  list(
    all = dformula,
    det = dformula_det,
    lag = dformula_lag,
    stoch = structure(
      dformula[channels_stoch],
      splines = attr(dformula, "splines")
    )
  )
}

evaluate_deterministic <- function(data, dd, dl, group_var, time_var) {
  n_time <- length(unique(data[[time_var]]))
  n_id <- length(unique(data[[group_var]]))
  fixed <- as.integer(attr(dl, "max_lag"))
  idx <- seq.int(1L, n_time * n_id, by = n_time)
  assign_initial_values(data, dd, dl, idx)
  n_det <- length(dd)
  n_lag <- length(dl)
  if (n_time > fixed && (n_det > 0 || n_lag > 0)) {
    cl <- get_quoted(dd)
    ro <- attr(dl, "rank_order")
    idx <- idx + fixed
    lag_lhs <- get_responses(dl)
    lag_rhs <- get_predictors(dl)
    if (n_det) {
      assign_deterministic(data, cl, idx)
    }
    for (i in seq.int(fixed + 2L, n_time)) {
      idx <- idx + 1L
      if (n_lag) {
        assign_lags(data, ro, idx, lag_lhs, lag_rhs)
      }
      if (n_det) {
        assign_deterministic(data, cl, idx)
      }
    }
  }
}

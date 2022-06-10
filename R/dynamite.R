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
  original_dformula <- dformula
  data <- parse_data(data, group_var, time_var)
  dformulas <- parse_lags(data, dformula, group_var, time_var)
  evaluate_deterministic(data, dformulas, group_var, time_var)
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
        y_obs_lags <- NULL
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
          if (length(y_obs_lags) > 0) {
            y_obs <- max(y_obs_lags)
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
            if (i > y_past_offset)
              y_past_idx <- y_past_idx + 1
              spec <- list(past = y_past[y_past_idx],
                           past_offset = max_lag,
                           resp_type = y_type)
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

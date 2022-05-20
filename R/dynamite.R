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
#' @param ... TODO
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
  original_dformula <- dformula
  if (!is.data.frame(data)) {
    stop_("Argument 'data' is not a data.frame object")
  }
  e <- new.env()
  e$data <- data
  vars <- parse_data(e, group, time)
  dformulas <- parse_lags(e, dformula, vars$group, vars$time)
  evaluate_deterministic(e, dformulas$det, vars$group, vars$time)
  stan <- prepare_stan_data(e, dformulas$stoch, vars$group, vars$time, priors)
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
      time_var = vars$time,
      group_var = vars$group,
      priors = dplyr::bind_rows(stan$priors)
      #spline = list(
      #  B = sampling_vars$Bs,
      #  D = sampling_vars$D
      #),
      #ord = data_names[!data_names %in% c(group_var, time_var)],
      #J = attr(model_matrix, "assign")
    ),
    class = "dynamitefit"
  )
  # Adds any object in the environment of this function to the return object
  # if its name is icluded in the debug argument
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

parse_data <- function(e, group, time) {
  if (missing(group) || is.null(group)) {
    group <- group_var <- NULL
  } else {
    group_var <- try_(group = group, type = "character")
    if (is.null(e$data[[group_var]])) {
      stop_("Grouping variable '", group_var, "' is not present in the data")
    }
  }
  if (missing(time)) {
    stop_("Argument 'time' is missing.")
  }
  time_var <- try_(time = time, type = "character")
  if (is.null(e$data[[time_var]])) {
    stop_("Time index variable '", time_var, "' is not present in the data")
  }
  e$data <- droplevels(e$data) # TODO document this in return value
  e$data <- e$data |>
    dplyr::mutate(dplyr::across(where(is.character), as.factor)) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(group_var, time_var))))
  time <- sort(unique(e$data[[time_var]]))
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
    time_groups <- e$data |>
      dplyr::group_by(.data[[group_var]]) |>
      dplyr::summarise(has_missing =
                         !identical(.data[[time_var]], full_time))
    if (any(time_groups$has_missing)) {
      full_data_template <- expand.grid(
        time = time,
        group = unique(e$data[[group_var]])
      )
      names(full_data_template) <- c(time_var, group_var)
      e$data <- full_data_template |>
        dplyr::left_join(e$data, by = c(group_var, time_var))
    }
  }
  list(group = group_var, time = time_var)
}

parse_lags <- function(e, dformula, group_var, time_var) {
  channels_det <- which_deterministic(dformula)
  channels_stoch <- which_stochastic(dformula)
  resp_all <- get_responses(dformula)
  n_rows <- nrow(e$data)
  n_channels <- length(resp_all)
  data_names <- names(e$data)
  lag_map <- extract_lags(get_predictors(dformula))
  lag_all <- attr(dformula, "lags")
  lags_channel <- list()
  lags_lhs <- NULL
  if (!is.null(lag_all)) {
    n_lag <- lag_all$k * n_channels
    lags_channel <- vector(mode = "list", length = n_lag)
    type <- lag_all$type
    lags_lhs <- character(n_lag)
    lags_rank <- integer(n_lag)
    for (i in seq_len(lag_all$k)) {
      for (j in seq_along(channels_stoch)) {
        k <- channels_stoch[j]
        idx <- (i - 1) * n_channels + j
        lags_lhs[idx] <- paste0(resp_all[k], "_lag_", i)
        if (i == 1) {
          lags_rhs <- paste0("lag_(", resp_all[k], ", 1)")
        } else {
          lags_rhs <- paste0("lag_(", resp_all[k], "_lag_", i - 1, ", 1)")
        }
        lags_channel[[idx]] <- dynamitechannel(
          formula = as.formula(paste0(lags_lhs[idx], " ~ ", lags_rhs)),
          family = deterministic_(),
          response = lags_lhs[idx],
          specials = list(rank = i)
        )
        attr(lags_channel[[idx]], "stoch_origin") <- (i == 1)
      }
    }
    for (j in seq_len(n_channels)) {
      if (identical(type, "varying")) {
        dformula[[j]] <- dynamiteformula_(
          formula = increment_varying(
            formula = dformula[[j]]$formula,
            x = lags_lhs,
            varying_idx = dformula[[j]]$varying
          ),
          family = dformula[[j]]$family
        )
      } else {
        dformula[[j]] <- dynamiteformula_(
          formula = increment_fixed(
            formula = dformula[[j]]$formula,
            x = lags_lhs
          ),
          family = dformula[[j]]$family
        )
      }
    }
    if (nrow(lag_map) > 0) {
      lag_map <- lag_map[!(lag_map$var %in% resp_all & lag_map$k <= lag_all$k), ]
    }
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
  if (n_lags > 0) {
    if (any(lag_map$k <= 0)) {
      stop_("Only positive shift values are allowed in lag()")
    }
    map_channel <- vector(mode = "list", length = n_lags)
    lag_resp <- unique(lag_map$var)
    idx <- 0
    for (y in lag_resp) {
      lag_idx <- which(lag_map$var == y)
      y_idx <- which(resp_all == y)
      if (length(y_idx) == 1) {
        y_past_offset <- 0
        y_past <- NULL
        y_stoch <- TRUE
        if (is_deterministic(dformula[[y_idx]]$family)) {
          y_past <- dformula[[y_idx]]$specials$past
          y_past_len <- length(y_past)
          lag_max <- max(lag_map$k[lag_idx])
          if (y_past_len < lag_max) {
            stop_("Deterministic channel '", y, "' requires ", lag_max, " ",
                  "initial values, but only ", y_past_len, " values ",
                  "have been specified")
          }
          if (y_past_len > lag_max) {
            y_past_offset <- 1
          }
          y_stoch <- FALSE
          proc_past <- seq(from = y_past_len, by = -1, length.out = lag_max)
          dformula[[y_idx]]$specials$past <-
            dformula[[y_idx]]$specials$past[-proc_past]
        }
        for (i in seq_along(lag_idx)) {
          j <- lag_idx[i]
          if (i == 1) {
            map_rhs <- paste0("lag_(", y, ", 1)")
          } else {
            map_rhs <- paste0("lag_(", y, "_lag", i - 1, ", 1)")
          }
          map_lhs <- paste0(y, "_lag", i)
          idx <- idx + 1
          map_channel[[idx]] <- dynamitechannel(
            formula = as.formula(paste0(map_lhs, " ~ ", map_rhs)),
            family = deterministic_(),
            response = map_lhs,
            specials = list(past = y_past[i + y_past_offset], rank = i)
          )
          attr(map_channel[[idx]], "stoch_origin") <- y_stoch && (i == 1)
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
            e$data <- e$data |>
              dplyr::group_by(.data[[group_var]]) |>
              dplyr::mutate({{lag_var}} := lag_(.data[[y]], lag_k)) |>
              dplyr::ungroup()
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
  dformula_det <- c(dformula[channels_det], lags_channel, map_channel)
  if (length(dformula_det) > 0) {
    rank <- get_ranks(dformula_det)
    u_rank <- sort(unique(rank))
    rank_list <- lapply(u_rank, function(x) which(rank == x))
    attr(dformula_det, "rank_list") <- rank_list
  }
  list(
    all = dformula,
    det = dformula_det,
    stoch = structure(
      dformula[channels_stoch],
      splines = attr(dformula, "splines")
    )
  )
}


# TODO Maybe having a continuous range of non-NA values for
# each individual is useful for some special case?
# data_mis <- data[ ,c(group_var, time_var)]
# data_mis$obs <- complete.cases(data)
# time_segments <- data_mis |>
#     dplyr::group_by(dplyr::across(dplyr::all_of(c(group_var)))) |>
#     dplyr::summarise(last_obs = which(obs)[sum(obs)],
#                      valid_missingness_pattern_ =
#                          all(obs[1:last_obs] == cummin(obs[1:last_obs])) ||
#                          all(obs[1:last_obs] == cummax(obs[1:last_obs])))
# if (any(!time_segments$valid_missingness_pattern_)) {
#     # TODO is there a better term or a way to convey this?
#     stop_("Observed time series must not contain gaps.")
# }

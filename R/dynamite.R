#' Estimate a Bayesian Dynamic Multivariate Panel Model
#'
#' @param formula \[`dynamiteformula`]\cr The model formula. See 'Details'.
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
dynamite <- function(formula, data, group, time,
                     priors = NULL, debug = NULL, ...) {
  # dots <- list(...) #Note: Use explicit debug-argument as otherwise it is passed to sampling with an error
  data <- droplevels(data) # TODO document this in return value
  if (missing(group)) {
    group <- group_var <- NULL
  } else {
    group_var <- substitute(group)
    if (!is.character(group_var)) {
      group_var <- deparse(group_var)
    }
    if (is.null(data[[group_var]])) {
      stop_("Grouping variable '", group_var, "' is not present in the data")
    }
  }
  if (missing(time)) {
    stop_("Argument 'time' is missing.")
  }
  time_var <- substitute(time)
  if (!is.character(time_var)) {
    time_var <- deparse(time_var)
  }
  if (is.null(data[[time_var]])) {
    stop_("Time index variable '", time_var, "' is not present in the data")
  }
  data <- data |>
    dplyr::mutate(dplyr::across(where(is.character), as.factor)) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(group_var, time_var))))
  time <- sort(unique(data[[time_var]]))
  full_time <- NULL
  # TODO convert Dates etc. to integers before this
  time_ivals <- diff(time)
  time_scale <- min(diff(time))
  if (any(time_ivals[!is.na(time_ivals)] %% time_scale > 0)) {
    stop_("Observations must occur at regular time intervals.")
  } else {
    full_time <- seq(time[1], time[length(time)], by = time_scale)
    time_groups <- data |>
      dplyr::group_by(!!as.symbol(group_var)) |>
      dplyr::summarise(has_missing =
                         !identical(!!as.symbol(time_var), full_time))
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
  group <- data[[group_var]]
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
  formula_det <- get_deterministic(formula)
  formula_stoch <- get_stochastic(formula)
  resp_det <- get_resp(formula_det)
  resp_all <- get_resp(formula)
  n_rows <- nrow(data)
  n_resp <- length(resp_all)
  data_names <- names(data)
  lag_map <- extract_lags(get_pred(formula))
  lag_all <- attr(formula, "lags")
  lags_uneval <- rep(TRUE, nrow(lag_map))
  lags_lhs <- NULL
  lags_rhs <- NULL
  #lags_def <- NULL
  lags_rank <- NULL
  if (!is.null(lag_all)) {
    type <- lag_all$type
    n_lag <- lag_all$k * n_resp
    lags_lhs <- character(n_lag)
    lags_rhs <- character(n_lag)
    #lags_def <- character(n_lag)
    lags_rank <- integer(n_lag)
    for (i in seq_len(lag_all$k)) {
      for (j in seq_along(formula_stoch)) {
        idx <- (i - 1) * n_resp + j
        lags_lhs[idx] <- paste0(resp_all[j], "_lag_", i)
        if (i == 1) {
          lags_rhs[idx] <- paste0("lag_(", resp_all[j], ", 1)")
        } else {
          lags_rhs[idx] <- paste0("lag_(", resp_all[j], "_lag_", i - 1, ", 1)")
        }
        lags_rank[idx] <- i
        #lags_def[idx] <- paste0(lags_lhs[idx], " ~ ", lags_rhs[idx])
      }
    }
    for (j in seq_len(n_resp)) {
      c(formula[[j]]$predictors) <- lags_lhs
      c(formula[[j]][[type]]) <- which(formula[[j]]$predictors %in% lags_rhs)
    }
    lags_uneval[lag_map$def %in% resp_all & lag_map$k <= lag_all$k] <- FALSE
  }
  map_lhs <- NULL
  map_def <- NULL
  map_rank <- NULL
  if (any(lags_uneval)) {
    idx_uneval <- which(lags_uneval)
    n_uneval <- sum(lag_map$k[idx_uneval])
    map_lhs <- character(sum(n_uneval))
    map_rhs <- character(sum(n_uneval))
    map_def <- character(sum(n_uneval))
    map_rank <- integer(sum(n_uneval))
    idx <- 0
    for (i in idx_uneval) {
      if (lag_map$k[i] <= 0) {
        stop_("Only positive shift values are allowed in lag()")
      }
      if (is_as_is(lag_map$src[i])) {
        lag_map$scr[i] <- gsub("lag", "lag_", lag_map$src[i])
      }
      if (is_as_is(lag_map$def[i])) {
        lag_map$def[i] <- gsub("I\\((.*)\\)", "\\1", lag_map$def[i])
      }
      for (j in 1:lag_map$k[i]) {
        idx <- idx + 1
        if (j == 1) {
          map_rhs[idx] <- paste0("lag_(", lag_map$def[i], ", ", lag_map$k[i], ")")
        } else {
          map_rhs[idx] <- paste0("lag_(", lag_map$def[i], "_lag_", lag_map$k[i] - 1, ", 1)")
        }
        map_lhs[idx] <- paste0(lag_map$def[i], "_lag_", lag_map$k[i])
        map_rank[idx] <- j
        #map_def[idx] <- paste0(map_lhs[idx], " ~ ", map_rhs)
        if (j == lag_map$k[i]) {
          for (l in seq_len(n_resp)) {
            formula[[l]]$formula <- gsub_formula(
              pattern = lag_map$src[i],
              replacement = map_lhs[idx],
              formula = formula[[l]]$formula,
              fixed = TRUE
            )
          }
        }
      }
    }
  }
  # for (j in seq_len(n_resp)) {
  #   if (length(formula[[j]]$predictors) > 0) {
  #     formula[[j]]$formula <- reformulate(
  #       termlabels = formula[[j]]$predictors,
  #       response = resp_all[j],
  #       intercept = attr(terms(formula[[j]]$formula), "intercept")
  #     )
  #   } else {
  #     formula[[j]]$formula <- as.formula(paste0(resp_all[j], "~ 1"))
  #   }
  # }
  # EVALUATION ORDER
  # 1st: initial values of deterministic channels
  # 2nd: lag-terms in order of precedence
  # 3rd: deterministic channels

  #deter_init <- has_initial(formula_deter)
  #if (length(deter_init) > 0) {
  #  #TODO can this be dplyrized?
  #  for (i in deter_init) {
  #    for (g in seq_along(group)) {
  #      data[g,]
  #    }
  #  }
  #}
  n_time <- length(full_time)
  n_id <- length(unique(group))
  id_offset <- seq(0, n_time * (n_id - 1), by = n_time)
  # Fist time index
  det_init <- has_initial(formula_det)
  if (any(!det_init)) {
    data[1 + id_offset, resp_det] <-
      full_model.matrix_pseudo(get_form(formula_det[!det_init]), data[1 + id_offset, ])
  }
  if (any(det_init)) {
    data[1 + id_offset, resp_det[det_init]] <-
      full_model.matrix_pseudo(get_initial(formula_det[det_init]), data[1 + id_offset, ])
  }

  # Lags (both deterministic and stochastic)
  lhs_new <- c(lags_lhs, map_lhs)
  n_new <- length(lhs_new)
  if (n_new > 0) {
    formula_new <- vector(mode = "list", length = n_new)
    u_lhs <- which(!duplicated(lhs_new))
    lhs_new <- c(lags_lhs, map_lhs)[u_lhs]
    rhs_new <- c(lags_rhs, map_rhs)[u_lhs]
    rank_new <- c(lags_rank, map_rank)[u_lhs]
    idx <- 0
    for (r in 1:max(rank_new)) {
      new_idx <- which(rank_new == r)
      new_start <- idx + 1
      for (i in new_idx) {
        idx <- idx + 1
        formula_new[[idx]] <- as.formula(paste0(" ~ ", rhs_new[i]))
      }
      new_end <- idx
      if (n_time > 1) {
        id_offset_vec <- rep(id_offset, each = 2)
        model_idx <- seq(2, 2 * n_id, by = 2)
        for (i in 2:n_time) {
          past_idx <- (i-1):i + id_offset_vec
          data_idx <- i + id_offset
          # TODO deterministic here!!
          data[data_idx, lhs_new[new_idx]] <-
            full_model.matrix_pseudo(formula_new[new_start:new_end],
                                     data[past_idx, ])[model_idx,]
        }
      }
    }
  }
  # Deterministic channels
  if (n_time > 1) {
    for (i in 2:n_time) {
      data[i + id_offset, resp_det] <-
        full_model.matrix_fast(get_form(formula_det), data[i + id_offset, ])
    }
  }
  responses <- data[, resp_all, drop = FALSE]
  # Needs sapply/lapply instead of apply to keep factors as factors
  attr(responses, "resp_class") <- sapply(responses, class)
  model_matrix <- full_model.matrix(get_form(formula), data)
  resp_levels <- lapply(responses, levels)
  # TODO: simplify I(lag(variable, 1)) to something shorter, e.g. lag_1(variable)?
  # TODO: shorten variable and or channel name if they are very long?
  # TODO: NOTE! you can use lag_map to get the variables within the complicated definition for formatting
  u_names <- colnames(model_matrix)
  coef_names <- lapply(seq_along(resp_all), function(i) {
    x <- paste0(resp_all[i], "_", u_names[attr(model_matrix, "assign")[[i]]])
    if (is_categorical(formula[[i]]$family)) {
      levels_ <- resp_levels[[i]][-length(resp_levels[[i]])]
      # for prior names, there's probably more elegant way...
      #Need to keep in mind as.data.frame function, and fixed vs varying
      simplified <- list(names = x, levels = levels_)
      x <- paste0(x, "_", rep(levels_, each = length(x)))
      attr(x, "simplified") <- simplified
    }
    x
  })
  specials <- evaluate_specials(formula, data)
  converted <- convert_data(formula, responses, specials, group,
                            full_time, model_matrix, coef_names, priors)
  model_vars <- converted$model_vars
  sampling_vars <- converted$sampling_vars
  model_priors <- converted$priors
  model_code <- create_blocks(formula, indent = 2L, vars = model_vars)
  model_vars[grep("_prior_distr_", names(model_vars))] <- NULL
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
    rstan::sampling(model, data = sampling_vars, ...)
  }
  # TODO return the function call for potential update method?
  out <- structure(
    list(
      stanfit = stanfit,
      coef_names = coef_names,
      # TODO what else do we need to return?
      time = time,
      time_var = time_var,
      group_var = group_var,
      levels = resp_levels,
      specials = specials,
      # TODO: extract only D for as.data.frame and J for predict
      model_vars = model_vars,
      data = data,
      spline = list(
        B = sampling_vars$Bs,
        D = sampling_vars$D
      ),
      priors = dplyr::bind_rows(model_priors),
      prediction_basis = list(
        formula = formula,
        # fixed = fixed,
        # past = model_matrix[(n_rows - fixed):n_rows,],
        # start = model_matrix[1:fixed,], # Needed for some posterior predictive checks?
        ord = data_names[!data_names %in% c(group_var, time_var)],
        J = attr(model_matrix, "assign")
      )
    ),
    class = "dynamitefit"
  )
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

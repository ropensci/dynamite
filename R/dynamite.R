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
  channels_det <- which_deterministic(formula)
  channels_stoch <- which_stochastic(formula)
  #resp_det <- get_responses(formula_det)
  resp_all <- get_responses(formula)
  #form_det <- get_formulas(formula_det)
  n_rows <- nrow(data)
  n_channels <- length(resp_all)
  data_names <- names(data)
  lag_map <- extract_lags(get_predictors(formula))
  lag_all <- attr(formula, "lags")
  lags_uneval <- rep(TRUE, nrow(lag_map))
  lags_channel <- list()
  lags_lhs <- NULL
  #lags_rhs <- NULL
  # <- NULL
  if (!is.null(lag_all)) {
    n_lag <- lag_all$k *
    lags_channel <- vector(mode = "list", length = n_lag)
    type <- lag_all$type
    lags_lhs <- character(n_lag)
    #lags_rhs <- character(n_lag)
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
        #lags_rank[idx] <- i
        lags_channel[[idx]] <- dynamitechannel(
          formula = as.formula(paste0(lags_lhs[idx], " ~ ", lags_rhs)),
          family = deterministic_(),
          response = lags_lhs[idx],
        )
        attr(lags_channel[[idx]], "stoch_origin") <- (i == 1)
      }
    }
    for (j in seq_len(n_channels)) {
      c(formula[[j]]$predictors) <- lags_lhs
      c(formula[[j]][[type]]) <- which(formula[[j]]$predictors %in% lags_rhs)
    }
    lags_uneval[lag_map$def %in% resp_all & lag_map$k <= lag_all$k] <- FALSE
  }
  map_lhs <- NULL
  map_rank <- NULL
  if (any(lags_uneval)) {
    idx_uneval <- which(lags_uneval)
    n_uneval <- sum(lag_map$k[idx_uneval])
    map_channel <- vector(mode = "list", length = n_uneval)
    idx <- 0
    for (i in idx_uneval) {
      if (lag_map$k[i] <= 0) {
        stop_("Only positive shift values are allowed in lag()")
      }
      y <- lag_map$def[i]
      # if (is_as_is(lag_map$src[i])) {
      #   lag_map$scr[i] <- gsub("lag", "lag_", lag_map$src[i])
      # }
      # if (is_as_is(lag_map$def[i])) {
      #   lag_map$def[i] <- gsub("I\\((.*)\\)", "\\1", lag_map$def[i])
      # }
      y_idx <- which(resp_all == y)
      y_past <- NULL
      y_stoch <- TRUE
      if (is_deterministic(formula[[y_idx]]$family)) {
        y_past <- formula[[y_idx]]$specials$past
        y_past_len <- length(y_past)
        if (y_past_len < lag_map$k[i]) {
          stop_("Deterministic channel '", y, "' requires ", lag_map$k[i], " ",
                "initial values, but only ", y_past_len, " values ",
                "have been specified")
        }
        y_stoch <- FALSE
        proc_past <- seq(from = y_past_len, by = -1, length.out = lag_map$k[i])
        formula[[y_idx]]$specials$past <-
          formula[[y_idx]]$specials$past[-proc_past]
      }
      for (j in 1:lag_map$k[i]) {
        idx <- idx + 1
        if (j == 1) {
          map_rhs <- paste0("lag_(", y, ", ", j, ")")
        } else {
          map_rhs <- paste0("lag_(", y, "_lag_", j - 1, ", 1)")
        }
        map_lhs <- paste0(y, "_lag_", j)
        map_channel[[idx]] <- dynamitechannel(
          formula = as.formula(paste0(map_lhs, " ~ ", map_rhs)),
          family = deterministic_(),
          response = map_lhs,
          specials = list(past = y_past[j], rank = j)
        )
        attr(map_channel[[idx]], "stoch_origin") <- y_stoch && (j == 1)
        for (l in seq_len(n_channels)) {
          formula[[l]]$formula <- gsub_formula(
            pattern = lag_map$src[i],
            replacement = map_lhs,
            formula = formula[[l]]$formula,
            fixed = TRUE
          )
        }
      }
    }
  }
  formula_det <- c(formula[channels_det], lags_channel, map_channel)
  resp_det <- get_responses(formula_det)
  rank_det <- get_ranks(formula_det)
  n_time <- length(full_time)
  n_id <- length(unique(group))
  id_offset <- seq(0, n_time * (n_id - 1), by = n_time)
  det_init <- has_past(formula_det)
  if (any(det_init)) {
    idx <- which(det_init)
    for (i in idx) {
      data[1 + id_offset, formula_det[[i]]$response] <-
        formula_det[[i]]$specials$past
    }
  }
  if (any(!det_init)) {
    det_from_stoch <- sapply(formula_det, function(y){
      isTRUE(attr(y, "stoch_origin"))
    })
    det_a <- !det_init & det_from_stoch
    det_b <- !det_init & !det_from_stoch
    if (any(det_a)) {
      data_eval <- data[rep(1 + id_offset, each = 2),]
      eval_idx <- seq(1, n_id, by = 2)
      data_eval[eval_idx,] <- NA
      data[1 + id_offset, resp_det[det_a]] <-
        full_model.matrix_pseudo(get_formulas(formula_det[det_a]),
                                 data_eval)[eval_idx + 1, ]
    }
    if (any(det_b)) {
      data[1 + id_offset, resp_det[det_b]] <-
        full_model.matrix_pseudo(get_formulas(formula_det[det_b]),
                                 data[1 + id_offset, ])
    }
  }
  if (length(formula_det) > 0 && n_time > 1) {
    id_offset_vec <- rep(id_offset, each = 2)
    model_idx <- seq(2, 2 * n_id, by = 2)
    u_rank <- sort(unique(rank_det))
    for (i in 2:n_time) {
      past_idx <- (i - 1):i + id_offset_vec
      data_idx <- i + id_offset
      for (j in u_rank) {
        idx <- which(rank_det == j)
        data[data_idx, get_responses(formula_det[idx])] <-
          full_model.matrix_pseudo(get_formulas(formula_det[idx]),
                                   data[past_idx, ])[model_idx,]
      }
    }
  }
  responses <- data[, resp_all, drop = FALSE]
  # Needs sapply/lapply instead of apply to keep factors as factors
  attr(responses, "resp_class") <- sapply(responses, class)
  model_matrix <- full_model.matrix(formula[channels_stoch], data)
  resp_levels <- lapply(responses, levels)
  u_names <- colnames(model_matrix)
  coef_names <- lapply(seq_along(channels_stoch), function(i) {
    j <- channels_stoch[j]
    x <- paste0(resp_all[j], "_", u_names[attr(model_matrix, "assign")[[i]]])
    if (is_categorical(formula[[j]]$family)) {
      levels_ <- resp_levels[[j]][-length(resp_levels[[j]])]
      # for prior names, there's probably more elegant way...
      #Need to keep in mind as.data.frame function, and fixed vs varying
      simplified <- list(names = x, levels = levels_)
      x <- paste0(x, "_", rep(levels_, each = length(x)))
      attr(x, "simplified") <- simplified
    }
    x
  })
  specials <- evaluate_specials(formula[channels_stoch], data)
  converted <- convert_data(formula[channels_stoch], responses, specials, group,
                            full_time, model_matrix, coef_names, priors)
  model_vars <- converted$model_vars
  sampling_vars <- converted$sampling_vars
  model_priors <- converted$priors
  model_code <- create_blocks(formula[channels_stoch], indent = 2L, vars = model_vars)
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

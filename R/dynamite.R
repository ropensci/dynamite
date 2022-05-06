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
  }
  if (missing(time)) {
    stop_("Argument 'time' is missing.")
  }
  time_var <- substitute(time)
  if (!is.character(time_var)) {
    time_var <- deparse(time_var)
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
      dplyr::summarise(has_missing = !identical({{ time_var }}, full_time))
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
  resp_all <- get_resp(formula)
  n_rows <- nrow(data)
  n_resp <- length(resp_all)
  lag_map <- extract_lags(unlist(get_pred(formula)))
  unprocessed_lags <- rep(TRUE, nrow(lag_map))
  data_names <- names(data)
  if (!is.null(lag_all <- attr(formula, "lags"))) {
    type <- lag_all$type
    pred_lag <- character(lag_all$k * n_resp)
    for (i in seq_len(lag_all$k)) {
      for (j in seq_len(n_resp)) {
        ix <- (i - 1) * n_resp + j
        pred_lag[ix] <- paste0("I(lag_(", resp_all[j], ", ", i, "))")
      }
    }
    for (j in seq_len(n_resp)) {
      c(formula[[j]]$predictors) <- pred_lag
      c(formula[[j]][[type]]) <- which(formula[[j]]$predictors %in% pred_lag)
    }
    unprocessed_lags[lag_map$def %in% resp_all & lag_map$k <= lag_all$k] <- FALSE
  }
  if (any(unprocessed_lags)) {
    for (i in which(unprocessed_lags)) {
      if (lag_map$k[i] <= 0) {
        stop_("Only positive shift values are allowed in lag().")
      }
      if (is_as_is(lag_map$src[i])) {
        lag_map$scr[i] <- gsub("lag", "lag_", lag_map$src[i])
      }
      if (is_as_is(lag_map$def[i])) {
        lag_map$def[i] <- gsub("I\\((.*)\\)", "\\1", lag_map$def[i])
      }
      pred_lag <- paste0("I(lag_(", lag_map$def[i], ", ", lag_map$k[i], "))")
      for (j in seq_len(n_resp)) {
        formula[[j]]$predictors <- gsub(lag_map$src[i], pred_lag,
                                        formula[[j]]$predictors, fixed = TRUE)
      }
    }
  }
  for (j in seq_len(n_resp)) {
    formula[[j]]$formula <- reformulate(
      termlabels = formula[[j]]$predictors,
      response = resp_all[j],
      intercept = attr(terms(formula[[j]]$formula), "intercept")
    )
  }
  responses <- data[, resp_all, drop = FALSE]
  # Needs sapply/lapply instead of apply to keep factors as factors
  attr(responses, "resp_class") <- sapply(responses, class)
  model_matrix <- full_model.matrix(formula, data)
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

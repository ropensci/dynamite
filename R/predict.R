#' Predict method for a Bayesian Time-Varying Coefficients Model
#'
#' @param object Object of class `dynamitefit`.
#' @param newdata Data frame used in predictions. If `NULL` (default), the
#'   data used in model estimation is used for predictions as well.
#' @param mode TODO: Think about the names, if we predict missing observations
#'   does the term counterfactual make sense?
#' @param type Type of prediction, `"response"` (default), `"mean"`
#'   or `"link"`.
#' @param n_draws Number of posterior samples to use, default is all.
#' @param n_fixed Number of observations per individual that should be assumed
#'  fixed and will not be predicted
#' @param ... Ignored.
#' @export
predict.dynamitefit <- function(object, newdata = NULL,
                                mode = c("counterfactual", "forecast"),
                                type = c("response", "mean", "link"),
                                n_draws = NULL, n_fixed = NULL, ...) {
  mode <- match.arg(mode)
  type <- match.arg(type)
  do.call(paste0("predict.dynamitefit_", mode),
          list(object = object, newdata = newdata,
               type = type, n_draws = n_draws, n_fixed = n_fixed))
}

predict.dynamitefit_counterfactual <- function(object, newdata, type,
                                               n_draws, n_fixed) {
  if (is.null(n_draws)) {
    n_draws <- ndraws(object)
  }
  fixed <- object$dformulas$max_lag
  if (is.null(n_fixed)) {
    n_fixed <- fixed
  }
  # TODO check that n_draws is positive and <= actual number of draws
  # TODO impute predictor values
  group_var <- object$group_var
  time_var <- object$time_var
  if (is.null(newdata)) {
    newdata_null <- TRUE
  }
  newdata <- check_newdata(newdata, object$data, group_var, time_var)
  group <- unique(newdata[[group_var]])
  time <- unique(newdata[[time_var]])
  n_time <- length(time)
  n_id <- length(group)
  n_new <- nrow(newdata)
  resp_stoch <- get_responses(object$dformulas$stoch)
  resp_det <- get_responses(object$dformulas$det)
  for (resp in resp_stoch) {
    if (is.null(newdata[[resp]])) {
      stop_("Response variable '", resp, "' not found in 'newdata'")
    }
  }
  if (n_fixed < fixed) {
    stop_("The model implies at least ", fixed, " fixed time points, ",
          "but only ", n_fixed, " were specified")
  }
  non_na <- newdata |>
    dplyr::group_by(.data[[group_var]]) |>
    dplyr::summarise(obs = stats::complete.cases(
      dplyr::across(
        dplyr::all_of(resp_stoch)
      )
    ), .groups = "keep")
  fixed_obs <- non_na |>
    dplyr::summarise(
      first_obs = which(.data$obs)[1],
      horizon = all(.data$obs[.data$first_obs:(.data$first_obs + n_fixed - 1)])
    )
  lacking_obs <- is.na(fixed_obs$horizon) | (fixed_obs$horizon < n_fixed)
  if (any(lacking_obs)) {
    groups_lacking <- unique(fixed_obs[lacking_obs, group_var])
    stop_("Insufficient non-NA observations in groups: ", cs(groups_lacking))
  }
  if (newdata_null) {
    predict_idx <- unlist(lapply(seq_len(n_id), function(i) {
      (fixed_obs$first_obs[i] + n_fixed):n_time + (i - 1) * n_time
    }))
    newdata[predict_idx, c(resp_stoch, resp_det)] <- NA
  }
  specials <- evaluate_specials(object$dformulas$stoch, newdata)
  if (type != "response") {
    # create separate column for each level of categorical variables
    for (i in seq_along(resp_stoch)) {
      resp <- resp_stoch[i]
      if (is_categorical(object$dformulas$stoch[[i]]$family)) {
        resp_levels <- object$levels[[resp]]
        newdata[, (glue::glue("{resp}_{resp_levels}")) := NA_integer_] # TODO: glued names to formula?
      } else {
        newdata[, (glue::glue("{resp}_store")) := newdata[[resp]]]
      }
    }
  }
  k <- n_id * n_draws
  newdata <- data.table::as.data.table(newdata)
  newdata <- newdata[rep(seq_len(n_new), n_draws), ]
  newdata[, ("draw") := rep(1:n_draws, each = n_new)]
  data.table::setkeyv(newdata, c("draw", group_var, time_var))
  n <- nrow(newdata)
  newdata_names <- names(newdata)
  full_group <- c("draw", group_var)
  assign_initial_values(newdata, object$dformulas$det, full_group)
  exprs <- lapply(get_formulas(object$dformulas$det), function(x) {
    formula2expression(formula_rhs(x))
  })
  ro <- attr(object$dformulas$det, "rank_order")
  samples <- rstan::extract(object$stanfit)
  J <- lapply(seq_along(resp_stoch), function(j) {
    object$stan$model_vars[[j]]$J
  })
  J_fixed <- lapply(seq_along(resp_stoch), function(j) {
    object$stan$model_vars[[j]]$J_fixed
  })
  J_varying <- lapply(seq_along(resp_stoch), function(j) {
    object$stan$model_vars[[j]]$J_varying
  })
  for (i in 2:n_time) {
    idx <- rep(seq(i - 1, n, by = n_time), each = 2) + rep(0:1, times = k)
    idx_i <- seq(i, n, by = n_time)
    assign_deterministic(newdata, resp_det, exprs, ro, idx, idx_i, full_group)
    model_matrix <- full_model.matrix_predict(
      object$dformulas$stoch,
      newdata[idx, ],
      object$stan$u_names
    )
    for (j in seq_along(resp_stoch)) {
      resp <- resp_stoch[j]
      na_i <- which(is.na(newdata[idx_i, .SD, .SDcols = resp]))
      idx_pred <- idx_i[na_i]
      if (length(na_i) > 0) {
        sim <- do.call(
          paste0("predict_", object$dformulas$stoch[[j]]$family),
          list(
            model_matrix = model_matrix[na_i, J[[j]], drop = FALSE],
            samples = samples, specials[[j]], resp,
            time = i - 1, type, n_draws = n_draws,
            J_fixed = J_fixed[[j]], J_varying = J_varying[[j]]
          )
        )
        if (is_categorical(object$dformulas$stoch[[j]]$family)) {
          resp_levels <- object$levels[[resp]]
          if (type != "response") {
            idx_resp <- which(newdata_names %in%
                                c(glue::glue("{resp}_{resp_levels}")))
            newdata[idx_pred, (newdata_names[idx_resp]) := sim$mean_or_link]
          }
          newdata[idx_pred, (resp) := resp_levels[sim$response]]
        } else {
          if (type != "response") {
            newdata[idx_pred, (glue::glue("{resp}_store")) := sim$mean_or_link]
          }
          newdata[idx_pred, (resp) := sim$response]
        }
      }
    }
  }
  if (type != "response") {
    for (i in seq_along(resp_stoch)) {
      resp <- resp_stoch[i]
      if (is_categorical(object$dformulas$stoch[[i]]$family)) {
        newdata[[resp]] <- NULL
      } else {
        newdata[[resp]] <- newdata[[c(glue::glue("{resp}_store"))]]
        newdata[[c(glue::glue("{resp}_store"))]] <- NULL
      }
    }
  }
  # TODO return either full newdata or just the responses
  # newdata[, resp_all, drop = FALSE]
  newdata
}

predict.dynamitefit_forecast <- function(object, newdata, type, n_draws) {
  stop_("Forecasting is not yet supported")
}

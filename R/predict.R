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
#'  fixed and will not be predicted (counting from the first non-NA observation)
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
  fixed <- attr(object$dformulas$lag, "max_lag")
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
  formulas_stoch <- get_formulas(object$dformulas$stoch)
  resp_stoch <- get_responses(object$dformulas$stoch)
  resp_det <- get_responses(object$dformulas$det)
  lag_lhs <- get_responses(object$dformulas$lag)
  lag_rhs <- get_predictors(object$dformulas$lag)
  n_det <- length(resp_det)
  n_lag <- length(lag_lhs)
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
    newdata[predict_idx, c(resp_stoch, resp_det, lag_lhs)] <- NA
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
  idx <- seq.int(1L, n, by = n_time)
  idx_par <- rep(1L:n_draws, each = n_id)
  assign_initial_values(newdata, object$dformulas$det, object$formulas$lag, idx)
  samples <- rstan::extract(object$stanfit)
  model_vars <- object$stan$model_vars
  n_resp <- length(resp_stoch)
  eval_envs <- vector(mode = "list", length = n_resp)
  for (j in seq_along(resp_stoch)) {
    resp <- resp_stoch[j]
    resp_family <- object$dformulas$stoch[[j]]$family
    eval_envs[[j]] <- new.env()
    eval_envs[[j]]$type <- type
    eval_envs[[j]]$newdata <- newdata # reference
    eval_envs[[j]]$J_fixed <- model_vars[[j]]$J_fixed
    eval_envs[[j]]$J_varying <- model_vars[[j]]$J_varying
    eval_envs[[j]]$K_fixed <- model_vars[[j]]$K_fixed
    eval_envs[[j]]$K_varying <- model_vars[[j]]$K_varying
    eval_envs[[j]]$k <- k
    eval_envs[[j]]$resp <- resp_stoch[j]
    eval_envs[[j]]$phi <- samples[[paste0("phi_", resp)]][idx_par]
    eval_envs[[j]]$sigma <- samples[[paste0("sigma_", resp)]][idx_par]
    eval_envs[[j]]$offset <- specials[[j]]$offset
    eval_envs[[j]]$trials <- specials[[j]]$trials
    if (is_categorical(resp_family)) {
      resp_levels <-  attr(class(object$stan$responses[,resp]), "levels")
      if (model_vars[[j]]$has_fixed_intercept) {
        eval_envs[[j]]$alpha <-
          samples[[paste0("alpha_", resp)]][idx_par, , drop=FALSE]
      } else if (model_vars[[j]]$has_varying_intercept) {
        eval_envs[[j]]$alpha <-
          samples[[paste0("alpha_", resp)]][idx_par, , , drop=FALSE]
      }
      eval_envs[[j]]$beta <-
        samples[[paste0("beta_", resp_stoch[j])]][idx_par, , , drop = FALSE]
      eval_envs[[j]]$delta <-
        samples[[paste0("delta_", resp_stoch[j])]][idx_par, , , , drop = FALSE]
      eval_envs[[j]]$nu <-
        samples[[paste0("nu_", resp)]][idx_par, , drop = FALSE]
      eval_envs[[j]]$resp_levels <- resp_levels
      eval_envs[[j]]$S <- dim(samples[[paste0("beta_", resp)]])[4] + 1
      eval_envs[[j]]$idx_resp <- which(newdata_names %in%
                          c(glue::glue("{resp}_{resp_levels}")))
    } else {
      if (model_vars[[j]]$has_fixed_intercept) {
        eval_envs[[j]]$alpha <-
          samples[[paste0("alpha_", resp)]][idx_par, drop = FALSE]
      } else if (model_vars[[j]]$has_varying_intercept) {
        eval_envs[[j]]$alpha <-
          samples[[paste0("alpha_", resp)]][idx_par, , drop = FALSE]
      }
      eval_envs[[j]]$beta <-
        samples[[paste0("beta_", resp_stoch[j])]][idx_par, , drop = FALSE]
      eval_envs[[j]]$delta <-
        samples[[paste0("delta_", resp_stoch[j])]][idx_par, , , drop = FALSE]
      eval_envs[[j]]$nu <- samples[[paste0("nu_", resp)]][idx_par]
    }
    eval_envs[[j]]$call <-
      generate_predict_call(resp, resp_family,
                            model_vars[[j]]$has_fixed,
                            model_vars[[j]]$has_varying,
                            model_vars[[j]]$has_fixed_intercept,
                            model_vars[[j]]$has_varying_intercept,
                            model_vars[[j]]$has_random_intercept,
                            model_vars[[j]]$has_offset)
  }
  cl <- get_quoted(object$dformulas$lag)
  ro <- attr(object$dformulas$lag, "rank_order")
  for (i in 2L:n_time) {
    idx <- idx + 1L
    if (n_lag > 0) {
      assign_lags(newdata, ro, idx, lag_lhs, lag_rhs)
    }
    if (n_det) {
      assign_deterministic(newdata, cl, idx)
    }
    model_matrix <- full_model.matrix_predict(
      formulas_stoch,
      newdata,
      idx,
      object$stan$u_names
    )
    for (j in seq_along(resp_stoch)) {
      e <- eval_envs[[j]]
      idx_na <- is.na(newdata[idx, .SD, .SDcols = resp])
      e$idx <- idx
      e$time <- i - 1
      e$idx_pred <- idx[which(idx_na)]
      e$model_matrix <- model_matrix
      e$a_time <- ifelse_(NCOL(e$alpha) == 1, 1, i - 1)
      if (any(idx_na)) {
        eval(e$call, envir = e)
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
  # for consistency with other output types
  as.data.frame(newdata)
}

predict.dynamitefit_forecast <- function(object, newdata, type, n_draws) {
  stop_("Forecasting is not yet supported")
}

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
#' @param ... Ignored.
#' @export
predict.dynamitefit <- function(object, newdata = NULL,
                                mode = c("counterfactual", "forecast"),
                                type = c("response", "mean", "link"),
                                n_draws = NULL, ...) {
  mode <- match.arg(mode)
  type <- match.arg(type)
  do.call(paste0("predict.dynamitefit_", mode),
          list(object = object, newdata = newdata,
               type = type, n_draws = n_draws))
}

predict.dynamitefit_counterfactual <- function(object, newdata,
                                               type, n_draws = NULL) {
  if (is.null(n_draws)) {
    n_draws <- ndraws(object)
  }
  # TODO needs to support different starting value? Although nothing is done for the time points without NA so perhaps not necessary?

  if (is.null(newdata)) {
    newdata <- object$data |>
      dplyr::arrange(
        dplyr::across(
          dplyr::all_of(c(object$group_var, object$time_var))
        )
      )
    group <- unique(newdata[[object$group_var]])
    time <- unique(newdata[[object$time_var]])
  } else {
    # TODO newdata must start from full_time[1]
    if (!(object$group_var %in% names(newdata))) {
      stop_("Grouping variable '", object$group_var, "' not found in 'newdata'")
    }
    group <- newdata[[object$group_var]]
    if (is.factor(group)) {
      # TODO is this necessary? only length of unique values matters
      group <- droplevels(group)
    }
    group <- unique(group)
    # TODO doesn't really matter at least at the moment
    if (!all(group %in% object$data[[object$group_var]])) {
      stop_("Grouping variable '", object$group_var, "' ",
            "contains new levels not found in the original data")
    }
    if (!(object$time_var %in% names(newdata))) {
      stop_("Time index variable '", object$time_var, "' ",
            "not found in 'newdata'")
    }
    newdata <- newdata |>
      dplyr::arrange(
        dplyr::across(
          dplyr::all_of(c(object$group_var, object$time_var))
        )
      )
    # TODO just use the original time points starting from start_time
    time <- unique(newdata[[object$time_var]])
    if (!all(time %in% object$time)) {
      stop_("Timing variable '", object$time_var, "' ",
            "contains time points not found in the original data")
    }
  }
  n_time <- length(time)
  n_id <- length(group)
  resp_stoch <- get_responses(object$dformulas$stoch)
  resp_det <- get_responses(object$dformulas$det)
  for (resp in resp_stoch) {
    if (is.null(newdata[[resp]])) {
      newdata[[resp]] <- NA
    }
  }
  specials <- evaluate_specials(object$dformulas$stoch, newdata)
  if (type != "response") {
    # create separate column for each level of categorical variables
    for (i in seq_along(resp_stoch)) {
      resp <- resp_stoch[i]
      if (is_categorical(object$dformulas$stoch[[i]]$family)) {
        resp_levels <- object$levels[[resp]]
        newdata[, c(glue::glue("{resp}_{resp_levels}"))] <- NA # TODO: glued names to formula?
      } else {
        newdata[[c(glue::glue("{resp}_store"))]] <- newdata[[resp]]
      }
    }
  }
  k <- n_id * n_draws
  e <- new.env()
  e$data <- data.frame(newdata, draw = rep(1:n_draws, each = nrow(newdata)))
  n <- nrow(e$data)
  assign_initial_values(e, object$dformulas$det, seq(1, n, by = n_time),
                        n_time, n_id)
  samples <- rstan::extract(object$stanfit)
  for (i in 2:n_time) {
    idx <- rep(seq(i - 1, n, by = n_time), each = 2) + rep(0:1, times = k)
    idx_i <- seq(i, n, by = n_time)
    model_matrix <- full_model.matrix_predict(
      object$dformulas$stoch,
      e$data[idx, ],
      object$u_names
    )
    for (j in seq_along(resp_stoch)) {
      resp <- resp_stoch[j]
      J <- object$stan$model_vars[[j]]$J
      if (any(is.na(e$data[idx_i, resp]))) { # TODO partial missingness?
        sim <- do.call(
          paste0("predict_", object$dformulas$stoch[[j]]$family),
          list(
            model_matrix = model_matrix[, J[[j]], drop = FALSE],
            samples = samples, specials[[j]], resp,
            time = i - 1, type, n_draws = n_draws
          )
        )
        if (is_categorical(object$dformulas$stoch[[j]]$family)) {
          resp_levels <- object$levels[[resp]]
          if (type != "response") {
            idx_resp <- which(names(newdata) %in%
                                c(glue::glue("{resp}_{resp_levels}")))
            e$data[idx_i, idx_resp] <- sim$mean_or_link
          }
          e$data[idx_i, resp] <- resp_levels[sim$response]
        } else {
          if (type != "response") {
            e$data[idx_i, c(glue::glue("{resp}_store"))] <- sim$mean_or_link
          }
          e$data[idx_i, resp] <- sim$response
        }
      }
    }
    assign_deterministic(e, object$dformulas$det, idx, idx_i, k)
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
  e$data
}

predict.dynamitefit_forecast <- function(object, newdata, type, n_draws) {
  stop_("Forecasting is not yet supported")
}

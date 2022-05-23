
#' Extract fitted values of dynamitefit
#'
#' TODO Note that these are conditional on the observed data i.e., we don't
#' simulate new lagged values for covariates, so we underestimate the uncertainty.
#' Better to use predict with type = "mean" TODO do it
#' @export
#' @param object An object of class \code{dynamitefit}.
#' @param newdata TODO
#' @param n_draws TODO
#' @param ... Ignored.
fitted.dynamitefit <- function(object, newdata = NULL, n_draws = NULL, ...) {
  if (is.null(n_draws)) {
    n_draws <- ndraws(object)
  }
  group_var <- object$group_var
  time_var <- object$time_var
  newdata <- check_newdata(newdata, object$data, group_var, time_var)
  group <- unique(newdata[[group_var]])
  time <- unique(newdata[[time_var]])
  n_time <- length(time)
  n_id <- length(group)
  resp_stoch <- get_responses(object$dformulas$stoch)
  resp_det <- get_responses(object$dformulas$det)
  fixed <- object$dformulas$lag_max
  n_time <- length(time)
  n_id <- length(group)
  if (n_time <= fixed) {
    stop_("Model definition implies ", fixed, " fixed time points, ",
          "but 'newdata' has only ", n_time, " time points.")
  }
  resp_all <- get_responses(basis$dformula)
  samples <- rstan::extract(object$stanfit)
  u_names <- unique(names(basis$start))
  model_matrix <- full_model.matrix_fast(
    basis$formula,
    newdata, u_names
  )
  # create separate column for each level of categorical variables
  for (i in seq_along(resp_all)) {
    resp <- resp_all[i]
    if (is_categorical(basis$formula[[i]]$family)) {
      resp_levels <- object$levels[[resp]]
      # TODO: glued names to formula?
      newdata[, c(glue::glue("{resp}_{resp_levels}"))] <- NA
      newdata[[resp]] <- NULL
    }
  }
  newdata <- data.frame(newdata, draw = rep(1:n_draws, each = nrow(newdata)))
  n <- nrow(newdata)
  for (i in (fixed + 1):n_time) {
    idx_i <- seq(i, n, by = n_time) # TODO fix for a case of unequal number of time points per group
    idx_i2 <- seq(i, nrow(model_matrix), by = n_time)
    for (j in seq_along(resp_all)) {
      resp <- resp_all[j]
      s <- matrix(samples[[paste0("beta_", resp)]][1:n_draws, i - fixed, ],
                  nrow = n_draws)

      if (is_categorical(basis$formula[[j]]$family)) {
        idx_resp <- which(names(newdata) %in%
                            c(glue::glue("{resp}_{resp_levels}")))
        newdata[idx_i, idx_resp] <- do.call(
          paste0("fitted_", basis$formula[[j]]$family),
          list(
            model_matrix = model_matrix[idx_i2, basis$J[[j]], drop = FALSE],
            samples = s
          )
        )
      } else {
        newdata[idx_i, resp] <- do.call(
          paste0("fitted_", basis$formula[[j]]$family),
          list(
            model_matrix = model_matrix[idx_i2, basis$J[[j]], drop = FALSE],
            samples = s
          )
        )
      }
    }
  }
  newdata
}

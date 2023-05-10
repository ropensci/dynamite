

# Wide format predictor matrix
parse_predictors_wide <- function(dformula, vars, time, group_var) {
  n_vars <- length(vars)
  n_time <- length(time)
  wide_vars <- c(t(outer(vars, time, FUN = "paste", sep = "__")))
  out <- matrix(
    0L,
    nrow = 1L + n_vars * n_time,
    ncol = 1L + n_vars * n_time,
    dimnames = list(c(group_var, wide_vars), c(group_var, wide_vars))
  )
  resp <- get_responses(dformula)
  idx <- 1L + seq_along(time) - n_time
  for (v in vars) {
    idx <- idx + n_time
    i <- which(resp == v)
    if (length(i) > 0L) {
      lag_map <- extract_lags(get_lag_terms(dformula[i]))
      lag_map <- lag_map[lag_map$var %in% vars, ]
      for (j in seq_len(nrow(lag_map))) {
        idx_lag <- 1L + seq_len(n_time) + n_time * (which(vars == lag_map$var[j]) - 1L)
        idx_arr <- cbind(idx, idx_lag - lag_map$k[j])
        idx_arr <- idx_arr[idx_arr[, 2L] >= min(idx_lag), , drop = FALSE]
        out[idx_arr] <- 1L
      }
      nonlags <- get_nonlag_terms(dformula[i])[[1L]]
      for (j in seq_along(nonlags)) {
        idx_nonlag <- 1L +
          seq_len(n_time) + n_time * (which(vars == nonlags[j]) - 1L)
        out[cbind(idx, idx_nonlag)] <- 1L
      }
    }
  }
  out + t(out)
}

#' Long-format Predictor Matrix for Imputation
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula.
#' @param vars \[`character()`]\cr Names of the variables with NA values in the
#'   data.
#' @param all_vars \[`character()`]\cr Names of all data variables.
#' @noRd
parse_predictors <- function(dformula, vars, all_vars) {
  n_vars <- length(vars)
  n_all <- length(all_vars)
  out <- matrix(
    0L,
    nrow = n_all,
    ncol = n_all,
    dimnames = list(all_vars, all_vars)
  )
  resp <- get_responses(dformula)
  for (v in vars) {
    i <- which(resp == v)
    if (length(i) > 0L) {
      lag_map <- rbind(
        extract_lags(get_lag_terms(dformula[i])),
        expand.grid(
          src = "",
          var = resp,
          k = attr(dformula, "lags")$k,
          present = FALSE
        )
      )
      for (j in seq_len(nrow(lag_map))) {
        lag <- paste0(lag_map$var[j], "_lag", lag_map$k[j])
        lead <- paste0(v, "_lead", lag_map$k[j])
        out[v, lag] <- 1L
        out[lag_map$var[j], lead] <- 1L #*
        #(lag_map$var[j] %in% vars)
        # out[lead, v] <- 1L
        # out[lead, lag] <- 1L
        # out[lag, v] <- 1L
        # out[lag, lead] <- 1L
      }
      nonlags <- get_nonlag_terms(dformula[i])[[1L]]
      for (j in seq_along(nonlags)) {
        out[v, nonlags[j]] <- 1L
        out[nonlags[j], v] <- 1L * (nonlags[j] %in% vars)
      }
    }
  }
  out
}

#' Estimate a Bayesian Dynamic Multivariate Panel Model With Multiple Imputation
#'
#' Applies multiple imputation using [mice::mice()] to the supplied `data`
#' and fits a dynamic multivariate panel model to each imputed data set using
#' [dynamite::dynamite()]. Posterior samples from each imputation run are
#' combined. The long format `data` is automatically converted to a
#' wide format before imputation to preserve the longitudinal structure, and
#' then converted back to long format for estimation.
#'
#' @family fitting
#' @inheritParams dynamite
#' @param mice_args \[`list()`]\cr
#'   Arguments passed to [mice::mice()] excluding `data`.
#' @export
dynamice <- function(dformula, data, time, group = NULL,
                     priors = NULL, backend = "rstan",
                     verbose = TRUE, verbose_stan = FALSE,
                     stanc_options = list("O1"),
                     threads_per_chain = 1L, grainsize = NULL,
                     debug = NULL, mice_args = list(), impute_method = "default", ...) {
  stopifnot_(
    requireNamespace("mice"),
    "Please install the {.pkg mice} package to use multiple imputation."
  )
  dynamite_check(
    dformula,
    data,
    time,
    group,
    priors,
    verbose,
    verbose_stan,
    stanc_options,
    threads_per_chain,
    grainsize,
    debug
  )
  backend <- try(match.arg(backend, c("rstan", "cmdstanr")), silent = TRUE)
  stopifnot_(
    !inherits(backend, "try-error"),
    "Argument {.arg backend} must be {.val rstan} or {.val cmdstanr}."
  )
  data <- droplevels(data)
  data <- data.table::as.data.table(data)
  stopifnot_(
    any(is.na(data)),
    "Argument {.arg data} does not contain missing values."
  )
  if (is.null(group)) {
    group <- ".group"
    data_names <- names(data)
    while (group %in% data_names) {
      group <- paste0(group, "_")
    }
    data[[group]] <- 1L
  }
  d <- match.call()$data
  data_name <- ifelse_(
    is.symbol(d),
    deparse1(d),
    "NULL"
  )
  value_vars <- setdiff(names(data), c(time, group))

  # Wide format imputation
  if (impute_method == "wide") {
    data_wide <- data.table::dcast(
      data = data,
      formula = as.formula(paste0(group, " ~ ", time)),
      value.var = value_vars,
      sep = "__"
    )
    if (length(value_vars) == 1L) {
      names(data_wide)[-1L] <- paste0(value_vars, "__", names(data_wide)[-1L])
    }
    wide_vars <- names(data_wide)[-1L]
    mice_args$data <- data_wide
    # mice does this automatically
    # complete_vars <- wide_vars[
    #   vapply(
    #     data_wide[, .SD, .SDcols = wide_vars],
    #     function(x) all(!is.na(x)),
    #     logical(1L)
    #   )
    # ]
    pred_mat <- parse_predictors_wide(
      dformula = dformula,
      vars = value_vars,
      time = sort(unique(data[[time]])),
      group_var = group
    )
    #pred_mat[complete_vars, ] <- 0L
    mice_args$predictorMatrix <- pred_mat
    nn <- colnames(pred_mat)
    n_time <- length(unique(data[[time]]))
    # mice_args$visitSequence <-
    #   c(nn[1], nn[1+2*n + 1:n],
    #     c(matrix(nn[1+1:(2*n)], 2, byrow=TRUE)))
    mice_args$method <- method#"norm"
    imputed <- do.call(mice::mice, args = mice_args)
  } else {
    # Long format imputation

    max_lag <- max(
      extract_lags(get_lag_terms(dformula))$k,
      onlyif(
        !is.null(attr(dformula, "lags")),
        attr(dformula, "lags")$k
      )
    )
    # Ensure that lags/leads exist for imputation by adding/incrementing lags()
    dformula_tmp <- dformula
    if (is.null(attr(dformula_tmp, "lags"))) {
      dformula_tmp <- dformula_tmp + lags(k = max_lag)
    } else {
      attr(dformula_tmp, "lags")$k <- max_lag
    }
    tmp <- dynamite(
      dformula = dformula_tmp,
      data = data,
      time = time,
      group = group,
      priors = get_priors(dformula, data, time, group),
      backend = backend,
      verbose = FALSE,
      debug = list(no_compile = TRUE, no_sampling = TRUE)
    )
    data_forward <- tmp$data
    data_rev <- data.table::copy(data)
    data.table::set(
      x = data_rev,
      j = time,
      value = rev(data[[time]])
    )
    tmp <- dynamite(
      dformula = dformula_tmp,
      data = data_rev,
      time = time,
      group = group,
      priors = get_priors(dformula, data, time, group),
      backend = backend,
      verbose = FALSE,
      debug = list(no_compile = TRUE, no_sampling = TRUE)
    )
    data_backward <- tmp$data
    data.table::set(
      x = data_backward,
      j = time,
      value = rev(data_backward[[time]])
    )
    data.table::setkeyv(data_backward, c(group, time))
    lag_stoch <- get_responses(tmp$dformulas$lag_stoch)
    lag_pred <- get_responses(tmp$dformulas$lag_pred)
    lag_det <- get_responses(tmp$dformulas$lag_det)
    lead_stoch <- gsub("_lag", "_lead", lag_stoch, fixed = TRUE)
    lead_pred <- gsub("_lag", "_lead", lag_pred, fixed = TRUE)
    lead_det <- gsub("_lag", "_lead", lag_det, fixed = TRUE)
    lags <- c(lag_stoch, lag_pred, lag_det)
    leads <- c(lead_stoch, lead_pred, lead_det)
    data_backward <- data_backward[, .SD, .SDcols = lags]
    colnames(data_backward) <- leads
    mice_args$data <- cbind_datatable(data_forward, data_backward)

    # this is done automatically by mice?
    # complete_vars <- vapply(
    #   mice_args$data,
    #   function(x) all(!is.na(x)),
    #   logical(1L)
    # )
    pred_mat <- parse_predictors(
      dformula = dformula_tmp,
      vars = value_vars,
      all_vars = colnames(mice_args$data)
    )
    # for comparison purposes do not impute with group or time
    # but these are needed for lag and lead "imputation"
    #pred_mat[, group] <- 1L
    #pred_mat[, time] <- 1L
    # pred_mat[complete_vars, ] <- 0L
    pred_mat[c(group, time), ] <- 0L
    pred_mat["x_lag1", "x"] <- 1L
    pred_mat["y_lag1", "y"] <- 1L
    pred_mat["x_lead1", "x"] <- 1L
    pred_mat["y_lead1", "y"] <- 1L
    pred_mat["x_lag1", group] <- 1L
    pred_mat["y_lag1", group] <- 1L
    pred_mat["x_lead1", group] <- 1L
    pred_mat["y_lead1", group] <- 1L
    pred_mat["x_lag1", time] <- 1L
    pred_mat["y_lag1", time] <- 1L
    pred_mat["x_lead1", time] <- 1L
    pred_mat["y_lead1", time] <- 1L

    method <- rep("", length = ncol(pred_mat))
    names(method) <- colnames(pred_mat)
    method[value_vars] <- "norm"
    # default is to impute only responses using non-imputed lead and lag
    if (impute_method == "impute_and_use_all") {
      method[lag_stoch] <- "lag"
      method[lead_stoch] <- "lead"
    }
    if (impute_method == "impute_and_use_only_lags") {
      # impute and use only lags
      method[lag_stoch] <- "lag"
      pred_mat[, "y_lead1"] <- 0L
      pred_mat[, "x_lead1"] <- 0L
    }
    if (impute_method == "use_only_lags") {
      # impute only responses, without using the leads
      pred_mat[, "y_lead1"] <- 0L
      pred_mat[, "x_lead1"] <- 0L
    }
    # if (impute_method == "completecases") {
    #   method[lag_stoch] <- "lag"
    #   method[lead_stoch] <- "lead"
    #   ignore <- rep(FALSE, length = nrow(data))
    #   ignore[rowSums(is.na(tmp$data)) > 0] <- TRUE
    # } else {
      ignore <- rep(FALSE, length = nrow(data))
      ignore[which(data[[time]] <= max_lag)] <- TRUE
    # }
    mice_args$predictorMatrix <- pred_mat
    #method[lag_stoch] <- paste0("~ I(lag_(", get_rhs(tmp$dformulas$lag_stoch), "))")
    #method[lag_pred] <- paste0("~ I(lag_(", get_rhs(tmp$dformulas$lag_pred), "))")
    #method[lag_det] <- paste0("~ I(lag_(", get_rhs(tmp$dformulas$lag_det), "))")
    #method[lead_stoch] <- paste0("~ I(lead_(", get_rhs(tmp$dformulas$lag_stoch), "))")
    #method[lead_pred] <- paste0("~ I(lead_(", get_rhs(tmp$dformulas$lag_pred), "))")
    #method[lead_det] <- paste0("~ I(lead_(", get_rhs(tmp$dformulas$lag_det), "))")
    mice_args$method <- method

    ## use only complete case rows for imputation model?
    # ignore <- rep(FALSE, length = nrow(mice_args$data))
    # ignore[rowSums(is.na(mice_args$data)) > 0] <- TRUE
    mice_args$ignore <- ignore

    imputed <- do.call(mice::mice, args = mice_args)
  }
  m <- imputed$m
  e <- new.env()
  sf <- vector(mode = "list", length = m)
  filenames <- character(m)
  model <- NULL
  tmp <- NULL
  outputdir <- paste0("C:/MyTemp/repos/modelling/imputation/", impute_method)
  do.call(file.remove, list(list.files(outputdir, full.names = TRUE)))
  for (i in seq_len(m)) {

    if (impute_method == "wide") {
      # Wide to long in wide format imputation
      melt_data <- data.table::as.data.table(mice::complete(imputed, action = i))
      # Need to construct melt call dynamically because of patterns
      melt_call_str <- ifelse_(
        length(value_vars) == 1L,
        paste0(
          "data.table::melt(",
          "melt_data, ",
          "id.vars = group, ",
          "variable.name = time, ",
          "value.name = value_vars",
          ")"
        ),
        paste0(
          "data.table::melt(",
          "melt_data, ",
          "id.vars = group, ",
          "variable.name = time, ",
          "measure.vars = patterns(",
          paste0(value_vars, " = '^", value_vars, "__'", collapse = ", "),
          "))"
        )
      )
      data_imputed <- eval(str2lang(melt_call_str))
      if (length(value_vars) == 1L) {
        pattern <- paste0("^", value_vars, "__")
        data.table::set(
          x = data_imputed,
          j = time,
          value = as.numeric(gsub(pattern, "", data_imputed[[time]]))
        )
      }
    } else {
      data_imputed <- mice::complete(imputed, action = i)
      data_imputed <- data_imputed[, c(group, time, value_vars)]
    }
    tmp <- dynamite(
      dformula = dformula,
      #data = data_long,
      data = data_imputed,
      time = time,
      group = group,
      priors = get_priors(dformula, data, time, group),
      backend = backend,
      verbose = FALSE,
      debug = list(
        dots = TRUE,
        model = TRUE,
        model_code = TRUE,
        no_compile = i > 1L,
        no_sampling = TRUE,
        stan_input = TRUE
      ),
      ...
    )
    onlyif(i == 1L, e$model <- tmp$model)
    dots <- tmp$dots
    if (identical(backend, "rstan")) {
      e$args <- c(
        list(object = e$model, data = tmp$stan_input$sampling_vars),
        dots
      )
      sf[[i]] <- with(e, {
        do.call(rstan::sampling, args)
      })
    } else {
      dots$chain_ids <- i
      e$args <- c(
        list(data = tmp$stan_input$sampling_vars),
        dots,
        output_dir = outputdir,
        threads_per_chain = onlyif(threads_per_chain > 1L, threads_per_chain)
      )
      sampling_out <- with(e, {
        do.call(model$sample, args)
      })
      #stanfit <- rstan::read_stan_csv(sampling_out$output_files())
      #tmp$stanfit <- stanfit
      filenames[i] <- sampling_out$output_files()
    }
  }
  if (identical(backend, "rstan")) {
    stanfit <- rstan::sflist2stanfit(sf)
  } else {
    stanfit <- rstan::read_stan_csv(filenames)
    stanfit@stanmodel <- methods::new("stanmodel", model_code = tmp$model_code)
  }
  # TODO return object? How is this going to work with update?
  structure(
    list(
      stanfit = stanfit,
      dformulas = tmp$dformulas,
      data = data,
      data_name = data_name,
      stan = tmp$stan,
      group_var = group,
      time_var = time,
      priors = priors,
      backend = backend,
      call = tmp$call, # TODO?
      imputed = imputed
    ),
    class = "dynamitefit"
  )
}

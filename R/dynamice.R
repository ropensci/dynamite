#' Estimate a Bayesian Dynamic Multivariate Panel Model With Multiple Imputation
#'
#' Applies multiple imputation using [mice::mice()] to the supplied `data`
#' and fits a dynamic multivariate panel model to each imputed data set using
#' [dynamite()]. Posterior samples from each imputation run are
#' combined. When using wide format imputation, the long format `data` is
#' automatically converted to a wide format before imputation to preserve the
#' longitudinal structure, and then converted back to long format for
#' estimation.
#'
#' @family fitting
#' @inheritParams dynamite
#' @param mice_args \[`list()`]\cr
#'   Arguments passed to [mice::mice()] excluding `data`.
#' @param impute_format \[`character(1)`]\cr Format of the data that will be
#'   passed to the imputation method. Should be either `"wide"` (the default)
#'   or `"long"` corresponding to wide format and long format imputation.
#' @param keep_imputed \[`logical(1)`]\cr Should the imputed datasets be
#'   kept in the return object? The default is `FALSE`. If `TRUE`, the
#'   imputations will be included in the `imputed` field in the return object
#'   that is otherwise `NULL`.
#' @param stan_csv_dir \[`character(1)`]\cr A directory path to output the
#'   Stan .csv files when `backend` is `"cmdstanr"`. The files are saved here
#'   via `$save_output_files()` to avoid garbage collection between sampling
#'   runs with different imputed datasets.
#' @export
#'
dynamice <- function(dformula, data, time, group = NULL,
                     priors = NULL, backend = "rstan",
                     verbose = TRUE, verbose_stan = FALSE,
                     stanc_options = list("O0"),
                     threads_per_chain = 1L, grainsize = NULL,
                     custom_stan_model = NULL, debug = NULL,
                     mice_args = list(), impute_format = "wide",
                     keep_imputed = FALSE, stan_csv_dir = tempdir(), ...) {
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
    custom_stan_model,
    debug
  )
  backend <- try(match.arg(backend, c("rstan", "cmdstanr")), silent = TRUE)
  stopifnot_(
    !inherits(backend, "try-error"),
    "Argument {.arg backend} must be either {.val rstan} or {.val cmdstanr}."
  )
  impute_format <- try(
    match.arg(impute_format, c("long", "wide")),
    silent = TRUE
  )
  stopifnot_(
    !inherits(impute_format, "try-error"),
    "Argument {.arg impute_format} must be either {.val long} or {.val wide}."
  )
  stopifnot_(
    checkmate::test_flag(x = keep_imputed),
    "Argument {.arg keep_imputed} must be a single {.cls logical} value."
  )
  stopifnot_(
    any(is.na(data)),
    "Argument {.arg data} does not contain missing values."
  )
  data <- droplevels(data)
  data <- data.table::as.data.table(data)
  if (is.null(group)) {
    group <- ".group"
    data_names <- names(data)
    while (group %in% data_names) {
      group <- paste0(group, "_")
    }
    data[, (group) := 1L]
  }
  d <- match.call()$data
  data_name <- ifelse_(
    is.symbol(d),
    deparse1(d),
    "NULL"
  )
  value_vars <- setdiff(names(data), c(time, group))
  imputed <- ifelse_(
    identical(impute_format, "long"),
    impute_long(dformula, data, time, group, backend, mice_args),
    impute_wide(dformula, data, time, group, mice_args)
  )
  m <- imputed$m
  e <- new.env()
  sf <- vector(mode = "list", length = m)
  filenames <- character(m)
  tmp <- NULL
  for (i in seq_len(m)) {
    data_imputed <- ifelse_(
      identical(impute_format, "long"),
      get_imputed_long(imputed, time, group, value_vars, i),
      get_imputed_wide(imputed, time, group, value_vars, i)
    )
    tmp <- dynamite(
      dformula = dformula,
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
        threads_per_chain = onlyif(threads_per_chain > 1L, threads_per_chain)
      )
      sampling_out <- with(e, {
        do.call(model$sample, args)
      })
      sampling_out$save_output_files(dir = stan_csv_dir)
      filenames[i] <- sampling_out$output_files()
    }
  }
  if (identical(backend, "rstan")) {
    stanfit <- rstan::sflist2stanfit(sf)
  } else {
    stanfit <- cmdstanr::as_cmdstan_fit(filenames, check_diagnostics = FALSE)
  }
  n_draws <- ifelse_(is.null(stanfit), 0L, get_ndraws(stanfit))
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
      permutation = sample(n_draws),
      imputed = onlyif(keep_imputed, imputed),
      call = tmp$call # TODO?
    ),
    class = "dynamitefit"
  )
}

#' Get i:th Imputed Data Set When Using Long Format Imputation
#'
#' @param imputed \[`mids`]\cr Output of [mice::mice()].
#' @param time \[`character(1)`]\cr Name of the time variable.
#' @param group \[`character(1)`]\cr Name of the group variable.
#' @param value_vars  \[`character()`]
#'   Names of data variables other than `time` and `group`.
#' @param i  \[`integer(1)`] Index of the imputed data set.
#' @noRd
get_imputed_long <- function(imputed, time, group, value_vars, i) {
  data_imputed <- mice::complete(imputed, action = i)
  data_imputed[, c(group, time, value_vars)]
}

#' Get i:th Imputed Data Set When Using Wide Format Imputation
#'
#' @param imputed \[`mids`]\cr Output of [mice::mice()].
#' @param time \[`character(1)`]\cr Name of the time variable.
#' @param group \[`character(1)`]\cr Name of the group variable.
#' @param value_vars  \[`character()`]\cr
#'   Names of data variables other than `time` and `group`.
#' @param i  \[`integer(1)`]\cr Index of the imputed data set.
#' @noRd
get_imputed_wide <- function(imputed, time, group, value_vars, i) {
  melt_data <- data.table::as.data.table(mice::complete(imputed, action = i))
  # Need to construct melt call dynamically because of patterns argument
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
      paste0(value_vars, " = '^", value_vars, "_'", collapse = ", "),
      "))"
    )
  )
  data_imputed <- eval(str2lang(melt_call_str))
  if (length(value_vars) == 1L) {
    pattern <- paste0("^", value_vars, "_")
    data.table::set(
      x = data_imputed,
      j = time,
      value = as.numeric(gsub(pattern, "", data_imputed[[time]]))
    )
  }
  data_imputed
}

#' Generate Imputed Data Sets in Long Format
#'
#' @inheritParams dynamite
#' @noRd
impute_long <- function(dformula, data, time, group, backend, mice_args) {
  value_vars <- setdiff(names(data), c(time, group))
  max_lag <- max(
    extract_lags(get_lag_terms(dformula))$k,
    onlyif(
      !is.null(attr(dformula, "lags")),
      attr(dformula, "lags")$k
    )
  )
  # Ensure that lags/leads exist for imputation by including all lags
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
  rhs_stoch <- get_rhs(tmp$dformulas$lag_stoch)
  rhs_pred <- get_rhs(tmp$dformulas$lag_pred)
  rhs_det <- get_rhs(tmp$dformulas$lag_det)
  lead_stoch <- gsub("_lag", "_lead", lag_stoch, fixed = TRUE)
  lead_pred <- gsub("_lag", "_lead", lag_pred, fixed = TRUE)
  lead_det <- gsub("_lag", "_lead", lag_det, fixed = TRUE)
  lags <- c(lag_stoch, lag_pred, lag_det)
  leads <- c(lead_stoch, lead_pred, lead_det)
  rhs <- c(rhs_stoch, rhs_pred, rhs_det)
  data_backward <- data_backward[, .SD, .SDcols = lags]
  colnames(data_backward) <- leads
  mice_args$data <- cbind_datatable(data_forward, data_backward)
  pred_mat <- parse_predictors_long(
    dformula = dformula_tmp,
    time_var = time,
    group_var = group,
    all_vars = colnames(mice_args$data)
  )
  pred_mat[lag_stoch, time] <- 1L
  pred_mat[lag_pred, time] <- 1L
  pred_mat[lag_det, time] <- 1L
  pred_mat[lead_stoch, time] <- 1L
  pred_mat[lead_pred, time] <- 1L
  pred_mat[lead_det, time] <- 1L
  pred_mat[lag_stoch, group] <- 1L
  pred_mat[lag_pred, group] <- 1L
  pred_mat[lag_det, group] <- 1L
  pred_mat[lead_stoch, group] <- 1L
  pred_mat[lead_pred, group] <- 1L
  pred_mat[lead_det, group] <- 1L
  pred_mat <- pred_mat[names(mice_args$data), names(mice_args$data)]
  if (n_unique(mice_args$data[[group]]) == 1L) {
    pred_mat[, group] <- 0L
  }
  mice_args$predictorMatrix <- pred_mat
  method <- rep("", length = ncol(pred_mat))
  names(method) <- colnames(pred_mat)
  method[value_vars] <- "norm"
  method[lag_stoch] <- "lag"
  method[lag_pred] <- "lag"
  method[lag_det] <- "lag"
  method[lead_stoch] <- "lead"
  method[lead_pred] <- "lead"
  method[lead_det] <- "lead"
  mice_args$method <- method
  blots_fun <- function(y) {
    list(
      # for some reason mice drops the grouping variable sometimes
      # we carry it via the blots just in case to compute the lags/leads
      group_val = mice_args$data[[group]],
      group_var = group,
      resp = y
    )
  }
  mice_args$blots <- c(
    mice_args$blots,
    stats::setNames(lapply(rhs, blots_fun), lags),
    stats::setNames(lapply(rhs, blots_fun), leads)
  )
  do.call(mice::mice, args = mice_args)
}

#' Generate Imputed Data Sets in Wide Format
#'
#' @inheritParams dynamite
#' @noRd
impute_wide <- function(dformula, data, time, group, mice_args) {
  value_vars <- setdiff(names(data), c(time, group))
  data_wide <- data.table::dcast(
    data = data,
    formula = as.formula(paste0(group, " ~ ", time)),
    value.var = value_vars,
    sep = "_"
  )
  if (length(value_vars) == 1L) {
    names(data_wide)[-1L] <- paste0(value_vars, "_", names(data_wide)[-1L])
  }
  mice_args$data <- data_wide
  n_time <- n_unique(data[[time]])
  pred_mat <- parse_predictors_wide(
    dformula = dformula,
    value_vars = value_vars,
    idx_time = seq_len(n_time),
    group_var = group
  )
  mice_args$predictorMatrix <- pred_mat
  do.call(mice::mice, args = mice_args)
}

#' Long-format Predictor Matrix for Imputation
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula.
#' @param time_var \[`character(1)`]\cr Name of the time variable.
#' @param group_var \[`character(1)`]\cr Name of the grouping variable.
#' @param all_vars \[`character()`]\cr Names of all data variables.
#' @noRd
parse_predictors_long <- function(dformula, time_var, group_var, all_vars) {
  value_vars <- setdiff(all_vars, c(time_var, group_var))
  pred_vars <- c(value_vars, time_var, group_var)
  n_vars <- length(value_vars)
  out <- matrix(
    0L,
    n_vars + 2L,
    n_vars + 2L,
    dimnames = list(pred_vars, pred_vars)
  )
  g <- get_dag(dformula, project = TRUE, covariates = TRUE, format = "lag")
  out[seq_len(n_vars), seq_len(n_vars)] <- g$A
  mb <- lapply(value_vars, function(y) get_markov_blanket(g, y))
  names(mb) <- value_vars
  for (y in value_vars) {
    out[y, mb[[y]]] <- 1L
  }
  out
}

#' Wide format Predictor Matrix for Imputation
#'
#' @param dformula \[`dynamiteformula`]\cr The model formula.
#' @param value_vars \[`character()`]\cr
#'   Names of data variables other than `time` and `group`.
#' @param idx_time \[`integer()`]\cr Time point indices.
#' @param group_var \[`character(1)`]\cr Name of the grouping variable.
#' @noRd
parse_predictors_wide <- function(dformula, value_vars, idx_time, group_var) {
  n_vars <- length(value_vars)
  n_time <- length(idx_time)
  wide_vars <- c(t(outer(value_vars, idx_time, FUN = "paste", sep = "_")))
  pred_vars <- c(group_var, wide_vars)
  out <- matrix(
    0L,
    nrow = 1L + n_vars * n_time,
    ncol = 1L + n_vars * n_time,
    dimnames = list(pred_vars, pred_vars)
  )
  resp <- get_responses(dformula)
  g <- get_dag(dformula, project = TRUE, covariates = TRUE, format = "default")
  mb <- lapply(resp, function(y) get_markov_blanket(g, paste0(y, "_{t}")))
  e <- new.env()
  for (t in idx_time) {
    e$t <- t
    for (i in seq_along(resp)) {
      mb_t <- vapply(
        mb[[i]],
        function(y) glue::glue(y, .envir = e),
        character(1L)
      )
      mb_t <- mb_t[mb_t %in% wide_vars]
      resp_ti <- paste0(resp[i], "_", t)
      out[resp_ti, mb_t] <- 1L
    }
  }
  out
}

#' Compute Lagging Values of an Imputed Response
#'
#' Function for computing the lagged values of an imputed response in `mice`.
#'
#' @inheritParams mice::mice.impute.norm
#' @param group_val \[`vector()`]\cr Values of the grouping variable.
#' @param group_var \[`character(1)`]\cr Name of the grouping variable.
#' @param resp \[`character(1)`]\cr Name of the response variable.
#' @keywords internal
#' @export
mice.impute.lag <- function(y, ry, x, wy = NULL, group_val,
                            group_var, resp, ...) {
  if (is.null(wy)) {
    wy <- !ry
  }
  if (!group_var %in% colnames(x)) {
    x_names <- colnames(x)
    x <- cbind(x, group_val)
    colnames(x) <- c(x_names, group_var)
  }
  imputed <- data.table::as.data.table(x)[,
    list(lag = lag_(resp)),
    by = group_var,
    env = list(resp = resp, lag_ = "lag_")
  ]$lag
  imputed[wy]
}

#' Compute Leading Values of an Imputed Response
#'
#' Function for computing the leading values of an imputed response in `mice`.
#'
#' @inheritParams mice::mice.impute.norm
#' @param group_val \[`vector()`]\cr Values of the grouping variable.
#' @param group_var \[`character(1)`]\cr Name of the grouping variable.
#' @param resp \[`character(1)`]\cr Name of the response variable.
#' @keywords internal
#' @export
mice.impute.lead <- function(y, ry, x, wy = NULL, group_val,
                             group_var, resp, ...) {
  if (is.null(wy)) {
    wy <- !ry
  }
  if (!group_var %in% colnames(x)) {
    x_names <- colnames(x)
    x <- cbind(x, group_val)
    colnames(x) <- c(x_names, group_var)
  }
  imputed <- data.table::as.data.table(x)[,
    list(lead = lead_(resp)),
    by = group_var,
    env = list(resp = resp, lead_ = "lead_")
  ]$lead
  imputed[wy]
}

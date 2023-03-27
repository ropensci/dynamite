#' Estimate a Bayesian Dynamic Multivariate Panel Model
#'
#' Fit a Bayesian dynamic multivariate panel model (DMPM) using Stan for
#' Bayesian inference. The \pkg{dynamite} package supports a wide range of
#' distributions and allows the user to flexibly customize the priors for the
#' model parameters. The dynamite model is specified using standard \R formula
#' syntax via [dynamite::dynamiteformula()]. For more information and examples,
#' see 'Details' and the package vignette.
#'
#' Any univariate unbounded continuous distributions supported by Stan can be
#' used as a prior for model parameters (the distribution is automatically
#' truncated to positive side for constrained parameters). In addition, any
#' univariate distribution bounded to the positive real line can be used as a
#' prior for parameters constrained to be positive.
#' See Stan function reference at
#' \url{https://mc-stan.org/users/documentation/} for details. For custom
#' priors, you should first get the default priors with [dynamite::get_priors()]
#' function, and then modify the `priors` column of the obtained data frame
#' before supplying it to the `dynamite` function.
#'
#' The default priors for regression coefficients are based on the standard
#' deviation of the covariates at the first non-fixed time point. In case this
#' is 0 or NA, it is transformed to (arbitrary) 0.5. The final prior is then
#' normal distribution with zero mean and two times this standard deviation.
#'
#' The prior for the correlation structure of the random intercepts is defined
#' via the Cholesky decomposition of the correlation matrix, as
#' `lkj_corr_cholesky(1)`. See
#' \url{https://mc-stan.org/docs/functions-reference/cholesky-lkj-correlation-distribution.html}
#' for details.
#'
#' The best-case scalability of `dynamite` in terms of data size should be
#' approximately linear in terms of number of time points and and number of
#' groups, but as wall-clock time of the MCMC algorithms provided by Stan can
#' depend on the discrepancy of the data and the model (and the subsequent
#' shape of the posterior), this can vary greatly.
#'
#' @export
#' @family fitting
#' @rdname dynamite
#' @param dformula \[`dynamiteformula`]\cr The model formula.
#'   See [dynamite::dynamiteformula()] and 'Details'.
#' @param data
#'   \[`data.frame`, `tibble::tibble`, or `data.table::data.table`]\cr
#'   The data that contains the variables in the
#'   model. Supported column types are `integer`, `logical`, `double`, and
#'   `factor`. Columns of type `character` will be converted to factors.
#'   Unused factor levels will be dropped. The `data` can contain missing
#'   values which will simply be ignored in the estimation in a case-wise
#'   fashion (per time-point and per channel). Input `data` is converted to
#'   channel specific matrix representations via [stats::model.matrix.lm()].
#' @param time \[`character(1)`]\cr A column name of `data` that denotes the
#'   time index of observations. If this variable is a factor, the integer
#'   representation of its levels are used internally for defining the time
#'   indexing.
#' @param group \[`character(1)`]\cr A column name of `data` that denotes the
#'   unique groups or `NULL` corresponding to a scenario without any groups.
#'   If `group` is `NULL`, a new column `.group` is created with constant
#'   value `1L` is created indicating that all observations belong to the same
#'   group. In case of name conflicts with `data`, see the `group_var` element
#'   of the return object to get the column name of the new variable.
#' @param priors \[`data.frame`]\cr An optional data frame with prior
#'   definitions. See [dynamite::get_priors()] and 'Details'.
#' @param backend \[`character(1)`]\cr Defines the backend interface to Stan,
#'   should be  either `"rstan"` (the default) or `"cmdstanr"`. Note that
#'   `cmdstanr` needs to be installed separately as it is not on CRAN. It also
#'   needs the actual `CmdStan` software. See https://mc-stan.org/cmdstanr/ for
#'   details.
#' @param verbose \[`logical(1)`]\cr All warnings and messages are suppressed
#'   if set to `FALSE`. Defaults to `TRUE`.
#' @param verbose_stan \[`logical(1)`]\cr This is the `verbose` argument for
#'   [rstan::sampling()]. Defaults to `FALSE`.
#' @param stanc_options \[`list()`]\cr This is the `stanc_options` argument
#'   passed to the compile method of a `CmdStanModel` object via
#'   [cmdstanr::cmdstan_model()] when `backend = "cmdstanr"`.
#'   Defaults to `list("O1")` to enable level one compiler optimizations.
#' @param debug \[`list()`]\cr A named list of form `name = TRUE` indicating
#'   additional objects in the environment of the `dynamite` function which are
#'   added to the return object. Additionally, values `no_compile = TRUE` and
#'   `no_sampling = TRUE` can be used to skip the compilation of the Stan code
#'   and sampling steps respectively. This can be useful for debugging when
#'   combined with `model_code = TRUE`, which adds the Stan model code to the
#'   return object.
#' @param ... For `dynamite()`, additional arguments to [rstan::sampling()] or
#'  [cmdstanr::sample()], such as `chains` and `cores` (`parallel_chains` in
#'  `cmdstanr`). For `summary()`, additional arguments to
#'  [dynamite::as.data.frame.dynamitefit()]. For `print()`, further arguments
#'  to the print method for tibbles (see [tibble::formatting]). Not used for
#'  `formula()`.
#' @return `dynamite` returns a `dynamitefit` object which is a list containing
#'   the following components:
#'
#'   * `stanfit`\cr A `stanfit` object, see [rstan::sampling()] for details.
#'   * `dformulas`\cr A list of `dynamiteformula` objects for internal use.
#'   * `data`\cr A processed version of the input `data`.
#'   * `data_name`\cr Name of the input data object.
#'   * `stan`\cr A `list` containing various elements related to Stan model
#'     construction and sampling.
#'   * `group_var`\cr Name of the variable defining the groups.
#'   * `time_var`\cr Name of the variable defining the time index.
#'   * `priors`\cr Data frame containing the used priors.
#'   * `backend`\cr Either `"rstan"` or `"cmdstanr"` indicating which
#'     package was used in sampling.
#'   * `call`\cr Original function call as an object of class `call`.
#'
#' @srrstats {G2.9} Potential loss of information is reported by `dynamite`.
#' @srrstats {RE1.1} Documented in `dformula` parameter.
#' @srrstats {RE1.4} Documented in `data` parameter.
#' @srrstats {RE2.0} Transformations are documented (in parameter `data`) and
#'   warnings are issued when appropriate.
#' @srrstats {RE2.1} Non-finite values are not allowed, NAs are appropriately
#'   considered.
#' @srrstats {RE3.0} Convergence diagnostics are delegated to Stan.
#' @srrstats {RE3.1} Stan messages can be suppressed.
#' @srrstats {RE4.0} `dynamite` returns a `dynamitefit` object.
#' @srrstats {RE4.1} `dynamitefit` object can be generated without sampling.
#' @srrstats {RE4.4} The model is specified via formula objects.
#' @srrstats {RE4.8} `dynamitefit` object contains the response variables
#' @srrstats {RE4.13} `dynamitefit` object contains the predictor variables
#' @srrstats {BS1.1} Data input is documented in the `data` parameter.,
#' @srrstats {BS1.3, BS1.3b} Computational parameters are delegated to Stan.
#' @srrstats {BS2.6} Checks for computational parameters are performed by Stan.
#' @srrstats {BS2.7} Starting values can be controlled via `...`.
#' @srrstats {BS2.9} Chains have different starting values by default.
#' @srrstats {BS2.12} Verbosity of output can be controlled via `...`.
#' @srrstats {BS2.13} Progress indicators can be suppressed without suppressing
#'   other output, via `...`.
#' @srrstats {BS2.14} Warnings can be suppressed by setting `verbose = FALSE`-
#' @srrstats {BS3.0} Documented in `data` parameter.
#' @srrstats {BS4.0} Stan is referenced.
#' @srrstats {BS5.0, BS5.1, BS5.2, BS5.3, BS5.5}
#'   Available from the resulting `dynamitefit` object.
#' @srrstats {RE5.0} The scalability of the algorithms is studied to some
#'   extent in the tests and noted here. As the computational algorithms are
#'   based on Stan, the  scalability of the package depends directly on the
#'   scalability of Stan.
#' @references
#' Santtu Tikka and Jouni Helske (2023). `dynamite`: An \R Package for Dynamic
#' Multivariate Panel Models. arXiv preprint,
#' <https://arxiv.org/abs/2302.01607/>.
#'
#' Jouni Helske and Santtu Tikka (2022). Estimating Causal Effects
#' from Panel Data with Dynamic Multivariate Panel Models. SocArxiv preprint,
#' <https://osf.io/preprints/socarxiv/mdwu5/>.
#' @examples
#' \donttest{
#' fit <- dynamite(
#'   dformula = obs(y ~ -1 + varying(~x), family = "gaussian") +
#'     lags(type = "varying") +
#'     splines(df = 20),
#'   gaussian_example,
#'   "time",
#'   "id",
#'   chains = 1,
#'   refresh = 0
#' )
#' }
#'
dynamite <- function(dformula, data, time, group = NULL,
                     priors = NULL, backend = "rstan",
                     verbose = TRUE, verbose_stan = FALSE,
                     stanc_options = list("O1"), debug = NULL, ...) {
  stopifnot_(
    !missing(dformula),
    "Argument {.arg dformula} is missing."
  )
  stopifnot_(
    !missing(data),
    "Argument {.arg data} is missing."
  )
  stopifnot_(
    !missing(time),
    "Argument {.var time} is missing."
  )
  stopifnot_(
    is.dynamiteformula(dformula),
    "Argument {.arg dformula} must be a {.cls dynamiteformula} object."
  )
  stopifnot_(
    length(which_stochastic(dformula)) > 0L,
    "Argument {.arg dformula} must contain at least one stochastic channel."
  )
  stopifnot_(
    is.data.frame(data),
    "Argument {.arg data} must be a {.cls data.frame} object."
  )
  stopifnot_(
    checkmate::test_string(
      x = group,
      null.ok = TRUE
    ),
    "Argument {.arg group} must be a single character string or {.code NULL}."
  )
  stopifnot_(
    is.null(group) || !is.null(data[[group]]),
    "Can't find grouping variable {.var {group}} in {.arg data}."
  )
  stopifnot_(
    checkmate::test_string(x = time),
    "Argument {.arg time} must be a single character string."
  )
  stopifnot_(
    !is.null(data[[time]]),
    "Can't find time index variable {.var {time}} in {.arg data}."
  )
  stopifnot_(
    checkmate::test_flag(x = verbose),
    "Argument {.arg verbose} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_flag(x = verbose_stan),
    "Argument {.arg verbose_stan} must be a single {.cls logical} value."
  )
  stopifnot_(
    is.null(debug) || is.list(debug),
    "Argument {.arg debug} must be a {.cls list} or NULL."
  )
  backend <- try(match.arg(backend, c("rstan", "cmdstanr")), silent = TRUE)
  stopifnot_(
    !inherits(backend, "try-error"),
    "Argument {.arg backend} must be {.val rstan} or {.val cmdstanr}."
  )
  if (is.null(group)) {
    group <- ".group"
    data_names <- names(data)
    while (group %in% data_names) {
      group <- paste0(group, "_")
    }
    data[[group]] <- 1L
  }
  data_name <- attr(dformula, "data_name")
  if (is.null(data_name)) {
    d <- match.call()$data
    data_name <- ifelse_(
      is.symbol(d),
      deparse1(d),
      ""
    )
  }
  data <- parse_data(dformula, data, group, time, verbose)
  dformula <- parse_past(dformula, data, group, time)
  dformulas <- parse_lags(dformula, data, group, time, verbose)
  evaluate_deterministic(dformulas, data, group, time)
  dformulas <- parse_components(dformulas, data, group, time)
  stan_out <- dynamite_stan(
    dformulas,
    data,
    data_name,
    group,
    time,
    priors,
    backend,
    verbose,
    verbose_stan,
    stanc_options,
    debug,
    ...
  )
  # extract elements for debug argument
  stan_input <- stan_out$stan_input
  model_code <- stan_out$model_code
  stanfit <- stan_out$stanfit
  dynamite_call <- match.call()
  if (!is.null(debug) && !is.null(debug$stanfit)) {
    dynamite_call$debug$stanfit <- NULL
    if (identical(length(dynamite_call$debug), 0L)) {
      dynamite_call$debug <- NULL
    }
  }
  out <- structure(
    list(
      stanfit = stanfit,
      dformulas = dformulas,
      data = data,
      data_name = data_name,
      stan = stan_input,
      group_var = group,
      time_var = time,
      priors = rbindlist_(stan_input$priors),
      backend = backend,
      call = dynamite_call
    ),
    class = "dynamitefit"
  )
  # Adds any object in the environment of this function to the return object
  # if its name is included in the debug argument.
  for (opt in setdiff(names(debug), names(out))) {
    got <- try(get(x = opt), silent = TRUE)
    out[[opt]] <- onlyif(!inherits(got, "try-error"), got)
  }
  out
}

#' Prepare Data for Stan and Construct a Stan Model for `dynamite`
#'
#' @inheritParams dynamite
#' @param dformulas \[`list()`]\cr Output of `parse_lags`.
#' @param data_name Name of the `data` object.
#' @noRd
dynamite_stan <- function(dformulas, data, data_name, group, time,
                          priors, backend, verbose, verbose_stan,
                          stanc_options, debug, ...) {
  stan_input <- prepare_stan_input(
    dformulas$stoch,
    data,
    group,
    time,
    priors,
    fixed = attr(dformulas$all, "max_lag"),
    verbose
  )
  model_code <- create_blocks(
    indent = 2L,
    backend = backend,
    cg = attr(dformulas$stoch, "channel_groups"),
    cvars = stan_input$channel_vars,
    cgvars = stan_input$channel_group_vars
  )
  sampling_info(dformulas, verbose, debug, backend)
  # if debug$stanfit exists (from the update method) then don't recompile
  model <- ifelse_(
    is.null(debug) || is.null(debug$stanfit),
    dynamite_model(
      compile = is.null(debug) || !isTRUE(debug$no_compile),
      model_code = model_code,
      backend = backend,
      stanc_options = stanc_options
    ),
    debug$stanfit
  )
  dots <- remove_redundant_parameters(stan_input, backend, verbose_stan, ...)
  dots <- check_stan_args(dots, verbose, backend)
  stanfit <- dynamite_sampling(
    sampling = !isTRUE(debug$no_compile) && !isTRUE(debug$no_sampling),
    backend = backend,
    model_code = model_code,
    model = model,
    sampling_vars = stan_input$sampling_vars,
    dots = dots
  )
  list(
    stan_input = stan_input,
    model_code = model_code,
    stanfit = stanfit
  )
}

#' Compile a Stan Model for `dynamite`
#'
#' @inheritParams dynamite
#' @param compile \[`logical(1)`]\cr Should the model be compiled?
#' @param model_code \[`logical(1)`]\cr The model code as a string.
#' @noRd
dynamite_model <- function(compile, model_code, backend, stanc_options) {
  if (compile) {
    if (backend == "rstan") {
      rstan::stan_model(model_code = model_code)
    } else {
      file <- cmdstanr::write_stan_file(model_code)
      cmdstanr::cmdstan_model(file, stanc_options = stanc_options)
    }
  } else {
    NULL
  }
}

#' Sample from the Stan Model for `dynamite`
#'
#' @param sampling \[`logical(1)`]\cr Should sampling be carried out?
#' @param backend \[`character(1)`]\cr `"rstan"` or `"cmdstanr"`.
#' @param model_code \[`logical(1)`]\cr The model code as a string.
#' @param model \[`stanmodel`]\cr The compiled Stan model.
#' @param sampling_vars \[`list()`]\cr Data for Stan sampling.
#' @param dots \[`list()`]\cr Additional arguments for `rstan` or `cmdstanr`.
#' @noRd
dynamite_sampling <- function(sampling, backend, model_code, model,
                              sampling_vars, dots) {
  out <- NULL
  if (sampling) {
    if (backend == "rstan") {
      out <- do.call(
        rstan::sampling,
        c(list(object = model, data = sampling_vars), dots)
      )
    } else {
      sampling_out <- do.call(
        model$sample,
        c(list(data = sampling_vars), dots)
      )
      out <- rstan::read_stan_csv(sampling_out$output_files())
      out@stanmodel <- methods::new("stanmodel", model_code = model_code)
    }
  }
  out
}

#' Print an Informative Message Related to Stan Model Compilation/Sampling
#'
#' @inheritParams dynamite_stan
#' @noRd
sampling_info <- function(dformulas, verbose, debug, backend) {
  if (!verbose || !is.null(debug) || isTRUE(debug$no_compile)) {
    return()
  }
  if (!stan_supports_categorical_logit_glm(backend) &&
    "categorical" %in% get_family_names(dformulas$all)) {
    warning_(
      c(
        "Efficient GLM variant of the categorical likelihood is not
         available in this version of {.pkg {backend}}.",
        `i` = "For more efficient sampling, please install a newer version
               of {.pkg {backend}}."
      )
    )
  }
}

#' Check Arguments Names Of `...` for Stan Sampling
#'
#' @inheritParams dynamite_stan
#' @param dots The `...` arguments of `dynamite` as a `list`
#' @noRd
check_stan_args <- function(dots, verbose, backend) {
  dots_names <- names(dots)
  args <- ifelse_(
    identical(backend, "rstan"),
    c(
      "pars", "chains", "iter", "warmup", "thin", "seed", "init", "check_data",
      "sample_file", "diagnostic_file", "verbose", "algorithm", "control",
      "include", "cores", "open_progress", "show_messages", "chain_id",
      "init_r", "test_grad", "append_samples", "refresh", "save_warmup",
      "enable_random_init"
    ),
    c(
      "seed", "refresh", "init", "save_latent_dynamics", "output_dir",
      "output_basename", "sig_figs", "chains", "parallel_chains", "chain_ids",
      "threads_per_chain", "opencl_ids", "iter_warmup", "iter_sampling",
      "save_warmup", "thin", "max_treedepth", "adapt_engaged", "adapt_delta",
      "step_size", "metric", "metric_file", "inv_metric", "init_buffer",
      "term_buffer", "window", "fixed_param", "show_messages", "diagnostics",
      "cores", "num_cores", "num_chains", "num_warmup", "num_samples",
      "validate_csv", "save_extra_diagnostics", "max_depth", "stepsize"
    )
  )
  valid_args <- dots_names %in% args
  invalid_args <- dots_names[!valid_args]
  dots <- dots[valid_args]
  if (verbose && any(!valid_args)) {
    warning_(
      "{cli::qty(invalid_args)}
       Argument{?s} {.arg {invalid_args}} passed to {backend} sampling function
       {cli::qty(invalid_args)}{?is/are} not recognized and will be ignored."
    )
  }
  dots
}

#' Count The Number of Random Effects in a `dynamiteformula`
#'
#' @inheritParams dynamite
#' @noRd
count_random_effects <- function(dformula, data) {
  cg <- attr(dformula, "channel_groups")
  n_cg <- n_unique(cg)
  M <- 0L
  for (i in seq_len(n_cg)) {
    cg_idx <- which(cg == i)
    j <- cg_idx[1L]
    family <- dformula[[j]]$family
    if (is_multivariate(family)) {
      if (is_multinomial(family)) {
        random_formula <- get_type_formula(dformula[[j]], type = "random")
        cats <- length(cg_idx) - 1L
        M <- M + ifelse_(
          is.null(random_formula),
          0L,
          cats * ncol(stats::model.matrix.lm(random_formula, data))
        )
      } else {
        for (k in cg_idx) {
          random_formula <- get_type_formula(dformula[[k]], type = "random")
          M <- M + ifelse_(
            is.null(random_formula),
            0L,
            ncol(stats::model.matrix.lm(random_formula, data))
          )
        }
      }
    } else {
      random_formula <- get_type_formula(dformula[[j]], type = "random")
      cats <- ifelse_(
        is_categorical(family),
        length(levels(data[[dformula[[j]]$response]])) - 1L,
        1L
      )
      M <- M + ifelse_(
        is.null(random_formula),
        0L,
        cats * ncol(stats::model.matrix.lm(random_formula, data))
      )
    }
  }
  M
}


#' Remove Redundant Parameters When Using `rstan`
#'
#' @param stan_input Output from `prepare_stan_input`.
#' @inheritParams dynamite
#' @noRd
remove_redundant_parameters <- function(stan_input, backend,
                                        verbose_stan, ...) {
  # don't save redundant parameters by default
  # could also remove omega_raw
  dots <- list(...)
  dots$verbose <- onlyif(backend == "rstan", verbose_stan)
  if (identical(dots$algorithm, "Fixed_param")) {
    return(dots)
  }
  M <- stan_input$sampling_vars$M
  if (is.null(dots$pars) && M > 0 && backend == "rstan") {
    dots$pars <- c("nu_raw", "nu", "L")
    dots$include <- FALSE
  }
  P <- stan_input$sampling_vars$P
  if (is.null(dots$pars) && P > 0 && backend == "rstan") {
    # many more, but depends on the channel names
    dots$pars <- c("omega_raw_psi", "L_lf")
    dots$include <- FALSE
  }
  dots
}

#' Access the Model Formula of a Dynamite Model
#'
#' The `formula` method returns the model definition as a quoted expression.
#'
#' @rdname dynamite
#' @param x \[`dynamitefit`\]\cr The model fit object.
#' @return `formula` returns a quoted expression.
#' @family formulas
#' @export
#' @examples
#' formula(gaussian_example_fit)
#'
formula.dynamitefit <- function(x, ...) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  formula_str <- vapply(
    get_originals(x$dformulas$all),
    deparse1,
    character(1L)
  )
  ch_stoch <- which_stochastic(x$dformulas$all)
  ch_det <- which_deterministic(x$dformulas$all)
  family_str <- vapply(
    get_families(x$dformulas$all),
    function(y) y$name,
    character(1L)
  )
  lag_def <- attr(x$dformulas$all, "lags")
  spline_def <- attr(x$dformulas$stoch, "splines")
  obs_str <- onlyif(
    length(ch_stoch) > 0L,
    obs_str <- paste0(
      glue::glue(
        "obs({formula_str[ch_stoch]}, family = '{family_str[ch_stoch]}')"
      ),
      collapse = " +\n"
    )
  )
  aux_str <- onlyif(
    length(ch_det) > 0L,
    aux_str <- paste0(
      glue::glue("aux({formula_str[ch_det]})"),
      collapse = " +\n"
    )
  )
  lags_str <- onlyif(
    !is.null(lag_def),
    glue::glue("lags(k = {lag_def$k}, type = {lag_def$type})")
  )
  spline_str <- onlyif(
    spline_def$has_splines,
    paste0(
      "splines(",
      "shrinkage = ", spline_def$shrinkage, ", ",
      "override = FALSE, ",
      "df = ", spline_def$bs_opts$df, ", ",
      "degree = ", spline_def$bs_opts$degree, ", ",
      "lb_tau = ", spline_def$lb, ", ",
      "noncentered = ", spline_def$noncentered, ")"
    )
  )
  lfactor_def <- attr(x$dformulas$stoch, "lfactor")
  lfactor_str <- onlyif(
    lfactor_def$has_lfactor,
    paste0(
      "lfactor(",
      "responses = ", lfactor_def$responses, ", ",
      "noncentered_psi = ", lfactor_def$noncentered_psi, ", ",
      "nonzero_lambda = ", lfactor_def$nonzero_lambda, ", ",
      "correlated = ", lfactor_def$correlated, ")"
    )
  )
  random_spec_def <- attr(x$dformulas$stoch, "random_spec")
  random_spec_str <- onlyif(
    random_spec_def$has_random_spec,
    paste0(
      "random_spec(",
      "correlated = ", random_spec_def$correlated, ", ",
      "noncentered = ", random_spec_def$noncentered, ")"
    )
  )
  str2lang(
    paste(
      c(
        obs_str,
        aux_str,
        lags_str,
        spline_str,
        lfactor_str,
        random_spec_str
      ),
      collapse = " + "
    )
  )
}

#' Is The Argument a `dynamitefit` Object
#'
#' @param x An \R object.
#' @noRd
is.dynamitefit <- function(x) {
  inherits(x, "dynamitefit")
}

#' Parse Data for Model Fitting
#'
#' @inheritParams dynamite
#' @param group_var \[`character(1)`] Grouping variable name.
#' @param time_var \[`character(1)`] Time index variable name.
#' @srrstats {G2.4d, G2.5} Factors and ordered factors are considered.
#' @srrstats {G2.6} Columns are preprocessed.
#' @srrstats {G2.11} Data checks rely on type, not class, except for factors.
#' @srrstats {G2.12} List columns are not supported.
#' @srrstats {G2.16} Non-finite values are not supported.
#' @noRd
parse_data <- function(dformula, data, group_var, time_var, verbose) {
  data <- droplevels(data)
  # droplevels creates a copy already so we can just convert
  data <- data.table::as.data.table(data)
  data_names <- names(data)
  stopifnot_(
    !is.character(data[[time_var]]),
    c(
      "Time index variable {.arg {time_var}} is of type {.cls character}:",
      `i` = "It must be of type {.cls numeric} or {.cls factor}."
    )
  )
  if (is.factor(data[[time_var]])) {
    if (verbose) {
      warning_(c(
        "Time index variable {.arg {time_var}} is a {.cls factor}:",
        `i` = "Converting the variable to {.cls integer} based on its levels."
      ))
    }
    data.table::set(data, j = time_var, value = as.integer(data[[time_var]]))
  }
  chr_cols <- names(data)[vapply(data, is.character, logical(1L))]
  if (length(chr_cols) > 0L) {
    data[, (chr_cols) := lapply(.SD, as.factor), .SDcols = chr_cols]
  }
  valid_types <- c("integer", "logical", "double")
  col_types <- vapply(data, typeof, character(1L))
  factor_cols <- vapply(data, is.factor, logical(1L))
  valid_cols <- (col_types %in% valid_types) | factor_cols
  stopifnot_(
    all(valid_cols),
    c(
      "Column{?s} {.var {data_names[!valid_cols]}} of {.arg data}
      {?is/are} invalid:",
      `x` = "Column type{?s} {.cls {col_types[!valid_cols]}}
             {?is/are} not supported."
    )
  )
  for (j in which(valid_cols & !factor_cols)) {
    type <- typeof(data[[j]])
    col <- data[[j]]
    val <- do.call(paste0("as.", type), args = list(col))
    data.table::set(data, j = j, value = val)
  }
  resp <- get_responses(dformula)
  ordered_factor_resp <- vapply(
    seq_along(resp),
    function(i) {
      is_categorical(dformula[[i]]$family) && is.ordered(data[[resp[i]]])
    },
    logical(1L)
  )
  if (any(ordered_factor_resp)) {
    rof <- resp[ordered_factor_resp]
    if (verbose) {
      warning_(c(
        "Response variable{?s} {.var {rof}} {?is/are} of class
        {.cls ordered factor} whose channel{?s} {?is/are} categorical:",
        `i` = "{.var {rof}} will be converted to {?an/} unordered factor{?s}."
      ))
    }
    for (i in seq_along(rof)) {
      class(data[[rof[i]]]) <- "factor"
    }
  }
  finite_cols <- vapply(
    data,
    function(x) all(is.finite(x) | is.na(x)),
    logical(1L)
  )
  stopifnot_(
    all(finite_cols),
    "Non-finite values were found in variable{?s}
    {.var {data_names[!finite_cols]}} of {.arg data}."
  )
  data <- fill_time(data, group_var, time_var)
  drop_unused(dformula, data, group_var, time_var)
  data.table::setkeyv(data, c(group_var, time_var))
  data
}

#' Evaluate Past Value Definitions of Each Deterministic Channel
#'
#' @inheritParams parse_data
#' @noRd
parse_past <- function(dformula, data, group_var, time_var) {
  past <- has_past(dformula)
  for (i in seq_along(dformula)) {
    if (past[i]) {
      y <- dformula[[i]]$response
      cl <- dformula[[i]]$specials$past
      if (identical(typeof(cl), "language")) {
        past_eval <- try(eval(cl), silent = TRUE)
        if (inherits(past_eval, "try-error")) {
          past_eval <- try(data[, eval(cl)], silent = TRUE)
          # past_eval <- try(data[, cl, env = list(cl = cl)], silent = TRUE)
          stopifnot_(
            !inherits(past_eval, "try-error"),
            c(
              "Unable to evaluate past definition of
              deterministic channel {.var {y}}:",
              `x` = attr(past_eval, "condition")$message
            )
          )
        }
      } else {
        past_eval <- cl
      }
      past_type <- dformula[[i]]$specials$past_type
      past_len <- length(past_eval)
      stopifnot_(
        identical(past_type, "init") || identical(past_len, nrow(data)),
        c(
          "Incompatible past definition of deterministic channel {.var {y}}:",
          `x` = "The definition evaluates to length {past_len} but the data has
                {nrow(data)} rows."
        )
      )
      dformula[[i]]$specials$past <- past_eval
    } else {
      dformula[[i]]$specials$past_type <- "init"
    }
  }
  dformula
}

#' Parse Additional Model Formula Components
#'
#' @inheritParams parse_data
#' @param dformulas \[`list()`]\cr Output of `parse_lags`.
#' @noRd
parse_components <- function(dformulas, data, group_var, time_var) {
  fixed <- attr(dformulas$all, "max_lag")
  resp <- get_responses(dformulas$stoch)
  families <- unlist(get_families(dformulas$stoch))
  attr(dformulas$stoch, "splines") <- parse_splines(
    spline_def = attr(dformulas$stoch, "splines"),
    resp = resp,
    times = seq.int(fixed + 1L, n_unique(data[[time_var]]))
  )
  M <- count_random_effects(dformulas$stoch, data)
  stopifnot_(
    n_unique(data[[group_var]]) > 1L || M == 0L,
    "Cannot estimate random effects using only one group."
  )
  attr(dformulas$stoch, "random_spec") <- parse_random_spec(
    random_spec_def = attr(dformulas$stoch, "random_spec"),
    M = M
  )
  attr(dformulas$stoch, "lfactor") <- parse_lfactor(
    lfactor_def = attr(dformulas$stoch, "lfactor"),
    resp = resp,
    families = families
  )
  if (!attr(dformulas$stoch, "lfactor")$has_lfactor) {
    for (i in seq_along(dformulas$stoch)) {
      y <- dformulas$stoch[[i]]$response
      empty <- length(dformulas$stoch[[i]]$random) == 0 &&
        length(dformulas$stoch[[i]]$varying) == 0 &&
        length(dformulas$stoch[[i]]$fixed) == 0 &&
        !dformulas$stoch[[i]]$has_random_intercept &&
        !dformulas$stoch[[i]]$has_fixed_intercept &&
        !dformulas$stoch[[i]]$has_varying_intercept
      stopifnot_(
        !empty,
        c(
          "Invalid formula for response variable {.var {y}}:",
          `x` = "There are no predictors, intercept terms, or latent factors."
        )
      )
    }
  }
  stopifnot_(
    n_unique(data[[group_var]]) > 1L ||
      attr(dformulas$stoch, "lfactor")$P == 0L,
    "Cannot estimate latent factors using only one group."
  )
  if (attr(dformulas$stoch, "lfactor")$has_lfactor) {
    nz <- which(attr(dformulas$stoch, "lfactor")$nonzero_lambda)
    if (length(nz) > 0L) {
      lresp <- attr(dformulas$stoch, "lfactor")$responses
      for (i in nz) {
        j <- which(resp %in% lresp[i])
        vicpt <- dformulas$stoch[[j]]$has_varying_intercept
        if (vicpt) {
          dformulas$stoch[[j]]$has_varying_intercept <- FALSE
          warning_(
            "The common time-varying intercept term of channel
            {.var {lresp[i]}} was removed as channel predictors
            contain latent factor specified with {.arg nonzero_lambda}
            as TRUE.")
        }
        ficpt <- dformulas$stoch[[j]]$has_fixed_intercept
        ricpt <- dformulas$stoch[[j]]$has_random_intercept
        if (ficpt && ricpt) {
          dformulas$stoch[[j]]$has_fixed_intercept <- FALSE
          warning_(
            "The common time-invariant intercept term of channel
            {.var {lresp[i]}} was removed as channel predictors
            contain random intercept and latent factor specified
            with {.arg nonzero_lambda} as TRUE.")
        }
      }

    }
  }
  dformulas
}

#' Parse B-spline Parameters
#'
#' @param spline_def A `splines` object.
#' @param n_channels Number of channels
#' @param times An `integer` vector of time indices
#' @noRd
parse_splines <- function(spline_def, resp, times) {
  out <- list()
  n_channels <- length(resp)
  if (!is.null(spline_def)) {
    out$has_splines <- TRUE
    out$shrinkage <- spline_def$shrinkage
    out$bs_opts <- spline_def$bs_opts
    out$bs_opts$x <- times
    if (is.null(out$bs_opts$Boundary.knots)) {
      out$bs_opts$Boundary.knots <- range(out$bs_opts$x)
    }
    out$Bs <- t(do.call(splines::bs, args = out$bs_opts))
    out$D <- nrow(out$Bs)
    out$lb <- spline_def$lb_tau
    if (length(out$lb) %in% c(1L, n_channels)) {
      out$lb <- rep(out$lb, length = n_channels)
    } else {
      stop_(
        "Length of the {.arg lb_tau} argument of {.fun splines} function
        is not equal to 1 or {n_channels}, the number of the channels."
      )
    }
    out$noncentered <- spline_def$noncentered
    if (length(out$noncentered) %in% c(1L, n_channels)) {
      out$noncentered <- rep(out$noncentered, length = n_channels)
    } else {
      stop_(
        "Length of the {.arg noncentered} argument of {.fun splines} function
        is not equal to 1 or {n_channels}, the number of the channels."
      )
    }
  } else {
    out <- list(
      has_splines = FALSE,
      lb = numeric(n_channels),
      noncentered = logical(n_channels),
      shrinkage = logical(1L)
    )
  }
  out
}

#' Parse Random Effect Definitions
#'
#' @param random_spec_def A `random_spec` object.
#' @param M Number of random effects.
#' @noRd
parse_random_spec <- function(random_spec_def, M) {
  out <- list()
  if (is.null(random_spec_def)) {
    out <- list(
      has_random_spec = FALSE,
      correlated = TRUE,
      noncentered = TRUE
    )
  } else {
    out <- random_spec_def
    out$has_random_spec <- TRUE
  }
  if (M < 2L) {
    out$correlated <- FALSE
  }
  out$M <- M
  out
}

#' Parse Latent Factor Definitions
#'
#' @param lfactor_def An `lfactor` object.
#' @param resp A `character` vector of response variable names.
#' @param families A `character` vector response family names.
#' @noRd
parse_lfactor <- function(lfactor_def, resp, families) {
  out <- list()
  if (!is.null(lfactor_def)) {
    valid_channels <- resp
    # default, use all channels except categorical
    if (is.null(lfactor_def$responses)) {
      lfactor_def$responses <- valid_channels
    } else {
      psi_channels <- lfactor_def$responses %in% resp
      stopifnot_(
        all(psi_channels),
        c(
          "Argument {.arg responses} of {.fun lfactor} contains variable{?s}
          {.var {cs(lfactor_def$responses[!psi_channels])}}:",
          `x` = "No such response variables in the model."
        )
      )
    }
    if (length(lfactor_def$responses) < 2L) {
      lfactor_def$correlated <- FALSE
    }
    out$has_lfactor <- TRUE
    out$responses <- lfactor_def$responses
    out$noncentered_psi <- lfactor_def$noncentered_psi
    n_channels <- length(resp)
    out$nonzero_lambda <- lfactor_def$nonzero_lambda
    stopifnot_(
      length(out$nonzero_lambda) %in% c(1L, n_channels),
      "Length of the {.arg nonzero_lambda} argument of {.fun lfactor} function
      is not equal to 1 or {n_channels}, the number of the channels."
    )
    out$nonzero_lambda <- rep(out$nonzero_lambda, length = n_channels)
    out$correlated <- lfactor_def$correlated
    out$P <- length(lfactor_def$responses)
  } else {
    n_channels <- length(resp)
    out <- list(
      has_lfactor = FALSE,
      responses = character(0L),
      noncentered_psi = FALSE,
      nonzero_lambda = logical(n_channels),
      correlated = FALSE,
      P = 0
    )
  }
  out
}

#' Adds NA Gaps to Fill In Missing Time Points in a Data Frame
#'
#' @inheritParams dynamite
#' @noRd
fill_time <- function(data, group_var, time_var) {
  time <- sort(unique(data[[time_var]]))
  stopifnot_(
    length(time) > 1L,
    "There must be at least two time points in the data."
  )
  # time_duplicated <- data[,
  #  any(duplicated(time_var)),
  #  by = group_var,
  #  env = list(time_var = time_var)
  # ]$V1
  split_data <- split(data, by = group_var)
  time_duplicated <- unlist(
    lapply(split_data, function(x) any(duplicated(x[[time_var]])))
  )
  d <- which(time_duplicated)
  stopifnot_(
    all(!time_duplicated),
    c(
      "Each time index must correspond to a single observation per group:",
      `x` = "{cli::qty(length(d))}Group{?s} {.var {d}} of {.var {group_var}}
             {cli::qty(length(d))}{?has/have} duplicate observations."
    )
  )
  time_ivals <- diff(time)
  time_scale <- min(diff(time))
  stopifnot_(
    all(time_ivals[!is.na(time_ivals)] %% time_scale == 0),
    "Observations must occur at regular time intervals."
  )
  full_time <- seq(time[1L], time[length(time)], by = time_scale)
  # time_missing <- data[,
  #  !identical(time_var, full_time),
  #  by = group_var,
  #  env = list(time_var = time_var, full_time = full_time)
  # ]$V1
  time_missing <- unlist(
    lapply(split_data, function(x) !identical(x[[time_var]], full_time))
  )
  if (any(time_missing)) {
    full_data_template <- data.table::as.data.table(
      expand.grid(
        time = full_time,
        group = unique(data[[group_var]])
      )
    )
    names(full_data_template) <- c(time_var, group_var)
    data <- data.table::merge.data.table(
      full_data_template,
      data,
      by = c(time_var, group_var),
      all.x = TRUE
    )
  }
  data
}

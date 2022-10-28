#' Extract Samples From a `dynamitefit` Object as a Data Table
#'
#' Provides a `data.table` representation of the posterior samples of the model
#' parameters. See [dynamite::as.data.frame.dynamitefit()] for details.
#'
#' @inheritParams as.data.frame.dynamitefit
#' @export
#' @examples
#' as.data.table(
#'   gaussian_example_fit,
#'   responses = "y",
#'   types = "beta",
#'   summary = FALSE
#' )
#'
as.data.table.dynamitefit <- function(x, row.names = NULL, optional = FALSE,
                                      responses = NULL, types = NULL,
                                      summary = TRUE, probs = c(0.05, 0.95),
                                      include_fixed = TRUE, ...) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  stopifnot_(
    !is.null(x$stanfit),
    "No Stan model fit is available."
  )
  stopifnot_(
    checkmate::test_character(
      x = responses,
      any.missing = FALSE,
      min.len = 1L,
      null.ok = TRUE
    ),
    "Argument {.arg responses} must be a {.cls character} vector."
  )
  stopifnot_(
    checkmate::test_character(
      x = types,
      any.missing = FALSE,
      min.len = 1L,
      null.ok = TRUE
    ),
    "Argument {.arg types} must be a {.cls character} vector."
  )
  stopifnot_(
    checkmate::test_flag(x = summary),
    "Argument {.arg summary} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_numeric(
      x = probs,
      lower = 0.0,
      upper = 1.0,
      any.missing = FALSE,
      min.len = 1L
    ),
    "Argument {.arg probs} must be a {.cls numeric} vector with values between
     0 and 1."
  )
  stopifnot_(
    checkmate::test_flag(x = include_fixed),
    "Argument {.arg include_fixed} must be a single {.cls logical} value."
  )
  if (is.null(responses)) {
    responses <- setdiff(unique(x$priors$response), "")
  } else {
    z <- responses %in% unique(x$priors$response)
    stopifnot_(
      all(z),
      "Model does not contain response variable{?s} {.var {responses[!z]}}."
    )
  }
  all_types <- c(
    "alpha", "beta", "delta", "tau", "tau_alpha", "xi",
    "sigma_nu", "corr_nu", "sigma", "phi", "nu", "omega", "omega_alpha"
  )
  if (is.null(types)) {
    types <- all_types[seq_len(11L)]
  } else {
    types <- onlyif(is.character(types), tolower(types))
    types <- try(match.arg(types, all_types, TRUE), silent = TRUE)
    stopifnot_(
      !inherits(types, "try-error"),
      "Argument {.arg type} contains unknown types."
    )
  }
  values <- function(type, response) {
    if (type %in% c("xi", "corr_nu")) {
      draws <- rstan::extract(
        x$stanfit,
        pars = type,
        permuted = FALSE
      )
    } else {
      draws <- rstan::extract(
        x$stanfit,
        pars = paste0(type, "_", response),
        permuted = FALSE
      )
    }
    category <- attr(x$stan$responses[[response]], "levels")[-1L]
    if (is.null(category)) {
      category <- NA
    }
    d <- do.call(
      what = paste0("as_data_table_", type),
      args = list(
        x = x,
        draws = draws,
        n_draws = prod(dim(draws)[1L:2L]),
        response = response,
        category = category,
        include_fixed = include_fixed
      )
    )
    n_d <- d[, .N]
    n_r <- nrow(draws)
    n_c <- ncol(draws)
    d[, response := rep(response, n_d)]
    d[, type := rep(type, n_d)]
    d[, .draw := rep_len(seq_len(n_r * n_c), n_d)]
    d[, .iteration := rep_len(seq_len(n_r), n_d)]
    d[, .chain := rep_len(rep(seq_len(n_c), each = n_r), n_d)]
    d
  }
  # avoid NSE notes from R CMD check
  .chain <- .draw <- .iteration <- NULL
  category <- group <- parameter <- response <- time <- type <- value <- NULL
  out_all <- NULL
  if ("xi" %in% types) {
    out_all <- data.table::data.table(
      type = "xi",
      response = "",
      parameter = "xi"
    )
  }
  if ("corr_nu" %in% types) {
    out_all <- rbind(
      out_all,
      data.table::data.table(
        type = "corr_nu",
        response = "",
        parameter = "corr_nu"
      )
    )
  }
  tmp <- data.table::as.data.table(
    expand.grid(
      type = types,
      response = responses,
      stringsAsFactors = FALSE
    )
  )
  tmp[, parameter := as.character(glue::glue("{tmp$type}_{tmp$response}"))]
  out <- data.table::rbindlist(list(out_all, tmp))
  rows <- apply(out, 1L, function(y) {
    any(
      grepl(
        paste0("^", y["parameter"]),
        x$stanfit@sim$pars_oi
      )
    )
  })
  out <- out[rows, c("response", "type")]
  n_pars <- nrow(out)
  stopifnot_(
    n_pars > 0L,
    paste0(
      "No parameters of type {.var ",
      paste(types, collapse = "}, {.var "),
      "} found for any of the response channels {.var ",
      paste(responses, collapse = "}, {.var "), "}."
    )
  )
  all_values <- vector(mode = "list", length = n_pars + 1L)
  # template for rbindlist
  all_values[[1L]] <- data.table::data.table(
    parameter = character(0L),
    value = numeric(0L),
    time = x$data[[x$time_var]][0L],
    category = character(0L),
    group = x$data[[x$group_var]][0L],
    response = character(0L),
    type = character(0L),
    .draw = integer(0L),
    .iteration = integer(0L),
    .chain = integer(0L)
  )
  all_values[seq.int(2L, n_pars + 1L)] <- .mapply(
    values,
    dots = list(type = out$type, response = out$response),
    MoreArgs = NULL
  )
  out <- data.table::rbindlist(all_values, fill = TRUE)
  if (summary) {
    pars <- unique(out$parameter)
    out <- out[,
      parameter := factor(parameter, levels = pars, ordered = TRUE),
      env = list(parameter = "parameter")
    ][,
      {
        mean = mean(value)
        sd = sd(value)
        tmp = quantile(value, na.rm = TRUE)
        q5 = tmp[1L]
        q95 = tmp[2L]
        list(mean = mean, sd = sd, q5 = q5, q95 = q95)
      },
      by = list(parameter, time, category, group, response, type)
    ][,
      parameter := as.character(parameter)
    ]
  }
  out
}

#' Construct a Data Table for a Parameter Type from a `dynamitefit` Object
#'
#' Arguments for all as_data_frame_type functions are documented here.
#'
#' @inheritParams as.data.frame.dynamitefit
#' @param draws \[`list()`]\cr A Stan fit draws object.
#' @param n_draws \[`integer(1)`]\cr Number of draws.
#' @param response \[`character(1)`]\cr Response variable name.
#' @param categories \[`character()`]\cr Levels of categorical responses.
#' @noRd
as_data_table_default <- function(type, draws, response) {
  data.table::data.table(
    parameter = paste0(type, "_", response),
    value = c(draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "xi" Parameter
#' @noRd
as_data_table_xi <- function(x, draws, ...) {
  data.table::data.table(
    parameter = "xi",
    value = c(draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "corr_nu" Parameter
#' @noRd
as_data_table_corr_nu <- function(x, draws, ...) {
  resp <- get_responses(x$dformulas$stoch)
  pairs <- apply(utils::combn(resp, 2L), 2L, paste, collapse = "_")
  data.table::data.table(
    parameter = paste0("corr_nu_", pairs),
    value = c(draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "nu" Parameter
#' @noRd
as_data_table_nu <- function(x, draws, n_draws, response, ...) {
  n_group <- dim(draws)[3L]
  data.table::data.table(
    parameter = paste0("nu_", response),
    value = c(draws),
    group = rep(sort(unique(x$data[[x$group_var]])), each = n_draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "alpha" Parameter
#' @noRd
as_data_table_alpha <- function(x, draws, n_draws,
                                response, category, include_fixed) {
  n_cat <- length(category)
  fixed <- x$stan$fixed
  all_time_points <- sort(unique(x$data[[x$time_var]]))
  if (x$stan$model_vars[[response]]$has_varying_intercept) {
    time_points <- ifelse_(
      include_fixed,
      all_time_points,
      all_time_points[seq.int(fixed + 1L, length(all_time_points))]
    )
    n_na <- include_fixed * fixed * n_draws
    n_time <- length(time_points)
    n_time2 <- n_time - include_fixed * fixed
    data.table::rbindlist(
      lapply(seq_len(n_cat), function(i) {
        data.table::data.table(
          parameter = paste0("alpha_", response),
          value = c(
            rep(NA, n_na),
            c(draws[, , (i - 1L) * n_time2 + seq_len(n_time2)])
          ),
          time = rep(time_points, each = n_draws),
          category = category[i]
        )
      })
    )
  } else {
    data.table::data.table(
      parameter = paste0("alpha_", response),
      value = c(draws),
      category = rep(category, each = n_draws)
    )
  }
}

#' @describeIn as_data_table_default Data Table for a "beta" Parameter
#' @noRd
as_data_table_beta <- function(x, draws, n_draws, response, category, ...) {
  var_names <- paste0(
    "beta_", response, "_",
    names(x$stan$model_vars[[response]]$J_fixed)
  )
  n_vars <- length(var_names)
  data.table::data.table(
    parameter = rep(var_names, each = n_draws),
    value = c(draws),
    category = rep(category, each = n_vars * n_draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "delta" Parameter
#' @noRd
as_data_table_delta <- function(x, draws, n_draws,
                                response, category, include_fixed) {
  n_cat <- length(category)
  fixed <- x$stan$fixed
  all_time_points <- sort(unique(x$data[[x$time_var]]))
  var_names <- paste0(
    "delta_", response, "_",
    names(x$stan$model_vars[[response]]$J_varying)
  )
  n_vars <- length(var_names)
  time_points <- ifelse_(
    include_fixed,
    all_time_points,
    all_time_points[seq.int(fixed + 1L, length(all_time_points))]
  )
  n_na <- include_fixed * fixed * n_draws
  n_time <- length(time_points)
  n_time2 <- n_time - include_fixed * fixed
  data.table::rbindlist(lapply(seq_len(n_cat), function(j) {
    data.table::rbindlist(lapply(seq_len(n_vars), function(i) {
      idx <- (j - 1L) * n_time2 * n_vars + (i - 1L) * n_time2 + seq_len(n_time2)
      data.table::data.table(
        parameter = var_names[i],
        value = c(
          rep(NA, n_na),
          c(draws[, , idx])
        ),
        time = rep(time_points, each = n_draws),
        category = rep(category[j], each = n_time * n_draws)
      )
    }))
  }))
}

#' @describeIn as_data_table_default Data Table for a "tau" Parameter
#' @noRd
as_data_table_tau <- function(x, draws, n_draws, response, ...) {
  var_names <- paste0(
    "tau_", response, "_",
    names(x$stan$model_vars[[response]]$J_varying)
  )
  data.table::data.table(
    parameter = rep(var_names, each = n_draws),
    value = c(draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "omega" Parameter
#' @noRd
as_data_table_omega <- function(x, draws, n_draws, response, category, ...) {
  n_cat <- length(category)
  D <- x$stan$sampling_vars$D
  var_names <- names(x$stan$model_vars[[response]]$J_varying)
  k <- length(var_names)
  data.table::data.table(
    parameter = rep(
      paste0(
        "omega_", rep(seq_len(D), each = n_cat * k), "_",
        rep(var_names, each = n_cat)
      ),
      each = n_draws
    ),
    value = c(draws),
    category = rep(category, each = n_draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "omega_alpha" Parameter
#' @noRd
as_data_table_omega_alpha <- function(x, draws, n_draws, category, ...) {
  n_cat <- length(category)
  D <- x$stan$sampling_vars$D
  data.table::data.table(
    parameter = rep(
      paste0("omega_alpha_", seq_len(D)),
      each = n_cat * n_draws
    ),
    value = c(draws),
    category = rep(category, each = n_draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "tau_alpha" Parameter
#' @noRd
as_data_table_tau_alpha <- function(draws, response, ...) {
  as_data_table_default("tau_alpha", draws, response)
}

#' @describeIn as_data_table_default Data Table for a "sigma" Parameter
#' @noRd
as_data_table_sigma <- function(draws, response, ...) {
  as_data_table_default("sigma", draws, response)
}

#' @describeIn as_data_table_default Data Table for a "sigma_nu" Parameter
#' @noRd
as_data_table_sigma_nu <- function(draws, response, ...) {
  as_data_table_default("sigma_nu", draws, response)
}

#' @describeIn as_data_table_default Data Table for a "phi" Parameter
#' @noRd
as_data_table_phi <- function(draws, response, ...) {
  as_data_table_default("phi", draws, response)
}

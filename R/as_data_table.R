#' Extract Samples From a `dynamitefit` Object as a Data Table
#'
#' Provides a `data.table` representation of the posterior samples of the model
#' parameters. See [dynamite::as.data.frame.dynamitefit()] for details.
#'
#' @export
#' @export as.data.table
#' @family output
#' @aliases as.data.table
#' @importFrom data.table as.data.table
#' @param keep.rownames \[`logical(1)`]\cr Not used.
#' @inheritParams as.data.frame.dynamitefit
#' @return A `data.table` containing either samples or summary statistics of
#'   the model parameters.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' as.data.table(
#'   gaussian_example_fit,
#'   responses = "y",
#'   types = "beta",
#'   summary = FALSE
#' )
#'
as.data.table.dynamitefit <- function(x, keep.rownames = FALSE,
                                      row.names = NULL, optional = FALSE,
                                      types = NULL, parameters = NULL,
                                      responses = NULL,
                                      times = NULL, groups = NULL,
                                      summary = FALSE, probs = c(0.05, 0.95),
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
      x = parameters,
      any.missing = FALSE,
      min.len = 1L,
      null.ok = TRUE
    ),
    "Argument {.arg parameters} must be a {.cls character} vector."
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
  types <- onlyif(is.character(types), tolower(types))
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
    checkmate::test_numeric(
      x = times,
      min.len = 1L,
      null.ok = TRUE
    ),
    "Argument {.arg times} must be a {.cls integer} vector."
  )
  stopifnot_(
    checkmate::test_vector(
      x = groups,
      min.len = 1L,
      null.ok = TRUE
    ),
    "Argument {.arg groups} must be a vector."
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
  if (!is.null(parameters)) {
    responses <- types <- NULL
  }
  all_responses <- unique(c(names(x$stan$responses), unlist(x$stan$responses)))
  if (is.null(responses)) {
    responses <- all_responses
  } else {
    valid_responses <- responses %in% all_responses
    stopifnot_(
      all(valid_responses),
      c(
        "Argument {.arg responses} contains invalid response variable names.",
        `x` = "Response variable{?s} {.val {responses[!valid_responses]}}
               {?is/are} not recognized.",
        `i` = "The response variable{?s} of the model
               {?is/are} {.val {all_responses}}."
      )
    )
  }
  if (is.null(types)) {
    types <- ifelse_(
      is.null(parameters),
      all_types[!grepl("omega", all_types, fixed = TRUE)],
      all_types
    )
  } else {
    match_types <- match(types, all_types)
    valid_types <- !is.na(match_types)
    stopifnot_(
      all(valid_types),
      c(
        "Argument {.arg types} contains invalid types.",
        `x` = "Type{?s} {.val {types[!valid_types]}} {?is/are} not recognized.",
        `i` = "Use {.fun get_parameter_types} to check available types."
      )
    )
  }
  values <- function(type, response, category) {
    ycat <- ifelse_(
      nzchar(category) & !is.na(category),
      paste0("_", category),
      ""
    )
    if (type %in% c("xi", "corr_nu", "corr_psi")) {
      draws <- rstan::extract(
        x$stanfit,
        pars = type,
        permuted = FALSE
      )
    } else {
      draws <- rstan::extract(
        x$stanfit,
        pars = paste0(type, "_", response, ycat),
        permuted = FALSE
      )
    }
    channel <- get_channel(x, response)
    idx <- which(names(x$stan$responses) %in% response)
    resps <- ifelse_(
      identical(length(idx), 0L),
      NULL,
      x$stan$responses[[idx]]
    )
    d <- do.call(
      what = paste0("as_data_table_", type),
      args = list(
        x = x,
        draws = draws,
        n_draws = prod(dim(draws)[1L:2L]),
        response = response,
        category = category,
        include_fixed = include_fixed,
        resps = resps
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
    if (!is.null(times)) {
      d <- d[time %in% times, , env = list(times = times)]
    }
    if (!is.null(groups)) {
      d <- d[group %in% groups, , env = list(groups = I(groups))]
    }
    d
  }
  # avoid NSE notes from R CMD check
  .chain <- .draw <- .iteration <- NULL
  category <- group <- parameter <- response <- NULL
  catstr <- time <- type <- value <- NULL
  out_all <- NULL
  if ("xi" %in% types) {
    out_all <- data.table::data.table(
      type = "xi",
      response = "",
      category = NA_character_,
      parameter = "xi"
    )
  }
  if ("corr_nu" %in% types) {
    out_all <- rbind(
      out_all,
      data.table::data.table(
        type = "corr_nu",
        response = "",
        category = NA_character_,
        parameter = "corr_nu"
      )
    )
  }
  if ("corr_psi" %in% types) {
    out_all <- rbind(
      out_all,
      data.table::data.table(
        type = "corr_psi",
        response = "",
        category = NA_character_,
        parameter = "corr_psi"
      )
    )
  }
  categories <- unique(
    c(
      NA_character_,
      ulapply(
        unlist(x$stan$responses),
        function(y) {
          channel <- get_channel(x, y)
          if (is_cumulative(channel$family)) {
            seq_len(channel$S - 1L)
          } else if (is_categorical(channel$family)) {
            channel$categories[-1L]
          } else {
            NA_character_
          }
        }
      )
    )
  )
  tmp <- data.table::as.data.table(
    expand.grid(
      type = types,
      response = responses,
      category = categories,
      stringsAsFactors = FALSE
    )
  )
  tmp[, catstr := ifelse(
    nzchar(category) & !is.na(category),
    glue::glue("_{tmp$category}"),
    ""
  )]
  tmp[, parameter := as.character(
    glue::glue("{tmp$type}_{tmp$response}{tmp$catstr}")
  )]
  tmp[, catstr := NULL]
  out <- data.table::rbindlist(list(out_all, tmp))
  rows <- apply(out, 1L, function(y) {
    any(
      grepl(
        paste0("^", y["parameter"], "$"),
        x$stanfit@sim$pars_oi
      )
    )
  })
  out <- out[rows, c("response", "category", "type")]
  n_pars <- nrow(out)
  stopifnot_(
    n_pars > 0L,
    "No parameters of type {.var {types}} were found for any of the response
     channels {.var {responses}}."
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
    dots = list(
      type = out$type,
      response = out$response,
      category = out$category
    ),
    MoreArgs = NULL
  )
  out <- data.table::rbindlist(all_values, fill = TRUE)
  if (!is.null(parameters)) {
    data.table::setkey(out, "parameter")
    valid_pars <- parameters %in% unique(out$parameter)
    stopifnot_(
      all(valid_pars),
      c(
        "Argument {.arg parameters} contains invalid parameter names.",
        `x` = "Parameter{?s} {.val {parameters[!valid_pars]}} {?is/are}
               not recognized.",
        `i` = "Use {.fun get_parameter_names} to check available parameters."
      )
    )
    out <- out[parameters]
  }
  if (summary) {
    pars <- unique(out$parameter)
    out <- out[,
      parameter := factor(parameter, levels = pars, ordered = TRUE)
    ][,
      {
        mean <- mean(value)
        sd <- sd(value)
        tmp <- quantile(value, probs = probs, na.rm = TRUE)
        names(tmp) <- paste0("q", 100 * probs)
        c(list(mean = mean, sd = sd), tmp)
      },
      by = list(parameter, time, group, category, response, type)
    ][,
      parameter := as.character(parameter)
    ]
    pnames <- c("time", "group", "category", "response", "type")
    cnames <- setdiff(colnames(out), pnames)
    data.table::setcolorder(out, neworder = c(cnames, pnames))
  }
  out
}

#' Get Channel Variables or Channel Group Variables
#'
#' @param x A `dynamitefit` object.
#' @param response The response variable name.
#' @noRd
get_channel <- function(x, response) {
  if (is.null(x$stan$channel_vars[[response]])) {
    x$stan$channel_group_vars[[response]]
  } else {
    x$stan$channel_vars[[response]]
  }
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
as_data_table_default <- function(type, draws, response, ...) {
  data.table::data.table(
    parameter = paste0(type, "_", response),
    value = c(draws)
  )
}

#' Shrinkage feature removed at least for now.
#'
#' @describeIn as_data_table_default Data Table for a "xi" Parameter
#' @noRd
# as_data_table_xi <- function(x, draws, n_draws, ...) {
#   D <- x$stan$model_vars$D
#   data.table::data.table(
#     parameter = rep(
#       paste0("xi_d", seq_len(D - 1L)),
#       each = n_draws
#     ),
#     value = c(draws)
#   )
# }

#' @describeIn as_data_table_default Data Table for a "corr_nu" Parameter
#' @noRd
as_data_table_corr_nu <- function(x, draws, n_draws, ...) {
  vars <- ulapply(
    x$stan$channel_vars,
    function(y) {
      if (y$has_random || y$has_random_intercept) {
        icpt <- ifelse_(
          y$has_random_intercept,
          "alpha",
          NULL
        )
        vars <- paste0(y$y, "_", c(icpt, names(y$J_random)))
        ifelse_(
          is_categorical(y$family),
          paste0(
            rep(vars, y$S - 1L),
            "_",
            rep(y$categories[-1L], each = y$K_random)
          ),
          vars
        )
      }
    }
  )
  pairs <- apply(utils::combn(vars, 2L), 2L, paste, collapse = "__")
  data.table::data.table(
    parameter = rep(paste0("corr_nu_", pairs), each = n_draws),
    value = c(draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "nu" Parameter
#' @noRd
as_data_table_nu <- function(x, draws, n_draws, response, category, ...) {
  icpt <- ifelse_(
    get_channel(x, response)$has_random_intercept,
    "alpha",
    NULL
  )
  var_names <- paste0(
    "nu_", response, "_",
    c(icpt, names(get_channel(x, response)$J_random))
  )
  n_vars <- length(var_names)
  groups <- sort(unique(x$data[[x$group_var]]))
  n_group <- length(groups)
  data.table::data.table(
    parameter = rep(var_names, each = n_draws * n_group),
    value = c(draws),
    group = rep(groups, each = n_draws),
    category = category
  )
}

#' @describeIn as_data_table_default Data Table for a "alpha" Parameter
#' @noRd
as_data_table_alpha <- function(x, draws, n_draws,
                                response, category, include_fixed, ...) {
  fixed <- x$stan$fixed
  all_time_points <- sort(unique(x$data[[x$time_var]]))
  if (get_channel(x, response)$has_varying_intercept) {
    time_points <- ifelse_(
      include_fixed,
      all_time_points,
      all_time_points[seq.int(fixed + 1L, length(all_time_points))]
    )
    n_na <- include_fixed * fixed * n_draws
    n_time <- length(time_points)
    n_time2 <- n_time - include_fixed * fixed
    data.table::data.table(
      parameter = paste0("alpha_", response),
      value = c(
        rep(NA, n_na),
        c(draws[, , seq_len(n_time2)])
      ),
      time = rep(time_points, each = n_draws),
      category = category
    )
  } else {
    data.table::data.table(
      parameter = paste0("alpha_", response),
      value = c(draws),
      category = category
    )
  }
}

#' @describeIn as_data_table_default Data Table for a "beta" Parameter
#' @noRd
as_data_table_beta <- function(x, draws, n_draws, response, category, ...) {
  var_names <- paste0(
    "beta_", response, "_",
    names(get_channel(x, response)$J_fixed)
  )
  n_vars <- length(var_names)
  data.table::data.table(
    parameter = rep(var_names, each = n_draws),
    value = c(draws),
    category = category
  )
}

#' @describeIn as_data_table_default Data Table for a "delta" Parameter
#' @noRd
as_data_table_delta <- function(x, draws, n_draws,
                                response, category, include_fixed, ...) {
  fixed <- x$stan$fixed
  all_time_points <- sort(unique(x$data[[x$time_var]]))
  var_names <- paste0(
    "delta_", response, "_",
    names(get_channel(x, response)$J_varying)
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
  data.table::rbindlist(lapply(seq_len(n_vars), function(i) {
    idx <- (i - 1L) * n_time2 + seq_len(n_time2)
    data.table::data.table(
      parameter = var_names[i],
      value = c(
        rep(NA, n_na),
        c(draws[, , idx])
      ),
      time = rep(time_points, each = n_draws),
      category = category
    )
  }))
}

#' @describeIn as_data_table_default Data Table for a "tau" Parameter
#' @noRd
as_data_table_tau <- function(x, draws, n_draws, response, category, ...) {
  var_names <- paste0(
    "tau_", response, "_",
    names(get_channel(x, response)$J_varying)
  )
  data.table::data.table(
    parameter = rep(var_names, each = n_draws),
    value = c(draws),
    category = category
  )
}

#' @describeIn as_data_table_default Data Table for a "omega" Parameter
#' @noRd
as_data_table_omega <- function(x, draws, n_draws, response, category, ...) {
  n_cat <- length(category)
  D <- x$stan$model_vars$D
  var_names <- paste0(
    "omega_", response, "_",
    names(get_channel(x, response)$J_varying)
  )
  k <- length(var_names)
  params_ord <- rep(
    paste0(rep(var_names, each = D), "_d", rep(seq_len(D), k)),
    each = n_draws
  )
  tmp <- data.table::data.table(
    parameter = rep(
      paste0(var_names, "_d", rep(seq_len(D), each = k)),
      each = n_draws
    ),
    value = c(draws),
    category = category
  )
  tmp[order(match(parameter, params_ord))]
}

#' @describeIn as_data_table_default Data Table for a "omega_alpha" Parameter
#' @noRd
as_data_table_omega_alpha <-function(x, draws, n_draws, response,
  category, ...) {
  D <- x$stan$model_vars$D
  data.table::data.table(
    parameter = rep(
      paste0("omega_alpha_", response, "_d", seq_len(D)),
      each = n_draws
    ),
    value = c(draws),
    category = category
  )
}

#' @describeIn as_data_table_default Data Table for a "tau_alpha" Parameter
#' @noRd
as_data_table_tau_alpha <- function(draws, response, category, ...) {
  data.table::data.table(
    parameter = paste0("tau_alpha_", response),
    value = c(draws),
    category = category
  )
}

#' @describeIn as_data_table_default Data Table for a "sigma" Parameter
#' @noRd
as_data_table_sigma <- function(draws, response, ...) {
  as_data_table_default("sigma", draws, response)
}

#' @describeIn as_data_table_default Data Table for a "sigma_nu" Parameter
#' @noRd
as_data_table_sigma_nu <- function(x, draws, n_draws, response, category, ...) {
  icpt <- ifelse_(
    get_channel(x, response)$has_random_intercept,
    "alpha",
    NULL
  )
  var_names <- paste0(
    "sigma_nu_", response, "_",
    c(icpt, names(get_channel(x, response)$J_random))
  )
  data.table::data.table(
    parameter = rep(var_names, each = n_draws),
    value = c(draws),
    category = category
  )
}

#' @describeIn as_data_table_default Data Table for a "phi" Parameter
#' @noRd
as_data_table_phi <- function(draws, response, ...) {
  as_data_table_default("phi", draws, response)
}

#' @describeIn as_data_table_default Data Table for a "lambda" Parameter
#' @noRd
as_data_table_lambda <- function(x, draws, n_draws, response, ...) {
  n_group <- dim(draws)[3L]
  data.table::data.table(
    parameter = paste0("lambda_", response),
    value = c(draws),
    group = rep(sort(unique(x$data[[x$group_var]])), each = n_draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "sigma_lambda" Parameter
#' @noRd
as_data_table_sigma_lambda <- function(draws, response, ...) {
  as_data_table_default("sigma_lambda", draws, response)
}

#' @describeIn as_data_table_default Data Table for a "psi" Parameter
#' @noRd
as_data_table_psi <- function(x, draws, n_draws, response,
                              category, include_fixed, ...) {
  fixed <- x$stan$fixed
  all_time_points <- sort(unique(x$data[[x$time_var]]))
  time_points <- ifelse_(
    include_fixed,
    all_time_points,
    all_time_points[seq.int(fixed + 1L, length(all_time_points))]
  )
  n_na <- include_fixed * fixed * n_draws
  n_time <- length(time_points)
  n_time2 <- n_time - include_fixed * fixed
  data.table::data.table(
    parameter = paste0("psi_", response),
    value = c(
      rep(NA, n_na),
      c(draws[, , seq_len(n_time2)])
    ),
    time = rep(time_points, each = n_draws),
    category = category
  )
}

#' @describeIn as_data_table_default Data Table for a "tau_psi" Parameter
#' @noRd
as_data_table_tau_psi <- function(draws, response, ...) {
  as_data_table_default("tau_psi", draws, response)
}

#' @describeIn as_data_table_default Data Table for a "omega_psi" Parameter
#' @noRd
as_data_table_omega_psi <- function(x, draws, n_draws, response,
                                    category, ...) {
  D <- x$stan$model_vars$D
  data.table::data.table(
    parameter = rep(
      paste0("omega_psi_", response, "_d", seq_len(D)),
      each =  n_draws
    ),
    value = c(draws),
    category = category
  )
}

#' @describeIn as_data_table_default Data Table for a "corr_psi" Parameter
#' @noRd
as_data_table_corr_psi <- function(x, draws, n_draws, ...) {
  resp <- attr(x$dformulas$stoch, "lfactor")$responses
  pairs <- apply(utils::combn(resp, 2L), 2L, paste, collapse = "__")
  data.table::data.table(
    parameter = rep(paste0("corr_psi_", pairs), each = n_draws),
    value = c(draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "corr" Parameter
#' @noRd
as_data_table_corr <- function(x, draws, n_draws, resps, ...) {
  pairs <- apply(utils::combn(resps, 2L), 2L, paste, collapse = "__")
  data.table::data.table(
    parameter = rep(paste0("corr_", pairs), each = n_draws),
    value = c(draws)
  )
}

#' @describeIn as_data_table_default Data Table for a "cutpoints" Parameter
#' @noRd
as_data_table_cutpoints <- function(x, draws, response,
                                    n_draws, include_fixed, ...) {
  channel <- get_channel(x, response)
  S <- channel$S
  fixed <- x$stan$fixed
  all_time_points <- sort(unique(x$data[[x$time_var]]))
  if (channel$has_varying_intercept) {
    time_points <- ifelse_(
      include_fixed,
      all_time_points,
      all_time_points[seq.int(fixed + 1L, length(all_time_points))]
    )
    n_na <- include_fixed * fixed * n_draws
    n_time <- length(time_points)
    n_time2 <- n_time - include_fixed * fixed
    data.table::rbindlist(lapply(seq_len(S - 1L), function(i) {
      idx <- (i - 1L) * n_time2 + seq_len(n_time2)
      data.table::data.table(
        parameter = paste0("cutpoints_", response),
        value = c(
          rep(NA, n_na),
          c(draws[, , idx])
        ),
        time = rep(time_points, each = n_draws),
        category = i
      )
    }))
  } else {
    data.table::data.table(
      parameter = paste0("cutpoints_", response),
      category = rep(seq_len(S - 1L), each = n_draws),
      value = c(draws)
    )
  }
}

# Parameter types ---------------------------------------------------------

all_types <- c(
  "alpha",
  "beta",
  "corr",
  "corr_nu",
  "corr_psi",
  "cutpoints",
  "delta",
  "lambda",
  "nu",
  "omega",
  "omega_alpha",
  "omega_psi",
  "phi",
  "psi",
  "sigma_nu",
  "sigma",
  "sigma_lambda",
  "tau",
  "tau_alpha",
  "tau_psi",
  "xi"
)

fixed_types <- c(
  "alpha",
  "beta",
  "corr",
  "corr_nu",
  "corr_psi",
  "cutpoints",
  "lambda",
  "nu",
  "omega",
  "omega_alpha",
  "omega_psi",
  "phi",
  "sigma",
  "sigma_lambda",
  "sigma_nu",
  "tau",
  "tau_alpha",
  "tau_psi",
  "xi"
)

varying_types <- c(
  "alpha",
  "cutpoints",
  "delta",
  "psi"
)

default_types <- c(
  "alpha",
  "beta",
  "cutpoints",
  "delta",
  "lambda",
  "nu",
  "psi"
)

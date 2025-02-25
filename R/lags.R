#' Add Lagged Responses as Predictors to Each Channel of a \pkg{dynamite} Model
#'
#' Adds the lagged value of the response of each channel specified via
#' [dynamiteformula()] as a predictor to each channel. The added predictors
#' can be either time-varying or time-invariant.
#'
#' @export
#' @family formulas
#' @param k \[`integer()`]\cr
#'   Values lagged by `k` units of time of each observed response variable
#'   will be added as a predictor for each channel. Should be a positive
#'   (unrestricted) integer.
#' @param type \[`integer(1)`]\cr Either
#'   `"fixed"` or `"varying"` which indicates whether the coefficients of the
#'   added lag terms should vary in time or not.
#' @return An object of class `lags`.
#' @srrstats {G2.3a} Uses match.arg
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' obs(y ~ -1 + varying(~x), family = "gaussian") +
#'   lags(type = "varying") + splines(df = 20)
#'
#' # A two-channel categorical model with time-invariant predictors
#' # here, lag terms are specified manually
#' obs(x ~ z + lag(x) + lag(y), family = "categorical") +
#'   obs(y ~ z + lag(x) + lag(y), family = "categorical")
#'
#' # The same categorical model as above, but with the lag terms
#' # added using 'lags'
#' obs(x ~ z, family = "categorical") +
#'   obs(y ~ z, family = "categorical") +
#'   lags(type = "fixed")
#'
lags <- function(k = 1L, type = c("fixed", "varying", "random")) {
  type <- onlyif(is.character(type), tolower(type))
  type <- try(match.arg(type, c("fixed", "varying", "random")), silent = TRUE)
  stopifnot_(
    !inherits(type, "try-error"),
    "Argument {.arg type} must be either
    {.val fixed}, {.val varying}, or {.val random}."
  )
  stopifnot_(
    checkmate::test_integerish(
      x = k,
      lower = 1L,
      any.missing = FALSE,
      min.len = 1L,
      unique = TRUE,
    ),
    "Argument {.arg k} must be an {.cls integer} vector with positive values."
  )
  structure(
    list(k = as.integer(k), type = type),
    class = "lags"
  )
}

#' Create a Lagged Version of a Vector
#'
#' @param x \[`vector()`]\cr A vector of values.
#' @param k \[`integer(1)`]\cr Number of positions to lag by.
#' @noRd
lag_ <- function(x, k = 1L) {
  lag_idx <- seq_len(length(x) - k)
  out <- x
  out[seq_len(k)] <- NA
  out[k + lag_idx] <- x[lag_idx]
  out
}

#' Create a Leading Version of a Vector
#'
#' @param x \[`vector()`]\cr A vector of values.
#' @param k \[`integer(1)`]\cr Number of positions to lead by.
#' @noRd
lead_ <- function(x, k = 1L) {
  rev(lag_(rev(x), k))
}

#' Adds Default Shift Values to Terms of the Form `lag(y)`
#'
#' @param x A `language` object.
#' @noRd
complete_lags <- function(x) {
  if (identical(length(x), 1L)) {
    return(x)
  }
  if (identical(deparse1(x[[1L]]), "lag")) {
    xlen <- length(x)
    if (identical(xlen, 2L)) {
      x <- str2lang(
        paste0("lag(", str_quote(deparse1(x[[2L]])), ", ", "1)")
      )
    } else if (identical(xlen, 3L)) {
      k <- verify_lag(x[[3L]], deparse1(x))
      x <- str2lang(
        paste0("lag(", str_quote(deparse1(x[[2L]])), ", ", k, ")")
      )
    } else {
      stop_(c(
        "Invalid lag definition {.code {deparse1(x)}}:",
        `x` = "Too many arguments supplied to {.fun lag}."
      ))
    }
  } else {
    for (i in seq_along(x)) {
      x[[i]] <- complete_lags(x[[i]])
    }
  }
  x
}

#' Extract Lagged Variables from a Language Object
#'
#' @param x A `language` object
#' @noRd
find_lags <- function(x) {
  if (!is.recursive(x)) {
    return(character(0L))
  }
  if (is.call(x)) {
    if (identical(as.character(x[[1L]]), "lag")) {
      return(deparse1(x))
    } else {
      ulapply(x[-1L], find_lags)
    }
  }
}

#' Extract the Order of Lagged Variables from a Language Object
#'
#' @param x A `language` object
#' @noRd
find_lag_orders <- function(x) {
  if (!is.recursive(x)) {
    return(
      data.frame(var = character(0L), order = integer(0L))
    )
  }
  if (is.call(x)) {
    if (identical(as.character(x[[1L]]), "lag")) {
      if (length(x) == 2L) {
        return(data.frame(var = deparse1(x[[2L]]), order = 1L))
      } else {
        return(data.frame(var = deparse1(x[[2L]]), order = x[[3L]]))
      }
    } else {
      rbindlist_(lapply(x[-1L], find_lag_orders))
    }
  }
}

#' Extract Non-lag Variables from a Language Object
#'
#' @param x A `language` object
#' @noRd
find_nonlags <- function(x) {
  if (!is.recursive(x)) {
    if (is.name(x)) {
      return(as.character(x))
    }
  }
  if (is.call(x)) {
    if (!identical(as.character(x[[1L]]), "lag")) {
      ulapply(x[-1L], find_nonlags)
    } else {
      character(0L)
    }
  }
}

#' Extract Lag Definitions
#'
#' Extract variables and shifts of lagged terms of the form `lag(var, k)`
#' and return them as a data frame for post processing.
#'
#' @param x \[`list`]\cr A list of `character` vectors.
#' @noRd
extract_lags <- function(x) {
  lag_terms <- unlist(x)
  lag_regex <- regexec(
    pattern = paste0(
      "(?<src>lag\\((?<var>[^\\+\\)\\,]+?)",
      "(?:,\\s*(?<k>[0-9]+)){0,1}\\))"
    ),
    text = lag_terms,
    perl = TRUE
  )
  lag_matches <- regmatches(lag_terms, lag_regex)
  if (length(lag_matches) > 0L) {
    lag_map <- do.call("rbind", args = lag_matches)
    lag_map <- as.data.frame(lag_map[, -1L, drop = FALSE])
    lag_map$var <- str_unquote(lag_map$var)
    lag_map$k <- as.integer(lag_map$k)
    lag_map$k[is.na(lag_map$k)] <- 1L
    lag_map$present <- TRUE
    lag_map <- unique(lag_map)
    lag_var <- sort(unique(lag_map$var))
    expanded <- vector(mode = "list", length = length(lag_var))
    for (i in seq_along(lag_var)) {
      v <- lag_var[i]
      tmp <- lag_map[lag_map$var == v, ]
      tmp <- tmp[order(tmp$k), ]
      full <- data.frame(
        src = "",
        var = v,
        k = seq.int(1L, max(tmp$k)),
        present = FALSE
      )
      expanded[[i]] <- full[full$k %in% tmp$k, ] <- tmp
    }
    lag_map <- rbindlist_(expanded)
  } else {
    data.frame(
      src = character(0L),
      var = character(0L),
      k = integer(0L),
      present = logical(0L)
    )
  }
}

#' Verify that `k` in `lag(y, k)` Represents a Valid Shift Value Expression
#'
#' @param k \[`language`]\cr The shift value definition.
#' @param lag_str \[`character(1)`]\cr The full lag term definition.
#' @noRd
verify_lag <- function(k, lag_str) {
  k_str <- deparse1(k)
  k_coerce <- try(eval(k), silent = TRUE)
  stopifnot_(
    !inherits(k_coerce, "try-error"),
    "Invalid shift value expression {.code {k_str}}."
  )
  k_coerce <- tryCatch(
    expr = as.integer(k_coerce),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  k_len <- length(k_coerce)
  stopifnot_(
    !is.null(k_coerce) && !identical(k_len, 0L) && !any(is.na(k_coerce)),
    "Unable to coerce shift value to {.cls integer} in {.code {lag_str}}."
  )
  stopifnot_(
    identical(k_len, 1L),
    c(
      "Shift value must be a single {.cls integer} in {.fun lag}:",
      `x` = "Multiple shift values were found in {.code {lag_str}}."
    )
  )
  stopifnot_(
    k_coerce > 0L,
    c(
      "Shift value must be positive in {.fun lag}:",
      `x` = "Nonpositive shift value was found in {.code {lag_str}}."
    )
  )
  k_coerce
}

#' Parse Lag and Lags Definitions of a `dynamiteformula` Object
#'
#' Also processes the checks on random and latent factor components and adds
#' them to dformulas.
#'
#' @return A `list` with the following components:
#'   * `all` A complete `dynamiteformula` for all channels of the model.
#'   * `det` A `dynamiteformula` for all deterministic channels.
#'   * `stoch` A `dynamiteformula` for all stochastic channels.
#'   * `lag_pred` A `dynamiteformula` for lagged predictors, meaning those
#'     variables that are not response variables of any channel
#'   * `lag_det` A `dynamiteformula` for lags of deterministic channels.
#'   * `lag_stoch` A `dynamiteformula` for lags of stochastic channels.
#'
#' @inheritParams parse_data
#' @noRd
parse_lags <- function(dformula, data, group_var, time_var, verbose) {
  channels_det <- which_deterministic(dformula)
  channels_stoch <- which_stochastic(dformula)
  resp_all <- get_responses(dformula)
  resp_stoch <- resp_all[channels_stoch]
  n_channels <- length(resp_all)
  for (i in seq_len(n_channels)) {
    fix_rhs <- complete_lags(formula_rhs(dformula[[i]]$formula))
    dformula[[i]]$formula <- as.formula(
      paste0(str_quote(resp_all[i]), "~", deparse1(fix_rhs))
    )
  }
  data_names <- names(data)
  non_lags <- unique(unlist(get_nonlag_terms(dformula[channels_stoch])))
  valid_resp <- c(resp_all, data_names)
  mis_vars <- which(!non_lags %in% valid_resp)
  stopifnot_(
    identical(length(mis_vars), 0L),
    "Can't find variable{?s} {.var {non_lags[mis_vars]}} in {.arg data}."
  )
  lag_map <- extract_lags(get_lag_terms(dformula))
  gl <- parse_global_lags(dformula, lag_map, resp_stoch, channels_stoch)
  dformula <- gl$dformula
  max_lag <- gl$max_lag
  mis_lags <- which(!lag_map$var %in% c(resp_all, data_names))
  stopifnot_(
    identical(length(mis_lags), 0L),
    c(
      "Unable to construct lagged values of
      {.var {cs(lag_map$var[mis_lags])}}:",
      `x` = "Can't find such variables in {.var data}."
    )
  )
  stoch_k <- lag_map[lag_map$var %in% resp_stoch, "k"]
  max_lag <- ifelse_(
    length(stoch_k) > 0L,
    max(max_lag, stoch_k),
    max_lag
  )
  sl <- parse_singleton_lags(
    dformula,
    data,
    group_var,
    lag_map,
    valid_resp,
    verbose
  )
  dformula <- sl$dformula
  dformula_det <- dformula[channels_det]
  dformula_lag_pred <- sl$channels[sl$pred]
  dformula_lag_stoch <- c(
    gl$channels[gl$stoch],
    sl$channels[sl$stoch & !sl$pred]
  )
  dformula_lag_det <- c(
    gl$channels[!gl$stoch],
    sl$channels[!sl$stoch & !sl$pred]
  )
  ro_lp <- sl$rank[sl$pred]
  cg <- attr(dformula, "channel_groups")
  attr(dformula_lag_pred, "rank_order") <- order(ro_lp)
  attr(dformula_lag_pred, "original_response") <- sl$resp[sl$pred]
  attr(dformula_lag_pred, "original_shift") <- ro_lp
  attr(dformula_lag_det, "rank_order") <-
    order(c(gl$rank[!gl$stoch], sl$rank[!sl$stoch & !sl$pred]))
  attr(dformula, "max_lag") <- max_lag
  list(
    all = dformula,
    det = dformula_det,
    stoch = structure(
      dformula[channels_stoch],
      channel_groups = rank_(cg[channels_stoch]),
      splines = attr(dformula, "splines"),
      random_spec = attr(dformula, "random_spec"),
      lfactor = attr(dformula, "lfactor")
    ),
    lag_pred = dformula_lag_pred,
    lag_det = dformula_lag_det,
    lag_stoch = dformula_lag_stoch
  )
}

#' Parse and Add Lags Defined via `lags` to a `dynamiteformula` Object
#'
#' @inheritParams parse_global_lags
#' @param channels_stoch \[`logical()`]\cr A vector indicating which channels
#'   are stochastic.
#' @param increment \[`logical()`]\cr  A vector indicating whether to add
#'   the new lag term or not (e.g.,, whether it was already present or not).
#' @param type \[`character(1)`]\cr
#'   Either `"fixed"`, `"varying"`, or `"random"`.
#' @param lhs \[`character()`]\cr A vector of the new lagged variable names.
#' @noRd
parse_new_lags <- function(dformula, channels_stoch, increment, type, lhs) {
  for (i in seq_along(channels_stoch)) {
    j <- channels_stoch[i]
    if (any(increment[[i]])) {
      dformula[[j]] <- formula_specials(
        x = increment_formula(
          formula = dformula[[j]]$formula,
          specials = dformula[[j]]$specials,
          x = lhs[increment[[j]]],
          type = type,
          varying_idx = dformula[[j]]$varying,
          fixed_idx = dformula[[j]]$fixed,
          random_idx = dformula[[j]]$random,
          varying_icpt = dformula[[j]]$has_varying_intercept,
          fixed_icpt = dformula[[j]]$has_fixed_intercept,
          random_icpt = dformula[[j]]$has_random_intercept
        ),
        original = dformula[[j]]$original,
        family = dformula[[j]]$family
      )
    }
  }
  dformula
}

#' Parse Lags That Actually Appear in a `dynamiteformula` Object
#'
#' @inheritParams parse_singleton_lags
#' @param y \[`list()`]\cr Output of `prepare_lagged_response`.
#' @param i \[`integer(1)`]\cr Row index of lag_map.
#' @param lhs \[`character(1)`]\cr Name of the new lagged response.
#' @noRd
parse_present_lags <- function(dformula, lag_map, y, i, lhs) {
  k <- y$lag_idx[i]
  if (lag_map$present[k]) {
    for (j in seq_along(dformula)) {
      dformula[[j]]$formula <- as.formula(
        gsub(
          pattern = lag_map$src[k],
          replacement = str_quote(lhs),
          x = deparse1(dformula[[j]]$formula),
          fixed = TRUE
        )
      )
    }
  }
  dformula
}

#' Parse a `lags` Definition in a `dynamiteformula` Object
#'
#' @inheritParams parse_lags
#' @param lag_map \[`data.frame`]\cr Output of `extract_lags`.
#' @param resp_stoch \[`character()`]\cr
#'   A vector of stochastic response variable names.
#' @param channels_stoch \[`logical()`]\cr
#'   A vector indicating which channels are stochastic.
#' @noRd
parse_global_lags <- function(dformula, lag_map, resp_stoch, channels_stoch) {
  lags_def <- attr(dformula, "lags")
  idx <- 0L
  k <- lags_def$k
  type <- lags_def$type
  max_lag <- ifelse_(is.null(lags_def), 0L, max(k))
  n_stoch <- length(resp_stoch)
  n_lag <- max_lag * n_stoch
  channels <- vector(mode = "list", length = n_lag)
  dterms <- lapply(dformula[channels_stoch], function(y) {
    attr(terms(y$formula), "term.labels")
  })
  stoch <- logical(n_lag)
  rank <- integer(n_lag)
  lhs <- character(n_lag)
  include <- logical(n_lag)
  increment <- replicate(n_stoch, logical(n_lag), simplify = FALSE)
  for (i in seq_len(max_lag)) {
    for (j in seq_len(n_stoch)) {
      y <- resp_stoch[j]
      idx <- idx + 1L
      stoch[idx] <- identical(i, 1L)
      lhs[idx] <- paste0(y, "_lag", i)
      rank[idx] <- i
      rhs <- ifelse_(stoch[idx], y, paste0(y, "_lag", i - 1L))
      new_term <- logical(n_stoch)
      lag_term <- paste0("lag(", y, ", ", i, ")")
      for (l in seq_len(n_stoch)) {
        new_term[l] <- !lag_term %in% dterms[[l]]
        increment[[l]][idx] <- (i %in% k) && new_term[l]
      }
      include[idx] <- any(new_term)
      channels[[idx]] <- dynamitechannel(
        formula = as.formula(paste0(lhs[idx], " ~ ", rhs)),
        family = deterministic_(),
        response = lhs[idx]
      )
    }
  }
  list(
    dformula = parse_new_lags(dformula, channels_stoch, increment, type, lhs),
    channels = channels[include],
    max_lag = max_lag,
    rank = rank[include],
    stoch = stoch[include]
  )
}

#' Parse Manual Lag Terms in a `dynamiteformula` Object
#'
#' @inheritParams parse_lags
#' @param lag_map \[`data.frame`]\cr Output of `extract_lags`.
#' @param valid_resp \[`character()`]\cr
#'   A vector of valid LHS variables that can  appear in the model formulas.
#' @noRd
parse_singleton_lags <- function(dformula, data, group_var,
                                 lag_map, valid_resp, verbose) {
  # avoid NSE notes from R CMD check
  group <- NULL
  n_lag <- nrow(lag_map)
  resp_all <- get_responses(dformula)
  channels <- vector(mode = "list", length = n_lag)
  resp_lag <- character(n_lag)
  pred <- logical(n_lag)
  stoch <- logical(n_lag)
  rank <- integer(n_lag)
  lag_var <- unique(lag_map$var)
  idx <- 0L
  for (resp in lag_var) {
    y <- prepare_lagged_response(
      dformula,
      lag_map,
      resp,
      resp_all,
      verbose
    )
    if (y$deterministic) {
      dformula[[y$idx]]$specials$past <- NULL
      dformula[[y$idx]]$specials$past_type <- NULL
    }
    for (i in seq_along(y$lag_idx)) {
      idx <- idx + 1L
      rhs <- ifelse_(identical(i, 1L), resp, paste0(resp, "_lag", i - 1L))
      lhs <- paste0(resp, "_lag", i)
      rank[idx] <- i
      stoch[idx] <- !y$deterministic
      pred[idx] <- !y$is_resp
      resp_lag[idx] <- resp
      spec <- NULL
      if (y$is_resp && !is.null(y$past_val)) {
        if (identical(y$past_type, "past")) {
          past_out <- lag_(y$past_val, i)
          na_idx <- data[,
            .I[base::seq_len(i)],
            by = group,
            env = list(i = i, group = group_var)
          ]$V1
          past_out[na_idx] <- NA
          spec <- list(
            past = past_out,
            resp_type = y$type
          )
        } else if (identical(y$past_type, "init")) {
          spec <- list(
            past = y$past_val[i],
            resp_type = y$type
          )
        }
      }
      channels[[idx]] <- dynamitechannel(
        formula = as.formula(paste0(str_quote(lhs), " ~ ", str_quote(rhs))),
        family = deterministic_(),
        response = lhs,
        specials = spec
      )
      dformula <- parse_present_lags(dformula, lag_map, y, i, lhs)
    }
  }
  list(
    dformula = dformula,
    channels = channels,
    pred = pred,
    rank = rank,
    resp = resp_lag,
    stoch = stoch
  )
}



#' Prepare a New Channel for a Lagged Response
#'
#' @inheritParams parse_singleton_lags
#' @param resp \[`character(1)`]\cr Name of the response variable being lagged
#' @param resp_all  \[`character()`]\cr Vector of all responses in the model.
#' @noRd
prepare_lagged_response <- function(dformula, lag_map,
                                    resp, resp_all, verbose) {
  y <- list()
  y$resp <- resp
  y$lag_idx <- which(lag_map$var == resp)
  y$src <- lag_map$src[y$lag_idx]
  y$idx <- which(resp_all == resp)
  y$is_resp <- length(y$idx) > 0L
  y$deterministic <- FALSE
  if (y$is_resp) {
    y$past_val <- dformula[[y$idx]]$specials$past
    y$past_idx <- NULL
    y$past_offset <- NULL
    y$deterministic <- is_deterministic(dformula[[y$idx]]$family)
    y$type <- dformula[[y$idx]]$specials$resp_type
    y$past_type <- dformula[[y$idx]]$specials$past_type
    if (y$deterministic && !identical(y$past_type, "past")) {
      y_max_lag <- max(lag_map$k[y$lag_idx])
      y_past_len <- length(y$past_val)
      if (y_past_len < y_max_lag) {
        if (verbose) {
          warning_(c(
            "Deterministic channel {.var {resp}} has a maximum lag of
            {y_max_lag} but you've supplied {cli::no(y_past_len)}
            initial value{?s//s}:",
            `i` = "This may result in NA values for {.var {resp}}."
          ))
        }
        y$past_val <- c(y$past_val, rep(NA, y_max_lag - y_past_len))
      }
    }
  }
  y
}

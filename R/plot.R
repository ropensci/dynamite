#' Traceplots and Density Plots for a `dynamitefit` Object
#'
#' Produces the traceplots and the density plots of the model parameters.
#'
#' @export
#' @family plotting
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param parameters \[`charecter()`]\ Parameter name(s) for which the plots
#'   should be drawn. Possible options can be found with the function
#'   [dynamite::get_parameter_names()]. The default is all parameters of a
#'   specific type for all responses, which can lead to too crowded a plot.
#' @param type \[`character(1)`]\cr Type of the parameter for which the plots
#'   should be drawn. Possible options can be found with the function
#'   [dynamite::get_parameter_types()]. Ignored if the argument `parameters`
#'   is supplied.
#' @param responses \[`character()`]\cr Response(s) for which the plots should
#'   be drawn. Possible options are `unique(x$priors$response)`. Default is
#'   all responses. Ignored if the argument `parameters` is supplied.
#' @param ... Not used..
#' @return A `ggplot` object.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.2, BS6.3, BS6.5} Implements the `plot`
#' method. Further plots can be easily constructed with the help of `as_draws`
#' combined with `ggplot2`, for example.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' plot(gaussian_example_fit, type = "beta")
#'
plot.dynamitefit <- function(x, parameters = NULL, type = NULL,
                             responses = NULL, ...) {
  stopifnot_(
    !is.null(parameters) || !is.null(type),
    "Both arguments {.arg parameters} and {.arg type} are missing, you should
    supply one of them."
  )
  stopifnot_(
    is.null(type) || checkmate::test_string(x = type, na.ok = FALSE),
    "Argument {.arg type} must be a single {.cls character} string."
  )
  out <- suppressWarnings(
    as_draws_df.dynamitefit(x, parameters = parameters,
      responses = responses, types = type)
  )
  # avoid NSE notes from R CMD check
  parameter <- density <- .iteration <- .chain <- NULL
  vars <- setdiff(names(out), c(".chain", ".iteration", ".draw"))
  p_list <- vector(mode = "list", length = 2L * length(vars))
  for (i in seq_along(vars)) {
    v <- vars[i]
    d <- out[, c(v, ".chain", ".iteration", ".draw")]
    p_list[[2L * (i - 1L) + 1L]] <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = !!rlang::sym(v), y = ggplot2::after_stat(density))
    ) +
      ggplot2::geom_density() +
      ggplot2::labs(x = v, y = "") +
      ggplot2::scale_x_continuous(expand = c(0.0, 0.0))
    p_list[[2L * i]] <- ggplot2::ggplot(
      d,
      ggplot2::aes(
        x = .iteration,
        y = !!rlang::sym(v),
        linetype = factor(.chain)
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::labs(x = "", y = v) +
      ggplot2::guides(linetype = ggplot2::guide_legend(title = "Chain"))
  }
  patchwork::wrap_plots(
    p_list,
    ncol = 2L,
    nrow = length(vars),
    byrow = TRUE
  )
}

#' Plot the Model Structure as a Directed Acyclic Graph (DAG)
#'
#' Plot a snapshot of the model structure at a specific time point with
#' a window of the highest-order lag dependency both into the past and the
#' future as a directed acyclic graph (DAG). Only response variables are
#' shown in the plot. This function can also produce a TikZ code of the DAG
#' to be used in reports and publications.
#'
#' @export
#' @family plotting
#' @param x \[`dynamiteformula`]\cr The model formula.
#' @param show_auxiliary \[`logical(1)`]\cr Should deterministic auxiliary
#'   responses be shown in the plot? If `FALSE`, the vertices corresponding
#'   to such responses will be projected out. The default is `TRUE`.
#' @param show_covariates \[`logical(1)`]\cr Should unmodeled covariates be
#'   shown in the plot? The defaults is `FALSE`.
#' @param tikz \[`logical(1)`]\cr Should the DAG be returned in TikZ format?
#'   The default is `FALSE` returning a `ggplot` object instead.
#' @param vertex_size \[`double(1)`]\cr The size (radius) of the vertex circles
#'   used in the `ggplot` DAG. (The vertical and horizontal distances between
#'   vertices in the grid are 1, for reference.)
#' @param label_size \[`double(1)`]\cr Font size (in points) to use for the
#'   vertex labels in the `ggplot` DAG.
#' @param ... Not used..
#' @return A `ggplot` object, or a `character` string if `tikz = TRUE`.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' multichanneL_formula <- obs(g ~ lag(g) + lag(logp), family = "gaussian") +
#'   obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
#'   obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
#'   aux(numeric(logp) ~ log(p + 1))
#' # A ggplot
#' plot(multichannel_formula)
#' # TikZ format
#' plot(multichannel_formula, tikz = TRUE)
#'
plot.dynamiteformula <- function(x, show_auxiliary = TRUE,
                                 show_covariates = FALSE, tikz = FALSE,
                                 vertex_size = 0.25, label_size = 18, ...) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamiteformula(x),
    "Argument {.arg x} must be a {.cls dynamiteformula} object."
  )
  stopifnot_(
    checkmate::test_flag(x = show_auxiliary),
    "Argument {.arg show_auxiliary} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_flag(x = show_covariates),
    "Argument {.arg show_covariates} must be a single {.cls logical} value."
  )
  stopifnot_(
    checkmate::test_flag(x = tikz),
    "Argument {.arg tikz} must be a single {.cls logical} value."
  )
  g <- get_dag(x, project = !show_auxiliary, covariates = show_covariates)
  ifelse_(
    tikz,
    plot_dynamiteformula_tikz(g),
    plot_dynamiteformula_ggplot(g, vertex_size, label_size)
  )
}

#' Create a TikZ plot of the Directed Acyclic Graph of a `dynamiteformula`
#'
#' @param g \[`list`]\cr Output of `get_dag`.
#' @noRd
plot_dynamiteformula_tikz <- function(g) {
  edgelist <- g$edgelist
  layout <- g$layout
  preamble <- paste_rows(
    "% Preamble",
    "\\usepackage{tikz}",
    "\\usetikzlibrary{positioning, arrows.meta, shapes.geometric}",
    "\\tikzset{%",
    "  semithick,",
    "  >={Stealth[width=1.5mm,length=2mm]},",
    "  obs/.style 2 args = {",
    "    name = #1, circle, draw, inner sep = 8pt, label = center:$#2$",
    "  }",
    "}",
    .parse = FALSE
  )
  n_v <- nrow(layout)
  layout$alias <- paste0("v", seq_len(n_v))
  vertices <- character(n_v)
  for (i in seq_along(vertices)) {
    vertices[i] <- paste0(
      "  \\node [obs = {",
      layout$alias[i], "}{",
      layout$var[i], "}] at (",
      layout$x[i], ", ",
      layout$y[i], ") {\\vphantom{0}};"
    )
  }
  n_e <- nrow(edgelist)
  edges <- character(n_e)
  for (i in seq_along(edges)) {
    from <- layout$alias[layout$var == edgelist$from[i]]
    to <- layout$alias[layout$var == edgelist$to[i]]
    x <- layout$x[layout$alias == from]
    y <- layout$y[layout$alias == from]
    xend <- layout$x[layout$alias == to]
    yend <- layout$y[layout$alias == to]
    vec <- c(xend - x, yend - y)
    curved <- (vec[1L] == 0 && abs(vec[2L]) > 1.0) ||
      (abs(vec[2L]) == 0 && vec[1L] > 1.0) ||
      (vec[1L] >= 2 && abs(vec[2L]) >= 2 && vec[1L] == abs(vec[2L]))
    edges[i] <- ifelse_(
      curved,
      paste0("  \\draw [->] (", from, ") to[bend right=45] (", to, ");"),
      paste0("  \\draw [->] (", from, ") -- (", to, ");")
    )
  }
  paste_rows(
    preamble,
    "% DAG",
    "\\begin{tikzpicture}",
    vertices,
    edges,
    "\\end{tikzpicture}",
    .parse = FALSE
  )
}

#' Create a `ggplot` of the Directed Acyclic Graph of a `dynamiteformula`
#'
#' @param g \[`list`]\cr Output of `get_dag`.
#' @noRd
plot_dynamiteformula_ggplot <- function(g, vertex_size, label_size) {
  v <- colnames(g$A)
  layout <- g$layout
  edgelist <- g$edgelist
  layout$var_expr <- gsub("(.+)_\\{(.+)\\}", "\\1\\[\\2\\]", layout$var)
  p <- ggplot2::ggplot(
    ggplot2::aes(x = x, y = y, label = var_expr),
    data = layout
  ) +
    ggplot2::theme_void()
  any_curved_x <- FALSE
  any_curved_y <- FALSE
  for (i in seq_len(nrow(edgelist))) {
    from <- edgelist$from[i]
    to <- edgelist$to[i]
    x <- layout$x[layout$var == from]
    y <- layout$y[layout$var == from]
    xend <- layout$x[layout$var == to]
    yend <- layout$y[layout$var == to]
    vec <- c(xend - x, yend - y)
    curved_x <- vec[1L] == 0 && abs(vec[2L]) > 1
    curved_y <- abs(vec[2L]) == 0 && vec[1L] > 1
    curved_diag <- vec[1L] >= 2 && abs(vec[2L]) >= 2 && vec[1L] == abs(vec[2L])
    any_curved_x <- any_curved_x || curved_x
    any_curved_y <- any_curved_y || curved_y
    arrow_fun <- ifelse_(
      curved_x || curved_y || curved_diag,
      ggplot2::geom_curve,
      ggplot2::geom_segment
    )
    rotation <- ifelse_(
      curved_x || curved_y,
      1/sqrt(2) * matrix(c(1, 1, -1, 1), 2, 2),
      diag(2L)
    )
    vec <- rotation %*% (vec / sqrt(sum(vec^2)))
    xend <- xend - vertex_size * vec[1L]
    yend <- yend - vertex_size * vec[2L]
    edge_data <- data.frame(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    )
    p <- p +
      arrow_fun(
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        data = edge_data,
        inherit.aes = FALSE,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.033, "npc"))
      )
  }
  p + ggforce::geom_circle(
    ggplot2::aes(x0 = x, y0 = y, r = vertex_size),
    fill = "white"
  ) +
    ggplot2::geom_text(
      vjust = 0.4,
      parse = TRUE,
      size = label_size / ggplot2::.pt
    ) +
    ggplot2::coord_equal(
      xlim = c(
        min(layout$x) - 0.15 - 0.25 * any_curved_x,
        max(layout$x) + 0.15 + 0.25 * any_curved_x
      ),
      ylim = c(
        min(layout$y) - 0.15 - 0.25 * any_curved_y,
        max(layout$y) + 0.15 + 0.25 * any_curved_y
      )
    )
}

#' Plot Time-varying Regression Coefficients of a Dynamite Model
#'
#' @export
#' @family plotting
#' @param x \[`dynamitefit`]\cr The model fit object
#' @param parameters \[`charecter()`]\ Parameter name(s) for which the plots
#'   should be drawn. Possible options can be found with function
#'   `get_parameter_names(types = "delta")`.
#' @param responses \[`character()`]\cr Response(s) for which the coefficients
#'   should be drawn. Possible options are elements of
#'   `unique(x$priors$response)`, and the default is this whole vector.
#    Ignored if the argument `parameters` is supplied.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity level for `geom_ribbon`.
#'   Default is 0.5.
#' @param scales \[`character(1)`] Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), plots also
#'   the time-varying alphas if such parameters exists in the model.
#' @return A `ggplot` object.
#' @srrstats {G2.3a} Uses match.arg.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' plot_deltas(gaussian_example_fit, level = 0.025, scales = "free") +
#'   ggplot2::theme_minimal()
#'
plot_deltas <- function(x, parameters = NULL, responses = NULL, level = 0.05,
                        alpha = 0.5, scales = c("fixed", "free"),
                        include_alpha = TRUE) {
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
    checkmate::test_number(
      x = level,
      lower = 0.0,
      upper = 1.0,
      na.ok = FALSE,
    ),
    "Argument {.arg level} must be a single
     {.cls numeric} value between 0 and 1."
  )
  stopifnot_(
    checkmate::test_number(
      x = alpha,
      lower = 0.0,
      upper = 1.0,
      na.ok = FALSE,
    ),
    "Argument {.arg alpha} must be a single
     {.cls numeric} value between 0 and 1."
  )
  if (!is.null(parameters)) {
    delta_names <- get_parameter_names(x, types = "delta")
    found_pars <- parameters %in% delta_names
    stopifnot_(all(found_pars),
      c(
        "Parameter{?s} {.var {parameters[!found_pars]}} not found or
        {?it is/they are} of wrong type:",
        `i` = 'Use {.fun get_parameter_names} with {.arg types = "delta"} to
        check suitable parameter names.'
      )
    )
  }
  coefs <- coef.dynamitefit(
    x,
    parameters = parameters,
    type = "delta",
    responses = responses,
    probs = c(level, 1 - level),
    include_alpha = include_alpha
  )
  stopifnot_(
    nrow(coefs) > 0L,
    "The model does not contain varying coefficients delta."
  )
  scales <- onlyif(is.character(scales), tolower(scales))
  scales <- try(match.arg(scales, c("fixed", "free")), silent = TRUE)
  stopifnot_(
    !inherits(scales, "try-error"),
    "Argument {.arg scales} must be either {.val fixed} or {.val free}."
  )
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the time-varying coefficients"
  )
  # avoid NSE notes from R CMD check
  time <- mean <- category <- parameter <- NULL
  if (any(!is.na(coefs$category))) {
    p <- ggplot2::ggplot(
      coefs,
      ggplot2::aes(
        time,
        mean,
        colour = category,
        fill = category
      )
    )
  } else {
    p <- ggplot2::ggplot(coefs, ggplot2::aes(time, mean))
  }
  p +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = !!rlang::sym(paste0("q", 100 * level)),
        ymax = !!rlang::sym(paste0("q", 100 * (1 - level)))
      ),
      alpha = alpha,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::scale_x_continuous(limits = range(coefs$time)) +
    ggplot2::facet_wrap("parameter", scales = scales) +
    ggplot2::labs(title = title, x = "Time", y = "Value")
}

#' Plot Time-invariant Regression Coefficients of a Dynamite Model
#'
#' @export
#' @family plotting
#' @inheritParams plot_deltas
#' @param parameters \[`charecter()`]\ Parameter name(s) for which the plots
#'   should be drawn. Possible options can be found with function
#'   `get_parameter_names(types = "beta")`.
#' @param include_alpha \[`logical(1)`]\cr If `TRUE` (default), plots also
#'   the time-invariant alphas if such parameters exists in the model.
#' @return A `ggplot` object.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' plot_betas(gaussian_example_fit, level = 0.1)
#'
plot_betas <- function(x, parameters = NULL, responses = NULL, level = 0.05,
                       include_alpha = TRUE) {
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
    checkmate::test_number(
      x = level,
      lower = 0.0,
      upper = 1.0,
      na.ok = FALSE,
    ),
    "Argument {.arg level} must be a single
     {.cls numeric} value between 0 and 1."
  )
  if (!is.null(parameters)) {
    beta_names <- get_parameter_names(x, types = "beta")
    found_pars <- parameters %in% beta_names
    stopifnot_(all(found_pars),
      c(
        "Parameter{?s} {.var {parameters[!found_pars]}} not found or
        {?it is/they are} of wrong type:",
        `i` = 'Use {.fun get_parameter_names} with {.arg types = "beta"} to
        check suitable parameter names.'
      )
    )
  }
  coefs <- coef.dynamitefit(
    x,
    parameters = parameters,
    type = "beta",
    responses = responses,
    probs = c(level, 1 - level),
    include_alpha = include_alpha
  )
  stopifnot_(
    nrow(coefs) > 0L,
    "The model does not contain fixed coefficients beta."
  )
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the time-invariant coefficients"
  )
  # avoid NSE notes from R CMD check
  time <- mean <- category <- parameter <- NULL
  if (any(!is.na(coefs$category))) {
    p <- ggplot2::ggplot(
      coefs,
      ggplot2::aes(
        mean,
        parameter,
        colour = category,
        group = category
      )
    )
  } else {
    p <- ggplot2::ggplot(coefs, ggplot2::aes(mean, parameter))
  }
  p +
    ggplot2::geom_pointrange(
      ggplot2::aes(
        xmin = !!rlang::sym(paste0("q", 100 * level)),
        xmax = !!rlang::sym(paste0("q", 100 * (1 - level)))
      ),
      position = ggplot2::position_dodge(0.5)
    ) +
    ggplot2::labs(title = title, x = "Value", y = "Parameter")
}

#' Plot Random effects of a Dynamite Model
#'
#' Note that as this function tries to draw a plot containing effects of all
#' groups, the plot will become messy with large number of groups.
#'
#' @export
#' @family plotting
#' @inheritParams plot_deltas
#' @param groups Group name(s) for which the plots should be drawn.
#'   Default is all groups.
#' @return A `ggplot` object.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' plot_nus(gaussian_example_fit)
#'
plot_nus <- function(x, parameters = NULL, responses = NULL, level = 0.05,
                     groups = NULL) {
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
    checkmate::test_number(
      x = level,
      lower = 0.0,
      upper = 1.0,
      na.ok = FALSE,
    ),
    "Argument {.arg level} must be a single
     {.cls numeric} value between 0 and 1."
  )
  if (!is.null(parameters)) {
    nu_names <- get_parameter_names(x, types = "nu")
    found_pars <- parameters %in% nu_names
    stopifnot_(all(found_pars),
      c(
        "Parameter{?s} {.var {parameters[!found_pars]}} not found or
        {?it is/they are} of wrong type:",
        `i` = 'Use {.fun get_parameter_names} with {.arg types = "nu"} to
        check suitable parameter names.'
      )
    )
  }
  coefs <- try(
    coef.dynamitefit(
      x,
      parameters = parameters,
      type = "nu",
      responses = responses,
      probs = c(level, 1 - level)
    ),
    silent = TRUE
  )
  stopifnot_(
    !inherits(coefs, "try-error"),
    "The model does not contain random effects nu."
  )
  coefs <- ifelse_(
    is.null(groups),
    coefs,
    coefs[coefs$group %in% groups, ]
  )
  # avoid NSE notes from R CMD check
  mean <- parameter <- NULL
  coefs$parameter <- glue::glue("{coefs$parameter}_{coefs$category}_{coefs$group}")
  # remove NAs from parameters which are not category specific
  coefs$parameter <- gsub("_NA", "", coefs$parameter)
  coefs$parameter <- factor(coefs$parameter, levels = coefs$parameter)
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the random intercepts"
  )
  ggplot2::ggplot(coefs, ggplot2::aes(mean, parameter)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(
        xmin = !!rlang::sym(paste0("q", 100 * level)),
        xmax = !!rlang::sym(paste0("q", 100 * (1 - level)))
      )
    ) +
    ggplot2::labs(title = title, x = "Value", y = "Parameter")
}
#' Plot Factor Loadings of a Dynamite Model
#'
#' @export
#' @family plotting
#' @inheritParams plot_deltas
#' @return A `ggplot` object.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#'
plot_lambdas <- function(x, responses = NULL, level = 0.05) {
  stopifnot_(
    checkmate::test_number(
      x = level,
      lower = 0.0,
      upper = 1.0,
      na.ok = FALSE,
    ),
    "Argument {.arg level} must be a single
     {.cls numeric} value between 0 and 1."
  )
  coefs <- try(
    coef.dynamitefit(
      x,
      type = "lambda",
      responses = responses,
      probs = c(level, 1 - level)
    ),
    silent = TRUE
  )
  stopifnot_(
    !inherits(coefs, "try-error"),
    "The model does not contain latent factor psi."
  )
  # avoid NSE notes from R CMD check
  time <- mean <- parameter <- NULL
  coefs$parameter <- glue::glue("{coefs$parameter}_{coefs$group}")
  coefs$parameter <- factor(coefs$parameter, levels = coefs$parameter)
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the latent factor loadings"
  )
  ggplot2::ggplot(coefs, ggplot2::aes(mean, parameter)) +
    ggplot2::geom_pointrange(ggplot2::aes(
      xmin = !!rlang::sym(paste0("q", 100 * level)),
      xmax = !!rlang::sym(paste0("q", 100 * (1 - level)))
    )) +
    ggplot2::labs(title = title, x = "Value", y = "Parameter")
}
#' Plot Latent Factors of a Dynamite Model
#'
#' @export
#' @family plotting
#' @param x \[`dynamitefit`]\cr The model fit object
#' @param responses  \[`character()`]\cr Response(s) for which the coefficients
#'   should be drawn. Possible options are elements of
#'   `unique(x$priors$response)`, and the default is this whole vector.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity level for `geom_ribbon`.
#'   Default is 0.5.
#' @param scales \[`character(1)`] Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @return A `ggplot` object.
#' @srrstats {G2.3a} Uses match.arg.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.3} Implements the `plot` method.
#'
plot_psis <- function(x, responses = NULL, level = 0.05, alpha = 0.5,
                      scales = c("fixed", "free")) {
  stopifnot_(
    checkmate::test_number(
      x = level,
      lower = 0.0,
      upper = 1.0,
      na.ok = FALSE,
    ),
    "Argument {.arg level} must be a single
     {.cls numeric} value between 0 and 1."
  )
  stopifnot_(
    checkmate::test_number(
      x = alpha,
      lower = 0.0,
      upper = 1.0,
      na.ok = FALSE,
    ),
    "Argument {.arg alpha} must be a single
     {.cls numeric} value between 0 and 1."
  )
  coefs <- coef.dynamitefit(
    x,
    type = "psi",
    responses = responses,
    probs = c(level, 1 - level)
  )
  stopifnot_(
    nrow(coefs) > 0L,
    "The model does not contain latent factor psi."
  )
  scales <- onlyif(is.character(scales), tolower(scales))
  scales <- try(match.arg(scales, c("fixed", "free")), silent = TRUE)
  stopifnot_(
    !inherits(scales, "try-error"),
    "Argument {.arg scales} must be either {.val fixed} or {.val free}."
  )
  title <- paste0(
    "Posterior mean and ",
    100 * (1 - 2 * level),
    "% intervals of the latent factors"
  )
  # avoid NSE notes from R CMD check
  time <- mean <- category <- NULL
  if (any(!is.na(coefs$category))) {
    p <- ggplot2::ggplot(
      coefs,
      ggplot2::aes(
        time,
        mean,
        colour = category,
        fill = category
      )
    )
  } else {
    p <- ggplot2::ggplot(coefs, ggplot2::aes(time, mean))
  }
  p +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = !!rlang::sym(paste0("q", 100 * level)),
        ymax = !!rlang::sym(paste0("q", 100 * (1 - level)))
      ),
      alpha = alpha
    ) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap("parameter", scales = scales) +
    ggplot2::labs(title = title, x = "Time", y = "Value")
}

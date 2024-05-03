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
#' multichannel_formula <- obs(g ~ lag(g) + lag(logp), family = "gaussian") +
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
  if (tikz) {
    g <- get_dag(
      x,
      project = !show_auxiliary,
      covariates = show_covariates,
      format = "default"
    )
    plot_dynamiteformula_tikz(g)
  } else {
    g <- get_dag(
      x,
      project = !show_auxiliary,
      covariates = show_covariates,
      format = "expression"
    )
    plot_dynamiteformula_ggplot(g, vertex_size, label_size)
  }
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
  # avoid NSE notes from R CMD check
  var_expr <- NULL
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
    curved <- curved_x || curved_y || curved_diag
    curvature <- ifelse_(curved_diag, 0.5, sqrt(0.5)^max(abs(vec)))
    angle <- pi / max(abs(vec) + 2)
    rotation <- ifelse_(
      curved,
      matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2),
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
      ifelse_(
        curved,
        ggplot2::geom_curve(
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          data = edge_data,
          inherit.aes = FALSE,
          arrow = ggplot2::arrow(length = ggplot2::unit(0.033, "npc")),
          curvature = curvature
        ),
        ggplot2::geom_segment(
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          data = edge_data,
          inherit.aes = FALSE,
          arrow = ggplot2::arrow(length = ggplot2::unit(0.033, "npc"))
        )
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

#' Plots for `dynamitefit` Objects
#'
#' Produces the traceplots and the density plots of the model parameters. Can
#' also be used to plot the time-varying and time-invariant parameters of the
#' model along with their posterior intervals. See the `plot_type` argument
#' for details on available plots.
#'
#' @export
#' @family plotting
#' @param x \[`dynamitefit`]\cr The model fit object.
#' @param plot_type \[`character(1)`]\cr What type of plot to draw? The default
#'   is `"default"` which draws posterior means and intervals of the parameters
#'   selected by `types` or `parameters`. If both `"types"` and
#'   `parameters` are `NULL`, all parameters are drawn up to the maximum
#'   specified by `n_params`. Option `"trace"` instead draws posterior
#'   densities and traceplots of the parameters.
#' @param types \[`character(1)`]\cr Types of the parameter for which the plots
#'   should be drawn. Possible options can be found with the function
#'   [dynamite::get_parameter_types()]. Ignored if the argument `parameters`
#'   is supplied.
#' @param parameters \[`charecter()`]\cr Parameter name(s) for which the plots
#'   should be drawn. Possible options can be found with the function
#'   [dynamite::get_parameter_names()]. The default is all parameters,
#'   limited by `n_params`.
#' @param responses \[`character()`]\cr Response(s) for which the plots should
#'   be drawn. Possible options are `unique(x$priors$response)`. Default is
#'   all responses. Ignored if the argument `parameters` is supplied.
#' @param times \[`double()`]\cr Time point(s) for which the plots should be
#'   drawn for time-varying parameters. By default, all time points are
#'   included, up to the maximum number of parameters specified by `n_params`
#'   starting from the first non-fixed time point.
#' @param groups \[`character(1)`]\cr Group name(s) for which the plots
#'   should be drawn for group-specific parameters.
#' @param level \[`numeric(1)`]\cr Level for posterior intervals.
#'   Default is 0.05, leading to 90% intervals.
#' @param alpha \[`numeric(1)`]\cr Opacity level for `geom_ribbon`.
#'   Default is 0.5.
#' @param scales \[`character(1)`]\cr Should y-axis of the panels be `"fixed"`
#'   (the default) or `"free"`? See [ggplot2::facet_wrap()].
#' @param n_params \[`integer()`]\cr A single value or a vector of length 2
#'   specifying the maximum number of parameters to plot. If a single value
#'   is provided, the same limit is used for all parameters. If a vector is
#'   supplied, the first element defines the maximum number of time-invariant
#'   parameters to plot and the second the maximum number of time-varying
#'   parameters to plot. The defaults values are 50 for time-invariant
#'   parameters and 3 for time-varying parameters. The default value is 5
#'   for `plot_type == "trace"`.
#' @param ... Not used..
#' @return A `ggplot` object.
#' @srrstats {BS6.1, RE6.0, RE6.1, BS6.2, BS6.3, BS6.5} Implements the `plot()`
#' method. Further plots can be easily constructed with the help of
#' `as_draws()` combined with `ggplot2`, for example.
#' @examples
#' data.table::setDTthreads(1) # For CRAN
#' plot(gaussian_example_fit, type = "beta")
#'
plot.dynamitefit <- function(x, plot_type = c("default", "trace"),
                             types = NULL, parameters = NULL,
                             responses = NULL, groups = NULL, times = NULL,
                             level = 0.05, alpha = 0.5,
                             scales = c("fixed", "free"),
                             n_params = NULL, ...) {
  stopifnot_(
    !missing(x),
    "Argument {.arg x} is missing."
  )
  stopifnot_(
    is.dynamitefit(x),
    "Argument {.arg x} must be a {.cls dynamitefit} object."
  )
  plot_type <- onlyif(is.character(plot_type), tolower(plot_type))
  plot_type <- try(
    match.arg(plot_type, c("default", "trace")),
    silent = TRUE
  )
  stopifnot_(
    !inherits(plot_type, "try-error"),
    "Argument {.arg type} must be either {.val default} or {.val trace}."
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
  scales <- onlyif(is.character(scales), tolower(scales))
  scales <- try(match.arg(scales, c("fixed", "free")), silent = TRUE)
  stopifnot_(
    !inherits(scales, "try-error"),
    "Argument {.arg scales} must be either {.val fixed} or {.val free}."
  )
  stopifnot_(
    checkmate::test_integerish(
      x = n_params,
      lower = 1.0,
      any.missing = FALSE,
      null.ok = TRUE,
      min.len = 1L,
      max.len = 2L
    ),
    "Argument {.arg n_params} must a
     positive integer vector of at most length 2."
  )
  n_params <- ifelse_(
    !is.null(n_params) && length(n_params) == 1L,
    c(n_params, n_params),
    n_params
  )
  if (!is.null(parameters)) {
    responses <- types <- NULL
  }
  if (identical(plot_type, "trace")) {
    plot_trace(
      x,
      types,
      parameters,
      responses,
      times,
      groups,
      n_params
    )
  } else {
    coefs <- coef.dynamitefit(
      x,
      types = types,
      parameters = parameters,
      responses = responses,
      times = times,
      groups = groups,
      probs = c(level, 1 - level)
    )
    p_fixed <- plot_fixed(
      coefs,
      level,
      alpha,
      scales,
      n_params[1L]
    )
    p_varying <- plot_varying(
      coefs,
      level,
      alpha,
      scales,
      n_params[2L]
    )
    stopifnot_(
      !is.null(p_fixed) || !is.null(p_varying),
      "The model does not contain any of the specified parameters."
    )
    if (is.null(p_fixed)) {
      p_varying
    } else if (is.null(p_varying)) {
      p_fixed
    } else {
      patchwork::wrap_plots(
        p_fixed,
        patchwork::free(p_varying),
        heights = c(0.33, 0.66),
        ncol = 1L,
        nrow = 2L,
        byrow = TRUE,
        guides = "keep"
      )
    }
  }
}

#' Traceplots and Density Plots for a `dynamitefit` Object
#'
#' @inheritParams plot.dynamitefit
#' @noRd
plot_trace <- function(x, types, parameters, responses,
                       times, groups, n_params) {
  out <- suppressWarnings(
    as_draws_df.dynamitefit(
      x,
      parameters = parameters,
      responses = responses,
      types = types,
      times = times,
      groups = groups
    )
  )
  # avoid NSE notes from R CMD check
  parameter <- density <- .iteration <- .chain <- NULL
  vars <- filter_params(out, n_params, 5L)
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

#' Plot Time-invariant Parameters of a Dynamite Model
#'
#' @inheritParams plot.dynamitefit
#' @noRd
plot_fixed <- function(coefs, level, alpha, scales, n_params) {
  coefs <- coefs[coefs$type %in% intersect(fixed_types, default_types), ]
  coefs <- coefs[is.na(coefs$time), ]
  if (nrow(coefs) == 0L) {
    return(NULL)
  }
  coefs <- filter_params(coefs, n_params, 50)
  n_coefs <- nrow(coefs)
  coefs$parameter <- glue::glue(
    "{coefs$parameter}_{coefs$category}_{coefs$group}"
  )
  # remove NAs from parameters which are not category specific
  coefs$parameter <- gsub("_NA", "", coefs$parameter)
  coefs$parameter <- factor(coefs$parameter, levels = coefs$parameter)
  title_spec <- "time-invariant parameters"
  if (n_unique(coefs$type) == 1L) {
    title_spec <- switch(
      coefs$type[1L],
      alpha = "time-invariant intercepts",
      beta = "time-invariant regression coefficients",
      cutpoints = "time-invariant cutpoints",
      nu = "random intercepts",
      lambda = "latent factor loadings",
      "time-invariant parameters"
    )
  }
  title <- glue::glue(
    "Posterior means and {100 * (1 - 2 * level)} ",
    "% intervals of the {title_spec}"
  )
  # avoid NSE notes from R CMD check
  time <- mean <- category <- parameter <- type <- NULL
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
    ggplot2::facet_wrap(~type, scales = "free") +
    ggplot2::geom_pointrange(
      ggplot2::aes(
        xmin = !!rlang::sym(paste0("q", 100 * level)),
        xmax = !!rlang::sym(paste0("q", 100 * (1 - level)))
      ),
      position = ggplot2::position_dodge(0.5)
    ) +
    ggplot2::labs(x = "Value", y = "Parameter") +
    ggplot2::ggtitle(title)
}


#' Plot Time-varying Parameters of a Dynamite Model
#'
#' @inheritParams plot.dynamitefit
#' @noRd
plot_varying <- function(coefs, level, alpha, scales, n_params) {
  coefs <- coefs[coefs$type %in% intersect(varying_types, default_types), ]
  coefs <- coefs[!is.na(coefs$time), ]
  if (nrow(coefs) == 0L) {
    return(NULL)
  }
  coefs <- filter_params(coefs, n_params, 3)
  n_coefs <- nrow(coefs)
  title_spec <- "time-varying parameters"
  if (n_unique(coefs$type) == 1L) {
    title_spec <- switch(
      coefs$type[1L],
      alpha = "time-varying intercepts",
      cutpoints = "time-invariant cutpoints",
      delta = "time-invariant regression coefficients",
      psi = "latent factors",
      "time-varying parameters"
    )
  }
  title <- glue::glue(
    "Posterior means and {100 * (1 - 2 * level)} %",
    "intervals of the {title_spec}"
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
    ggplot2::labs(x = "Time", y = "Value") +
    ggplot2::ggtitle(title)
}

#' Select only specific number of parameters for plotting
#'
#' @param x Output of `coef.dynamitefit()` or a `draws` object.
#' @param n_params Number of parameters to keep.
#' @param n_params_default Number of parameters to keep if `n_params` is `NULL`.
#' @noRd
filter_params <- function(x, n_params, n_params_default) {
  vars <- character(0L)
  keep_params <- logical(0L)
  is_draws <- inherits(x, "draws")
  n_params_set <- !is.null(n_params)
  n_params <- ifelse_(
    n_params_set,
    n_params,
    n_params_default
  )
  if (is_draws) {
    vars <- setdiff(names(x), c(".chain", ".iteration", ".draw"))
    n_vars <- length(vars)
    keep_params <- logical(n_vars)
    keep_params[seq_len(min(n_vars, n_params))] <- TRUE
    onlyif(
      !n_params_set && !all(keep_params),
      warning_(c(
        "Number of parameters to be plotted ({n_u_params}) exceeds the
       maximum number of parameters ({n_params}). The remaining parameters will
       not be plotted.",
       `i` = "Please increase {.arg n_params} to plot more parameters."
      ))
    )
  } else {
    params <- glue::glue("{x$parameter}_{x$group}")
    params <- gsub("_NA", "", params)
    keep_params <- logical(nrow(x))
    for (type in unique(x$type)) {
      type_params <- params[x$type == type]
      u_params <- unique(type_params)
      n_u_params <- length(u_params)
      max_params <- min(n_params, n_u_params)
      keep_type_params <- type_params %in% u_params[seq_len(n_params)]
      keep_params <- keep_params |
        params %in% u_params[seq_len(n_params)]
      onlyif(
        !n_params_set && !all(keep_type_params),
        warning_(c(
          "Number of parameters to be plotted ({n_u_params}) exceeds the
         maximum number of parameters ({max_params}) for parameters
         of type {.var type}. The remaining parameters of this type will
         not be plotted.",
         `i` = "Please increase {.arg n_params} to plot more parameters."
        ))
      )
    }
  }
  ifelse_(
    is_draws,
    vars[keep_params],
    x[keep_params, ]
  )
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.dynamitefit}
\alias{predict.dynamitefit}
\title{Predict Method for a \pkg{dynamite} Model}
\usage{
\method{predict}{dynamitefit}(
  object,
  newdata = NULL,
  type = c("response", "mean", "link"),
  funs = list(),
  impute = c("none", "locf", "nocb"),
  new_levels = c("none", "bootstrap", "gaussian", "original"),
  global_fixed = FALSE,
  n_draws = NULL,
  thin = 1,
  expand = TRUE,
  df = TRUE,
  drop = TRUE,
  ...
)
}
\arguments{
\item{object}{[\code{dynamitefit}]\cr The model fit object.}

\item{newdata}{[\code{data.frame}]\cr Data used in predictions. Predictions are
computed for missing (\code{NA}) values in the response variable columns, and
non-missing values are assumed fixed.
If \code{NULL} (default), the data used in model estimation is used for
predictions as well, after all values in the response variable columns
after the first \code{fixed} time point are converted to \code{NA} values.
Missing values in predictor columns can be imputed (argument \code{impute}).
There should be no new time points that were not present in the data that
were used to fit the model. New group levels can be included, but if the
model contains random effects, an option for the random
effects for the new levels must be chosen (argument \code{new_levels}).
If the grouping variable of the original data is missing, it is assumed
that all observations in \code{newdata} belong to the first group in the
original data. New group levels are not allowed for models using latent
factors.}

\item{type}{[\code{character(1)}]\cr Type of prediction,
\code{"response"} (default), \code{"mean"}, or \code{"link"}.}

\item{funs}{[\code{list()}]\cr A named list whose names should correspond to the
response variables of the model. Each element of \code{funs} should be a
a named \code{list} of functions that will be applied to the
corresponding predicted \code{type} of the channel over the individuals
for each combination of the posterior draws and time points.
In other words, the resulting predictions will be averages
over the individuals. The functions should take the corresponding
\code{type} variable values as their only argument.
If \code{funs} is empty, the full individual level values are returned
instead. Note that this argument can only be used
if there are multiple individuals (i.e., \code{group} was not \code{NULL} in the
\code{dynamite} call).}

\item{impute}{[\code{character(1)}]\cr Which imputation scheme to use for
missing exogenous predictor values. Currently supported options are
no imputation: \code{"none"} (default), last observation carried forward:
\code{"locf"}, and next observation carried backward: \code{"nocb"}.}

\item{new_levels}{[\code{character(1)}]\cr
Defines if and how to sample the random effects for observations whose
group level was not present in the original data. The options are:
\itemize{
\item \code{"none"} (the default) which will signal an error if new levels
are encountered.
\item \code{"bootstrap"} which will randomly draw from the posterior samples of
the random effects across all original levels.
\item \code{"gaussian"} which will randomly draw from a Gaussian
distribution using the posterior samples of the random effects
standard deviation (and correlation matrix if applicable).
\item \code{"original"} which will randomly match each new level to one of
the original levels. The posterior samples of the random effects of
the matched levels will then be used for the new levels.
}

This argument is ignored if the model does not contain random effects.}

\item{global_fixed}{[\code{logical(1)}]\cr If \code{FALSE} (the default),
the first non-fixed time point is counted from the the first non-NA
observation for each group member separately. Otherwise, the first
non-fixed time point is counted from the first time point globally.
If there are no groups, then the options are equivalent.}

\item{n_draws}{[\code{integer(1)}]\cr Number of posterior samples to use,
default is \code{NULL} which uses all samples without permuting (with chains
concatenated). If \code{n_draws}is smaller than \code{ndraws(object)}, a random
subset of \code{n_draws} posterior samples are used.}

\item{thin}{[\code{integer(1)}]\cr Use only every \code{thin} posterior sample.
This can be beneficial with when the model object contains
large number of samples. Default is \code{1} meaning that all samples are used.}

\item{expand}{[\code{logical(1)}]\cr If \code{TRUE} (the default), the output
is a single \code{data.frame} containing the original \code{newdata} and the
predicted values. Otherwise, a \code{list} is returned with two components,
\code{simulated} and \code{observed}, where the first contains only the
predicted values, and the second contains the original \code{newdata}.
Setting \code{expand} to \code{FALSE} can help conserve memory because \code{newdata}
is not replicated \code{n_draws} times in the output.
This argument is ignored if \code{funs} are provided.}

\item{df}{[\code{logical(1)}]\cr If \code{TRUE} (default) the output
consists of \code{data.frame} objects, and \code{data.table} objects otherwise.}

\item{drop}{[\code{logical(1)}]\cr If \code{TRUE} (default), the columns of \code{newdata}
that are not used by any model formula are dropped from the output.
If \code{FALSE}, all columns are kept.}

\item{...}{Ignored.}
}
\value{
A \code{data.frame} containing the predicted values or a \code{list} of two
\code{data.frames}. See the \code{expand} argument for details. Note that the
\code{.draw} column is not the same as \code{.draw} from \code{as.data.frame} and
\code{as_draws} methods as \code{predict} uses permuted samples. A mapping between
these variables can be done using information in
\code{object$stanfit@sim$permutation}.
}
\description{
Obtain counterfactual predictions for a \code{dynamitefit} object.
}
\details{
Note that forecasting (i.e., predictions for time indices beyond the last
time index in the original data) is not supported by the \pkg{dynamite}
package. However, such predictions can be obtained by augmenting the
original data with \code{NA} values before model estimation.
}
\examples{
data.table::setDTthreads(1) # For CRAN
out <- predict(gaussian_example_fit, type = "response", n_draws = 2L)
head(out)

# using summary functions
sumr <- predict(multichannel_example_fit, type = "mean",
  funs = list(g = list(m = mean, s = sd), b = list(sum = sum)),
  n_draws = 2L)
head(sumr$simulated)
\donttest{
# Please update your rstan and StanHeaders installation before running
# on Windows
if (!identical(.Platform$OS.type, "windows")) {
  # Simulate from the prior predictive distribution

  f <- obs(y ~ lag(y) + varying(~ -1 + x), "gaussian") +
    splines(df = 10, noncentered = TRUE)

  # Create data with missing observations
  # Note that due to the lagged term in the model,
  # we need to fix the first time point
  d <- data.frame(y = c(0, rep(NA, 49)), x = rnorm(50), time = 1:50)

  # Suppress warnings due to the lack of data
  suppressWarnings(
    priors <- get_priors(f, data = d, time = "time")
  )

  # Modify default priors which can produce exploding behavior when used
  # without data
  priors$prior <- c(
    "normal(0, 1)",
    "normal(0.6, 0.1)",
    "normal(-0.2, 0.5)",
    "normal(0.2, 0.1)",
    "normal(0.5, 0.1)"
  )

  # Samples from the prior conditional on the first time point and x
  fit <- dynamite(
    dformula = f,
    data = d,
    time = "time",
    verbose = FALSE,
    priors = priors,
    chains = 1
  )

  # Simulate new data
  pp <- predict(fit)

  ggplot2::ggplot(pp, ggplot2::aes(time, y_new, group = .draw)) +
    ggplot2::geom_line(alpha = 0.1) +
    ggplot2::theme_bw()
}
}

}
\seealso{
Obtaining predictions
\code{\link{fitted.dynamitefit}()}
}
\concept{prediction}

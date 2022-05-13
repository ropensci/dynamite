% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_priors.R
\name{get_data}
\alias{get_data}
\alias{get_data.dynamiteformula}
\alias{get_data.dynamitefit}
\title{Extract the Model Data of the Dynamite Model}
\usage{
get_data(x, data, group, time, ...)

\method{get_data}{dynamiteformula}(x, data, group, time, ...)

\method{get_data}{dynamitefit}(x, ...)
}
\arguments{
\item{x}{[\code{dynamiteformula} or \code{dynamitefit}]\cr The model formula or
existing \code{dynamitefit} object. See \code{\link[=dynamiteformula]{dynamiteformula()}} and \code{\link[=dynamite]{dynamite()}}.}

\item{data}{[\code{data.frame}]\cr The data frame containing the variables in
the model.}

\item{group}{[\code{character(1)}]\cr A column name of \code{data} that denotes the
unique groups.}

\item{time}{[\code{character(1)}]\cr A column name of \code{data} that denotes the
time axis.}

\item{...}{Ignored.}
}
\description{
Returns the input data to Stan model. Mostly useful for debugging.
}
\examples{

d <- data.frame(y = rnorm(10), x = 1:10, time = 1:10, id = 1)
str(get_data(obs(y ~ x, family = gaussian()),
  data = d, time = "time", group = "id"))
}
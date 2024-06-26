% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lfo.R
\name{plot.lfo}
\alias{plot.lfo}
\title{Diagnostic Plot for Pareto k Values from LFO}
\usage{
\method{plot}{lfo}(x, ...)
}
\arguments{
\item{x}{[\code{lfo}]\cr Output of the \code{lfo} method.}

\item{...}{Ignored.}
}
\value{
A \code{ggplot} object.
}
\description{
Plots Pareto k values per each time point (with one point per group),
together with a horizontal line representing the used threshold.
}
\examples{
data.table::setDTthreads(1) # For CRAN
\donttest{
# Please update your rstan and StanHeaders installation before running
# on Windows
if (!identical(.Platform$OS.type, "windows")) {
  # This gives warnings due to the small number of iterations
  plot(suppressWarnings(
    lfo(gaussian_example_fit, L = 20, chains = 1, cores = 1)
  ))
}
}

}

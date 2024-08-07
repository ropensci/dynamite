% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamite-package.R
\docType{data}
\name{gaussian_example}
\alias{gaussian_example}
\title{Simulated Data of a Gaussian Response}
\format{
A data frame with 3000 rows and 5 variables:
\describe{
\item{y}{The response variable.}
\item{x}{A continuous covariate.}
\item{z}{A binary covariate.}
\item{id}{Variable defining individuals (1 to 50).}
\item{time}{Variable defining the time point of the measurement (1 to 30).}
}
}
\source{
The data was generated via \code{gaussian_example.R} in
\url{https://github.com/ropensci/dynamite/tree/main/data-raw/}
}
\usage{
gaussian_example
}
\description{
Simulated data containing a Gaussian response variable \code{y} with two
covariates. The dataset was generated from a model with time-varying effects
of covariate \code{x} and the lagged value of the response variable, time-varying
intercept, and time-invariant effect of covariate \code{z}. The time-varying
coefficients vary according to a spline with 20 degrees of freedom.
}
\seealso{
Example models
\code{\link{categorical_example}},
\code{\link{categorical_example_fit}},
\code{\link{gaussian_example_fit}},
\code{\link{multichannel_example}},
\code{\link{multichannel_example_fit}}
}
\concept{examples}
\keyword{datasets}

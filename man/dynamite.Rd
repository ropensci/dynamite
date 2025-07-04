% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamite.R, R/print.R, R/summary.R
\name{dynamite}
\alias{dynamite}
\alias{formula.dynamitefit}
\alias{print.dynamitefit}
\alias{summary.dynamitefit}
\title{Estimate a Bayesian Dynamic Multivariate Panel Model}
\usage{
dynamite(
  dformula,
  data,
  time,
  group = NULL,
  priors = NULL,
  backend = "cmdstanr",
  verbose = TRUE,
  verbose_stan = FALSE,
  stanc_options = list("O0"),
  threads_per_chain = 1L,
  grainsize = NULL,
  custom_stan_model = NULL,
  debug = NULL,
  interval = 1L,
  ...
)

\method{formula}{dynamitefit}(x, ...)

\method{print}{dynamitefit}(x, full_diagnostics = FALSE, ...)

\method{summary}{dynamitefit}(object, ...)
}
\arguments{
\item{dformula}{[\code{dynamiteformula}]\cr The model formula.
See \code{\link[=dynamiteformula]{dynamiteformula()}} and 'Details'.}

\item{data}{[\code{data.frame}, \code{tibble::tibble}, or \code{data.table::data.table}]\cr
The data that contains the variables in the model in long format.
Supported column types are \code{integer}, \code{logical}, \code{double}, and
\code{factor}. Columns of type \code{character} will be converted to factors.
Unused factor levels will be dropped. The \code{data} can contain missing
values which will simply be ignored in the estimation in a case-wise
fashion (per time-point and per channel). Input \code{data} is converted to
channel specific matrix representations via \code{\link[stats:model.matrix]{stats::model.matrix.lm()}}.}

\item{time}{[\code{character(1)}]\cr A column name of \code{data} that denotes the
time index of observations. If this variable is a factor, the integer
representation of its levels are used internally for defining the time
indexing.}

\item{group}{[\code{character(1)}]\cr A column name of \code{data} that denotes the
unique groups or \code{NULL} corresponding to a scenario without any groups.
If \code{group} is \code{NULL}, a new column \code{.group} is created with constant
value \code{1L} is created indicating that all observations belong to the same
group. In case of name conflicts with \code{data}, see the \code{group_var} element
of the return object to get the column name of the new variable.}

\item{priors}{[\code{data.frame}]\cr An optional data frame with prior
definitions. See \code{\link[=get_priors]{get_priors()}} and 'Details'.}

\item{backend}{[\code{character(1)}]\cr Defines the backend interface to Stan,
should be  either \code{"cmdstanr"} (the default) or \code{"rstan"}. Note that
\code{cmdstanr} needs to be installed separately as it is not on CRAN. It also
needs the actual \code{CmdStan} software. See \url{https://mc-stan.org/cmdstanr/}
for details. Defaults to \code{"rstan"} if \code{"cmdstanr"} cannot be used.}

\item{verbose}{[\code{logical(1)}]\cr All warnings and messages are suppressed
if set to \code{FALSE}. Defaults to \code{TRUE}. Setting this to \code{FALSE} will also
disable checks for perfect collinearity in the model matrix.}

\item{verbose_stan}{[\code{logical(1)}]\cr This is the \code{verbose} argument for
\code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}}. Defaults to \code{FALSE}.}

\item{stanc_options}{[\code{list()}]\cr This is the \code{stanc_options} argument
passed to the compile method of a \code{CmdStanModel} object via
\code{cmdstan_model()} when \code{backend = "cmdstanr"}. Defaults to \code{list("O0")}.
To enable level one compiler optimizations, use \code{list("O1")}.
See \url{https://mc-stan.org/cmdstanr/reference/cmdstan_model.html}
for details.}

\item{threads_per_chain}{[\code{integer(1)}]\cr A Positive integer defining the
number of parallel threads to use within each chain. Default is \code{1}. See
\code{\link[rstan:rstan_options]{rstan::rstan_options()}} and
\url{https://mc-stan.org/cmdstanr/reference/model-method-sample.html}
for details.}

\item{grainsize}{[\code{integer(1)}]\cr A positive integer defining the
suggested size of the partial sums when using within-chain parallelization.
Default is number of time points divided by \code{threads_per_chain}.
Setting this to \code{1} leads the workload division entirely to the internal
scheduler. The performance of the within-chain parallelization can be
sensitive to the choice of \code{grainsize}, see Stan manual on reduce-sum for
details.}

\item{custom_stan_model}{[\code{character(1)}]\cr An optional character string
that either contains a customized Stan model code or a path to a \code{.stan}
file that contains the code. Using this will override the generated model
code. For expert users only.}

\item{debug}{[\code{list()}]\cr A named list of form \code{name = TRUE} indicating
additional objects in the environment of the \code{dynamite} function which are
added to the return object. Additionally, values \code{no_compile = TRUE} and
\code{no_sampling = TRUE} can be used to skip the compilation of the Stan code
and sampling steps respectively. This can be useful for debugging when
combined with \code{model_code = TRUE}, which adds the Stan model code to the
return object.}

\item{interval}{[\code{integer(1)}]\cr This arguments acts as an offset for
the evaluation of lagged observations when measurements are not available
at every time point. For example, if measurements are only available at
every second time point, setting \code{interval = 2} means that a lag of order
\code{k} will instead use the observation at \code{2 * k} time units in the past.
The default value is \code{1} meaning that there is a one-to-one correspondence
between the lag order and the time scale. For expert users only.}

\item{...}{For \code{dynamite()}, additional arguments to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}} or
the \verb{$sample()} method of the \code{CmdStanModel} object
(see \url{https://mc-stan.org/cmdstanr/reference/model-method-sample.html}),
such as \code{chains} and \code{cores}
(\code{chains} and \code{parallel_chains} in \code{cmdstanr}). For \code{summary()},
additional arguments to \code{\link[=as.data.frame.dynamitefit]{as.data.frame.dynamitefit()}}. For \code{print()},
further arguments to the print method for tibbles
(see \link[tibble:formatting]{tibble::formatting}). Not used for \code{formula()}.}

\item{x}{[\code{dynamitefit}]\cr The model fit object.}

\item{full_diagnostics}{By default, the effective sample size (ESS) and Rhat
are computed only for the time- and group-invariant parameters
(\code{full_diagnostics = FALSE}). Setting this to \code{TRUE} computes ESS and Rhat
values for all model parameters, which can take some time for complex models.}

\item{object}{[\code{dynamitefit}]\cr The model fit object.}
}
\value{
\code{dynamite} returns a \code{dynamitefit} object which is a list containing
the following components:
\itemize{
\item \code{stanfit}\cr A \code{stanfit} object, see \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}} for details.
\item \code{dformulas}\cr A list of \code{dynamiteformula} objects for internal use.
\item \code{data}\cr A processed version of the input \code{data}.
\item \code{data_name}\cr Name of the input data object.
\item \code{stan}\cr A \code{list} containing various elements related to Stan model
construction and sampling.
\item \code{group_var}\cr Name of the variable defining the groups.
\item \code{time_var}\cr Name of the variable defining the time index.
\item \code{priors}\cr Data frame containing the used priors.
\item \code{backend}\cr Either \code{"rstan"} or \code{"cmdstanr"} indicating which
package was used in sampling.
\item \code{permutation}\cr Randomized permutation of the posterior draws.
\item \code{call}\cr Original function call as an object of class \code{call}.
}

\code{formula} returns a quoted expression.

\code{print} returns \code{x} invisibly.

\code{summary} returns a \code{data.frame}.
}
\description{
Fit a Bayesian dynamic multivariate panel model (DMPM) using Stan for
Bayesian inference. The \pkg{dynamite} package supports a wide range of
distributions and allows the user to flexibly customize the priors for the
model parameters. The dynamite model is specified using standard \R formula
syntax via \code{\link[=dynamiteformula]{dynamiteformula()}}. For more information and examples,
see 'Details' and the package vignettes.

The \code{formula} method returns the model definition as a quoted expression.

Information on the estimated \code{dynamite} model can be obtained via
\code{print()} including the following: The model formula, the data,
the smallest effective sample sizes, largest Rhat and summary statistics of
the time-invariant and group-invariant model parameters.

The \code{summary()} method provides statistics of the posterior samples of the
model; this is an alias of \code{\link[=as.data.frame.dynamitefit]{as.data.frame.dynamitefit()}} with
\code{summary = TRUE}.
}
\details{
The best-case scalability of \code{dynamite} in terms of data size should be
approximately linear in terms of number of time points and and number of
groups, but as wall-clock time of the MCMC algorithms provided by Stan can
depend on the discrepancy of the data and the model (and the subsequent
shape of the posterior), this can vary greatly.
}
\examples{
data.table::setDTthreads(1) # For CRAN
\donttest{
# Please update your rstan and StanHeaders installation before running
# on Windows
if (!identical(.Platform$OS.type, "windows")) {
  fit <- dynamite(
    dformula = obs(y ~ -1 + varying(~x), family = "gaussian") +
      lags(type = "varying") +
      splines(df = 20),
    gaussian_example,
    "time",
    "id",
    chains = 1,
    refresh = 0
  )
}
}

data.table::setDTthreads(1) # For CRAN
formula(gaussian_example_fit)

data.table::setDTthreads(1) # For CRAN
print(gaussian_example_fit)

data.table::setDTthreads(1) # For CRAN
summary(gaussian_example_fit,
  types = "beta",
  probs = c(0.05, 0.1, 0.9, 0.95)
)

}
\references{
Santtu Tikka and Jouni Helske (2024). \pkg{dynamite}: An \R Package for
Dynamic Multivariate Panel Models. arXiv preprint,
\url{doi:10.48550/arXiv.2302.01607}.

Jouni Helske and Santtu Tikka (2022). Estimating Causal Effects
from Panel Data with Dynamic Multivariate Panel Models.
\emph{Advances in Life Course Research}, 60, 100617.
\url{doi:10.1016/j.alcr.2024.100617}.
}
\seealso{
Model fitting
\code{\link{dynamice}()},
\code{\link{get_priors}()},
\code{\link{update.dynamitefit}()}

Model formula construction
\code{\link{dynamiteformula}()},
\code{\link{lags}()},
\code{\link{lfactor}()},
\code{\link{random_spec}()},
\code{\link{splines}()}

Model outputs
\code{\link{as.data.frame.dynamitefit}()},
\code{\link{as.data.table.dynamitefit}()},
\code{\link{as_draws_df.dynamitefit}()},
\code{\link{coef.dynamitefit}()},
\code{\link{confint.dynamitefit}()},
\code{\link{get_code}()},
\code{\link{get_data}()},
\code{\link{get_parameter_dims}()},
\code{\link{get_parameter_names}()},
\code{\link{get_parameter_types}()},
\code{\link{ndraws.dynamitefit}()},
\code{\link{nobs.dynamitefit}()}
}
\concept{fitting}
\concept{formulas}
\concept{output}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamiteformula.R
\name{dynamiteformula}
\alias{dynamiteformula}
\alias{obs}
\alias{aux}
\alias{+.dynamiteformula}
\alias{print.dynamiteformula}
\title{Model Formula for \pkg{dynamite}}
\usage{
dynamiteformula(formula, family, link = NULL)

obs(formula, family, link = NULL)

aux(formula)

\method{+}{dynamiteformula}(e1, e2)

\method{print}{dynamiteformula}(x, ...)
}
\arguments{
\item{formula}{[\code{formula}]\cr An \R formula describing the model.}

\item{family}{[\code{character(1)}]\cr The family name. See 'Details' for the
supported families.}

\item{link}{[\code{character(1)}]\cr The name of the link function to use or
\code{NULL}. See details for supported link functions and default values of
specific families.}

\item{e1}{[\code{dynamiteformula}]\cr A model formula specification.}

\item{e2}{[\code{dynamiteformula}]\cr A model formula specification.}

\item{x}{[\code{dynamiteformula}]\cr The model formula.}

\item{...}{Ignored.}
}
\value{
A \code{dynamiteformula} object.
}
\description{
Defines a new observational or a new auxiliary channel for the model using
standard \R formula syntax. Formulas of individual response variables can be
joined together via \code{+}. See 'Details' and the package vignettes for more
information. The function \code{obs} is a shorthand alias for \code{dynamiteformula},
and \code{aux} is a shorthand alias for
\code{dynamiteformula(formula, family = "deterministic")}.
}
\details{
Currently the \pkg{dynamite} package supports the following
distributions for the observations:
\itemize{
\item Categorical: \code{categorical} (with a softmax link using the first category
as reference). See the documentation of the \code{categorical_logit_glm} in the
Stan function reference manual \url{https://mc-stan.org/users/documentation/}.
\item Multinomial: \code{multinomial} (softmax link, first category is reference).
\item Gaussian: \code{gaussian} (identity link, parameterized using mean and standard
deviation).
\item Multivariate Gaussian: \code{mvgaussian} (identity link, parameterized using
mean vector, standard deviation vector and the Cholesky decomposition of
the correlation matrix).
\item Poisson: \code{poisson} (log-link, with an optional known offset variable).
\item Negative-binomial: \code{negbin} (log-link, using mean and dispersion
parameterization, with an optional known offset variable). See the
documentation on \code{NegBinomial2} in the Stan function reference manual.
\item Bernoulli: \code{bernoulli} (logit-link).
\item Binomial: \code{binomial} (logit-link).
\item Exponential: \code{exponential} (log-link).
\item Gamma: \code{gamma} (log-link, using mean and shape parameterization).
\item Beta: \code{beta} (logit-link, using mean and precision parameterization).
\item Student t: \code{student} (identity link, parameterized using degrees of
freedom, location and scale)
}

The models in the \pkg{dynamite} package are defined by combining the
channel-specific formulas defined via \R formula syntax.
Each channel is defined via the \code{obs} function, and the channels are
combined with \code{+}. For example a formula
\code{obs(y ~ lag(x), family = "gaussian") + obs(x ~ z, family = "poisson")}
defines a model with two channels;
first we declare that \code{y} is a Gaussian variable depending on a previous
value of \code{x} (\code{lag(x)}), and then we add a second channel declaring \code{x} as
Poisson distributed depending on some exogenous variable \code{z}
(for which we do not define any distribution).

Number of trials for binomial channels should be defined via a \code{trials}
model component, e.g., \code{obs(y ~ x + trials(n), family = "binomial")},
where \code{n} is a data variable defining the number of trials. For multinomial
channels, the number of trials is automatically defined to be the sum
of the observations over the categories, but can also be defined using
the \code{trials} component, for example for prediction.

Multivariate channels are defined by providing a single formula for all
components or by providing component-specific formulas separated by a \code{|}.
The response variables that correspond to the components should be joined by
\code{c()}. For instance, the following would define \code{c(y1, y2)} as multivariate
gaussian with \code{x} as a predictor for the mean of the first component and
\code{x} and \code{z} as predictors for the mean of the second component:
\code{obs(c(y1, y2) ~ x | x + z, family = "mvgaussian")}. A multinomial channel
should only have a single formula.

In addition to declaring response variables via \code{obs}, we can also use
the function \code{aux} to define auxiliary channels which are deterministic
functions of other variables. The values of auxiliary variables are computed
dynamically during prediction, making the use of lagged values and other
transformations possible. The function \code{aux} also does not use the
\code{family} argument, which is automatically set to \code{deterministic} and is a
special channel type of \code{obs}. Note that lagged values of deterministic
\code{aux} channels do not imply fixed time points. Instead they must be given
starting values using a special function \code{init} that directly initializes
the lags to specified values, or by \code{past} which computes the initial values
based on an R expression. Both \code{init} and \code{past} should appear on the
right hand side of the model formula, separated from the primary defining
expression via \code{|}.

The formula within \code{obs} can also contain an additional special
function \code{varying}, which defines the time-varying part of the model
equation, in which case we could write for example
\code{obs(x ~ z + varying(~ -1 + w), family = "poisson")}, which defines a model
equation with a constant intercept and time-invariant effect of \code{z}, and a
time-varying effect of \code{w}. We also remove the duplicate intercept with \code{-1}
in order to avoid identifiability issues in the model estimation
(we could also define a time varying intercept, in which case we would write
\verb{obs(x ~ -1 + z + varying(~ w), family = "poisson)}). The part of the formula
not wrapped with \code{varying} is assumed to correspond to the fixed part of the
model, so \code{obs(x ~ z + varying(~ -1 + w), family = "poisson")} is actually
identical to
\code{obs(x ~ -1 + fixed(~ z) + varying(~ -1 + w), family = "poisson")} and
\code{obs(x ~ fixed(~ z) + varying(~ -1 + w), family = "poisson")}.

When defining varying effects, we also need to define how the these
time-varying regression coefficient behave. For this, a \code{splines} component
should be added to the model, e.g.,
\verb{obs(x ~ varying(~ -1 + w), family = "poisson) + splines(df = 10)} defines a
cubic B-spline with 10 degrees of freedom for the time-varying coefficient
corresponding to the \code{w}. If the model contains multiple time-varying
coefficients, same spline basis is used for all coefficients, with unique
spline coefficients and their standard deviation.

If the desired model contains lagged predictors of each response in each
channel, these can be quickly added to the model as either time-invariant
or time-varying predictors via \code{\link[=lags]{lags()}} instead of writing them manually
for each channel.

It is also possible to define group-specific (random) effects term
using the special syntax \code{random()} similarly as \code{varying()}. For example,
\code{random(~1)} leads to a model where in addition to the common intercept,
each individual/group has their own intercept with zero-mean normal prior and
unknown standard deviation analogously with the typical mixed models. An
additional model component \code{\link[=random_spec]{random_spec()}} can be used to define
whether the random effects are allowed to correlate within and across
channels and whether to use centered or noncentered parameterization for
the random effects.
}
\examples{
data.table::setDTthreads(1) # For CRAN
# A single gaussian response channel with a time-varying effect of 'x',
# and a time-varying effect of the lag of 'y' using B-splines with
# 20 degrees of freedom for the coefficients of the time-varying terms.
obs(y ~ -1 + varying(~x), family = "gaussian") +
  lags(type = "varying") +
  splines(df = 20)

# A two-channel categorical model with time-invariant predictors
# here, lag terms are specified manually
obs(x ~ z + lag(x) + lag(y), family = "categorical") +
  obs(y ~ z + lag(x) + lag(y), family = "categorical")

# The same categorical model as above, but with the lag terms
# added using 'lags'
obs(x ~ z, family = "categorical") +
  obs(y ~ z, family = "categorical") +
  lags(type = "fixed")

# A multichannel model with a gaussian, Poisson and a Bernoulli response and
# an auxiliary channel for the logarithm of 'p' plus one
obs(g ~ lag(g) + lag(logp), family = "gaussian") +
  obs(p ~ lag(g) + lag(logp) + lag(b), family = "poisson") +
  obs(b ~ lag(b) * lag(logp) + lag(b) * lag(g), family = "bernoulli") +
  aux(numeric(logp) ~ log(p + 1))

data.table::setDTthreads(1) # For CRAN
obs(y ~ x, family = "gaussian") + obs(z ~ w, family = "exponential")

data.table::setDTthreads(1) # For CRAN
x <- obs(y ~ x + random(~ 1 + lag(d)), family = "gaussian") +
  obs(z ~ varying(~w), family = "exponential") +
  aux(numeric(d) ~ log(y) | init(c(0, 1))) +
  lags(k = 2) +
  splines(df = 5) +
  random_spec(correlated = FALSE)
print(x)

}
\seealso{
Model formula construction
\code{\link{dynamite}()},
\code{\link{lags}()},
\code{\link{lfactor}()},
\code{\link{random_spec}()},
\code{\link{splines}()}
}
\concept{formulas}

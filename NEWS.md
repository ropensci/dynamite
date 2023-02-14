# dynamite 1.2.0

  * Added support for the multivariate gaussian distribution via `"mvgaussian"`
    family in `obs`. See the documentation of `dynamiteformula` for details
    on how to define multivariate channels.
  * Latent factors were not previously used in predict by error, this is now 
    fixed. However, due to identifiability constraints no new group levels are 
    allowed with models using latent factors.
  * Response variable names of the channels are now processed to avoid
    invalid variable names in the generated Stan code.
  * Optimized prediction code by removing redundant expressions and using
    better indexing.

# dynamite 1.1.2

  * The argument `verbose_stan` is now ignored when `backend = "cmdstanr"`.
  * The `stanc_options` argument for defining compiler options when 
    using `cmdstanr` can now be controlled via `dynamite`.
  * Optimized column binding of `data.table` objects in `predict` leading to 
    faster computation.
  * The `update` method now checks if the `backend` has changed from the 
    original model fit.
  * The `update` method now properly recompiles the model (if necessary) in 
    cases where the `update` is used for already updated `dynamitefit` object.
  * Fixed a bug in the default prior definitions of intercept for families using 
    log-link which lead to a prior mean -Inf if all observations at the first 
    time point were zero.
  * Fixed some issues in the code generation of latent factor components.

# dynamite 1.1.1

  * `plot_deltas` and other plotting functions now throw an error if you try to 
    plot parameters of incorrect type with them.

# dynamite 1.1.0

  * `dynamite` now supports general group-level random effects. New `random()` 
    works analogously with `varying()` inside `obs()`, and the new optional
    `random_spec()` component can be used to define whether the random effects 
    should be correlated or not and whether to use noncentered parameterization.
  * The package no longer depends on the `bayesplot` package. Instead, `ggplot2`
    and `patchwork` packages are used for the `plot` method.
  * Argument order of the `dynamite` function has been changed: `time` now 
    precedes `group` and `backend` now precedes `verbose`. This change is also 
    reflected in the `get_data`, `get_priors`, and `get_code` functions.
  * Vectorized priors and various indexing variables are now passed as data to 
    Stan instead of being hard-coded in the generated model code.
  * The package now supports contemporaneous dependencies between channels 
    such that the dependency structure is acyclic. For example, having 
    `y ~ x` and `x ~ z` simultaneously is valid, but adding `z ~ y` to these 
    would result in a cycle.
  * The output of `mcmc_diagnostics` is now clearer.
  * The default of the argument `summary` was changed to `FALSE` in 
    `as.data.frame` and `as.data.table` methods, whereas it is now hard-coded 
    to `TRUE` in `summary` method. The column ordering of the output of these 
    methods was also changed so that the estimates are before the extra 
    columns such as `time`.
  * The standard deviation of the default priors for spline coefficient 
    standard deviations is now scaled based on the data analogously with 
    regression coefficients.
  * Added argument `parameters` to `as.data.frame` and similar methods as well 
    for the plotting functions. 
  * Added functions `get_parameter_types` and `get_parameter_names` for 
    extracting model parameter types and names respectively.

# dynamite 1.0.2

  * Fixed a name clash issue in Stan code generation.

# dynamite 1.0.1

  * The package no longer depends on the `data.table` development version.
  * Removed the grunfeld example from vignette due to CRAN's size restrictions.
  * `multichannel_example` and the corresponding fit was modified: The standard 
    deviation parameter of the Gaussian channel used in the data generation was 
    decreased in order to make the example in the vignette more interesting.
  * The latent factor model was also modified by removing the `random()` 
    component in order to reduce the size of the model fit object.
  * Fixed the name extraction of the supplied data.
  * `plot_deltas` no longer unnecessarily warns about missing values.

# dynamite 1.0.0

  * Increased the version number to 1.0.0 to reflect the fact that the package 
    is now fully functional and has successfully passed the rOpenSci review.

# dynamite 0.0.3

  * `get_prior`, `get_code`, and `get_data` now support case without `group` 
    argument, as per issue #48.
  * Fixed some typos and other issues in the vignette raised by @nicholasjclark 
    during the rOpenSci review process.
  * Added an example on simulating from the prior predictive distribution to the 
    documentation of `predict`.
  * Declarations now occur before statements in the generated Stan code.
  * Added support for `cmdstanr` via argument `backend` in `dynamite`.
  * Added a link to the contributing guidelines to README.
  * The package no longer depends on the development version of `rstan`.
  * Dropped R version dependency from 4.1.0 to 3.5.0.
  * Moved `dplyr` and `tidyr` to 'Suggests'.
  * `categorical_logit` is now used instead of `categorical_logit_glm` on older
    `rstan` and `cmdstanr` versions.
  * Random intercepts with `random` now also support centered parametrization.
  * Added more comments to the generated Stan code.
  * Fixed the output of `formula.dynamitefit` so that it is now compatible with 
    the `update` method. Also added the required `call` object to the 
    `dynamitefit` object.
  * Added `loo` and `lfo` methods for the dynamite models which can be used 
    for approximate leave-one-out and leave-future-out cross validation.
  * Cleaned up NAMESPACE.
  * The `env` argument of `data.table` is now used to avoid possible variable
    name conflicts.
  * Breaking change: The shrinkage parameter which was previously named as 
    `lambda` is now `xi` in order to free `lambda` for factor loadings 
    parameter as is customary in factor analysis.
  * Added support for correlated latent dynamic factors (modeled as splines).
  * `get_code` applied to fitted model now correctly returns only the model 
    code and not the stanmodel object.
  * Fixed the `.draw` column of the `as.data.frame` output.

# dynamite 0.0.2

  * Improved the memory usage of `predict` and `fitted` by separating the 
    simulated values from the predictors that are independent of the posterior 
    draws.
  * Added support for summarized predictions via a new argument `funs`, this
    can further significantly reduce memory usage when individual level 
    predictions are not of interest.

# dynamite 0.0.1

  * First version of `dynamite`

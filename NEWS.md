# dynamite 1.3.0

  * Added support for Student's t-distribution via `"student"` family in `obs()`.
  * Added support for the multinomial distribution via `"multinomial"` family in `obs()`. A `trials()` term is now mandatory for multinomial channels.
  * The generated Stan code now automatically switches between the array keyword syntax and the deprecated syntax based on the backend Stan version (see https://mc-stan.org/docs/reference-manual/brackets-array-syntax.html for details).
  * The presence of variables used in `trials()` and `offset()` is now properly checked in the data.
  * The model components `trials()` and `offset()` now function correctly in `predict()` when they contain response variables of the model.
  * Fixed the calculation of the number of observations in `nobs()` for models that have multivariate channels.
  * Fixed an issue in `predict()` with models that contained multivariate channels with random effects.
  * Scenarios that have zero non-missing observations at specific time indices are now handled properly in the Stan code generation.
  * The names of additional arguments passed to `rstan::sampling()` and the `sample()` method of the `cmdstanr` Stan model via `...` in the call to `dynamite` are now checked and unrecognized arguments will be ignored.
  * Added a new function `get_parameter_dims()` that returns the parameter dimensions of the Stan model for `"dynamitefit"` and `"dynamiteformula"` objects.
  * Group-level random effects are now supported also for categorical and multinomial channels.
  * Added a new vignette that describes how the package can be used to simulate data from a dynamic multivariate panel model.
  * Added a new vignette that describes how the default priors of the model parameters are defined.

# dynamite 1.2.1

  * Removed argument `noncentered_lambda` from `lfactor()` as this did not work as intended.
  * Added next observation carried backward imputation scheme for fixed predictors in predict as option `"nocb"`.
  * Changed naming of `omega` parameters, they now include also the channel name.
  * Fixed an issue related to channels with latent factors that did not not have any other predictors.
  * Improved efficiency of sum-to-zero constraints based post by @aaronjg on the Stan forums.
  * Fixed several issues related to Stan code generation for the multivariate gaussian distribution.
  * The package no longer uses `gregexec()` internally which made it dependent on R version 4.1.0 or higher. 
  * Corrected R version dependency to 3.6.0 or higher based on the package dependencies.

# dynamite 1.2.0

  * Added support for the multivariate gaussian distribution via `"mvgaussian"` family in `obs()`. See the documentation of the `dynamiteformula()` function for details on how to define multivariate channels.
  * Latent factors were not previously used in predict by error, this is now fixed. However, due to identifiability constraints no new group levels are allowed with models using latent factors.
  * Response variable names of the channels are now processed to avoid invalid variable names in the generated Stan code. Note that these variables names should be used when defining priors and when using methods of the `"dynamitefit"` class. You can use the functions `get_priors()` and `get_parameter_names()` to see the names that are available, as before.
  * Optimized prediction code by removing redundant expressions and using better indexing.

# dynamite 1.1.2

  * The argument `verbose_stan` is now ignored when `backend = "cmdstanr"`.
  * The `stanc_options` argument for defining compiler options when using `cmdstanr` can now be controlled via `dynamite()`.
  * Optimized column binding of `"data.table"` objects in `predict()` leading to faster computation.
  * The `update()` method now checks if the `backend` has changed from the original model fit.
  * The `update()` method now properly recompiles the model (if necessary) in cases where `update()` is used for already updated `"dynamitefit"` object.
  * Fixed a bug in the default prior definitions of intercept for families using log-link which lead to `-Inf` prior mean if all observations at the first time point were zero.
  * Fixed some issues in the code generation of latent factor components.

# dynamite 1.1.1

  * `plot_deltas()` and other plotting functions now throw an error if the user tries to plot parameters of an incorrect type with them.

# dynamite 1.1.0

  * `dynamite()` now supports general group-level random effects. New `random()` works analogously with `varying()` inside `obs()`, and the new optional `random_spec()` component can be used to define whether the random effects should be correlated or not and whether to use noncentered parameterization.
  * The package no longer depends on the `bayesplot` package. Instead, `ggplot2` and `patchwork` packages are used for the `plot` method.
  * Argument order of the `dynamite()` function has been changed: `time` now precedes `group` and `backend` now precedes `verbose`. This change is also reflected in the `get_data()`, `get_priors()`, and `get_code()` functions.
  * Vectorized priors and various indexing variables are now passed as data to Stan instead of being hard-coded in the generated model code.
  * The package now supports contemporaneous dependencies between channels such that the dependency structure is acyclic. For example, having `y ~ x` and `x ~ z` simultaneously is valid, but adding `z ~ y` to these would result in a cycle.
  * The output of `mcmc_diagnostics()` is now clearer.
  * The default value of the `summary` argument was changed to `FALSE` in `as.data.frame()` and `as.data.table()` methods, whereas it is now hard-coded to `TRUE` in the `summary()` method. The column ordering of the output of these methods was also changed so that the estimate columns are placed before the extra columns such as `time`.
  * The standard deviation of the default priors for spline coefficient standard deviations is now scaled based on the data analogously with regression coefficients.
  * Added argument `parameters` to `as.data.frame()` and similar methods as well for the plotting functions. 
  * Added functions `get_parameter_types()` and `get_parameter_names()` for extracting model parameter types and names respectively.

# dynamite 1.0.2

  * Fixed a name clash issue in Stan code generation.

# dynamite 1.0.1

  * The package no longer depends on the development version of the `data.table` package.
  * Removed the Grunfeld example from vignette due to CRAN file size restrictions.
  * `multichannel_example` and the corresponding fit was modified: The standard deviation parameter of the Gaussian channel used in the data generation was decreased in order to make the example in the vignette more interesting.
  * The latent factor model was also modified by removing the `random()` component in order to reduce the size of the model fit object.
  * Fixed the name extraction of the supplied data.
  * `plot_deltas()` no longer unnecessarily warns about missing values.

# dynamite 1.0.0

  * Increased the version number to 1.0.0 to reflect the fact that the package is now fully functional and has successfully passed the rOpenSci review.

# dynamite 0.0.3

  * `get_prior()`, `get_code()`, and `get_data()` now support case without `group` argument, as per issue #48.
  * Fixed some typos and other issues in the vignette raised by @nicholasjclark during the rOpenSci review process.
  * Added an example on simulating from the prior predictive distribution to the documentation of `predict()`.
  * Declarations now occur before statements in the generated Stan code.
  * Added support for `cmdstanr` via argument `backend` in `dynamite`.
  * Added a link to the contributing guidelines to README.
  * The package no longer depends on the development version of `rstan`.
  * Dropped R version dependency from 4.1.0 to 3.5.0.
  * Moved `dplyr` and `tidyr` to 'Suggests'.
  * `categorical_logit()` is now used instead of `categorical_logit_glm()` on older `rstan` and `cmdstanr` versions.
  * Random intercepts with `random()` now also support centered parametrization.
  * Added more comments to the generated Stan code.
  * Fixed the output of `formula.dynamitefit()` so that it is now compatible with the `update()` method. Also added the required `"call"` object to the `"dynamitefit"` object.
  * Added `loo()` and `lfo()` methods for the dynamite models which can be used for approximate leave-one-out and leave-future-out cross validation.
  * Cleaned up NAMESPACE.
  * The `env` argument of `data.table()` is now used to avoid possible variable name conflicts.
  * Breaking change: The shrinkage parameter which was previously named as `lambda` is now `xi` in order to free `lambda` for factor loadings parameter as is customary in factor analysis.
  * Added support for correlated latent dynamic factors (modeled as splines).
  * `get_code()` applied to fitted model now correctly returns only the model code and not the `stanmodel` object.
  * Fixed the `.draw` column of the `as.data.frame()` output.

# dynamite 0.0.2

  * Improved the memory usage of `predict()` and `fitted()` by separating the simulated values from the predictors that are independent of the posterior draws.
  * Added support for summarized predictions via a new argument `funs`, this can further significantly reduce memory usage when individual level predictions are not of interest.

# dynamite 0.0.1

  * First version of `dynamite`

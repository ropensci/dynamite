# dynamite 1.6.2

  * The default Stan backend is now CmdStan via `cmdstanr`. If CmdStan or `cmdstanr` are not available, the package will default to `rstan`.

# dynamite 1.6.1

  * The `interval` argument of `dynamite()` is now coerced to `integer` for efficiency.
  * Fixed an issue with the `interval` argument with lagged predictors that were not response variables.

# dynamite 1.6.0

  * Added a new argument `drop` for `predict()` and `fitted()` which can be set to `FALSE` to keep redundant columns of `newdata` in the output. The default value `TRUE` matches the previous behavior of removing all redundant columns.
  * Added a new argument `interval` for controlling the time interval of observations in `dynamite()`. See the documentation for more details.
  * Predictions are now only generated for observations with no missing predictors in `predict()` and `fitted()`. This change does not affect the simulation results, but warnings related to missing or non-finite values should no longer be unnecesarily generated during the simulation. 
  * Disabled messages from `data.table` related to longer computations when using `dynamite()`, `dynamice()`, `predict()`, and `fitted()`.

# dynamite 1.5.7

  * Fixed an issue with `predict()` when the time index variable had attributes.
  * Fixed an issue with `predict()` where `character` columns were not converted to factors. Both `predict()` and `fitted()` now ensure that `character` columns are converted to factors with the same levels as the corresponding variables had in the original data.

# dynamite 1.5.6

  * Variable names and factor levels are now checked and modified if needed for compatibility with Stan. Previously only response variable names were checked. It is also now possible to have spaces in variable names by quoting them.
  * Fixed an example that made the package depend on R version 4.1.0.

# dynamite 1.5.5

  * The package vignettes are now prerendered as some of them took a long time to build.
  * The vignette builder has been switched to `quarto`.
  * The package no longer depends on `prder`, `pryr`, and `rmarkdown`.

# dynamite 1.5.4

  * Obtaining the model parameter dimensions via `get_parameter_dims()` no longer requires a compiled Stan model. This leads to a significant performance improvement when applied to `dynamiteformula` objects.
  * Model fitting using `cmdstanr` backend no longer relies on `rstan::read_stan_csv()` to construct the fit object. Instead, the resulting `CmdStanMCMC` object is used directly. This should provide a substantial performance improvement in some instances. For `dynamice()`, samples from different imputed datasets are combined using `cmdstanr::as_cmdstan_fit()` instead.

# dynamite 1.5.3

  * Restored and updated the main package vignette. The vignette now also contains a real data example and information on multiple imputation.
  * The package data `gaussian_simulation_fit` has been removed to accommodate CRAN package size requirements. The code to generate the data is still available in the `data_raw` directory.

# dynamite 1.5.2

  * The main package vignette has been temporarily removed as it contained out-of-date information. Please see the arXiv preprint for up-to-date information instead: https://arxiv.org/abs/2302.01607

# dynamite 1.5.1

  * The `type` argument of `coef()` and `plot()` has been replaced by `types` accepting multiple types simultaneously, similar to `as.data.table()` and `as.data.frame()`.
  * The functions `plot_betas()`, `plot_deltas()`, `plot_nus()`, `plot_lambdas()` and `plot_psis()` have been deprecated and are now provided via the default plot method by selecting the appropriate `types`.
  * A new argument `plot_type` has been added to control what type of plot will be drawn by the `plot()` method. The default value `"default"` draws the posterior means and posterior intervals of all parameters. The old functionality of drawing posterior densities and traceplots is provided by the option `"trace"`.
  * The `plot()` method has gained the argument `n_params` to limit the amount of parameters drawn at once (per parameter type).
  * Both time-varying and time-invariant parameters can now be plotted simultaneously.
  * Fixed an issue with `predict()` and `fitted()` for multinomial responses.
  * Priors of the cutpoint parameters of the `cumulative` family are now customizable.
  * Both `factor` and `ordered factor` responses are now supported for `categorical` and `cumulative` families. In addition, `ordered factor` columns of `data` are no longer converted to `factor` columns.
  * Arguments that have the different names but the same functionality between `rstan` and `cmdstanr` can now be used interchangeably for either backend, such as `iter` and `iter_samples`.
  * The latent factor component was reparametrized for additional robustness. User-visible changes are related to priors: Instead of prior on the standard deviations `sigma_lambda` and `tau_psi`, prior is now defined on `zeta`, the sum of these, as well as on `kappa`, which is the proportion of `zeta` attributable to `sigma_lambda`.

# dynamite 1.5.0

  * Estimation of dynamic multivariate panel models with multiple imputation is now available via the function `dynamice()` which uses the `mice` package.
  * `predict` and `fitted` functions no longer permutes the posterior samples when all samples are used i.e. when `n_draws = NULL` (default). This also corrects the standard error estimates of `loo()`, which were not correct earlier due to the mixing of chains.
  * Added an argument `thin` for `loo()`, `predict()` and `fitted()` methods.
  * Print method now only prints the run time for the fastest and the slowest chain instead of all chains.
  * A new exported function `hmc_diagnostics()` is now available.
  * Added a vignette on `get_code()` and `get_data()` functions and how they can be used to modify the generated Stan code and perform variational Bayes inference.
  * Contemporaneous dependencies are now allowed between different components of multivariate distributions, e.g., `obs(c(y, x) ~ x | 1, family = "mvgaussian")`.
  * Ordered probit and logit regressions are now available via `obs(., family = "cumulative", link = "probit")` and `obs(., family = "cumulative", link = "logit")`, respectively.

# dynamite 1.4.11

  * The package now depends on `data.table` version 1.15.0 or higher and the `ggforce` package.
  * Added a `plot` method for `dynamiteformula` objects. This method draws a directed acyclic graph (DAG) of the model structure as a snapshot in time with timepoints from the past and the future equal to the highest-order lag dependency in the model as a `ggplot` object. Alternatively, setting the argument `tikz = TRUE` returns the DAG as a `character` string in TikZ format. See the documentation for more details.

# dynamite 1.4.10

  * The formula interface now prohibits additional invalid `fixed()`, `varying()`, and `random()` definitions in `obs()`.
  * Fixed an error in Stan code generation if an offset term was included in the model formula.
  * Fixed an issue when using `character` type `group` variables.

# dynamite 1.4.9

  * Added option to input a custom model code for `dynamite` which can be used to tweak some aspects of the model (no checks on the compatibility with the post processing are made).
  * Changed the default optimization level for `cmdstanr` backend to `O0`, as the `O1` is not necessarily stable in all cases.
  * Added a new argument `full_diagnostics` to the `print()` method which can be used to control the computation of the ESS and Rhat values. By default, these are now computed only for the time- and group-invariant parameters (which are also printed).
  * The `print()` method now also warns about possible divergences, treedepth saturation, and low E-BMFI.
  * Fixed an error related to `predict()` code generation.

# dynamite 1.4.8

  * Made several performance improvements to data parsing.
  * `dynamite()` will now retain the original column order of `data` in all circumstances.

# dynamite 1.4.7

  * Added a note on priors vignette regarding default priors for `tau` parameters.
  * Fixed `mcmc_diagnostics()` function so that HMC diagnostics are checked also for models run with the `cmdstanr` backend.
  
# dynamite 1.4.6

  * Fixed the construction of latent factors for categorical responses.

# dynamite 1.4.5

  * The `get_data()` method for `dynamitefit` objects now correctly uses the previously defined priors instead of the default ones.
  * Fixed a bug in indexing of random effect terms.
  * Limited the number of parallel threads used by the `data.table` package to 1 in examples, tests, and vignettes for CRAN.

# dynamite 1.4.4

  * Example of the `lfo()` method now uses a single chain and core to avoid a compatibility issue with CRAN.
  * Fixed `plot_nus()` for categorical responses.
  * Fixed an issue which caused an error in error message of `predict()` and `fitted()` methods when `newdata` contained duplicate time points within group.
  * Fixed an issue (#72) which caused NA ELPD value in `lfo()` in case of missing data.

# dynamite 1.4.3

  * Fixed an issue with `formula.dynamitefit()` with models defined using `lags()` with a vector `k` argument with more than one value.
  * Fixed an issue in the `lfo()` method which resulted wrong ELPD estimates in panel data setting.
  * Fixed an issue in the `lfo()` method which in case of lagged responses caused the ELPD computations to skip last time points.

# dynamite 1.4.2

  * Added further checks and fixes for backwards compatibility with Stan.
  * Fixed code generation for intercept-only categorical model.
  * Fixed code generation in the transformed data block to be backwards compatible with Stan.

# dynamite 1.4.1

  * Fixed an issue in `dynamite()` data parsing that caused substantial memory usage in some instances.
  * Fixed an issue with Stan code generation for categorical responses.
  * Fixed an issue with `formula.dynamitefit()` with models that had multinomial channels.
  * Fixed an issue with `formula.dynamitefit()` when the `df` argument of `splines()` was `NULL`.
  * Formulas with `trials()` and `offset()` terms are now properly parsed when using `lags()`.
  * Removed experimental shrinkage feature.

# dynamite 1.4.0

  * `dynamite()` now supports parallel computation via the reduce-sum functionality of Stan.
  * Fixed an issue in `predict()` that resulted in redundant `NAs produced` warnings.
  * Fixed an issue with `formula.dynamitefit()` with models that had multivariate channels.

# dynamite 1.3.3

  * Fixed a partial argument name issue in the internal `update()` method used by `lfo()`.

# dynamite 1.3.2

  * Fixed the regularization of the default priors so that they match with the priors vignette.
  * Fixed an issue with the `update()` method for model fit objects without a group variable.
  * Fixed an issue with the `update()` method in `lfo()`.
  * Fixed an issue with `"tau"` and `"tau_alpha"` type parameters with the `as_draws()` method for categorical responses.
  * Fixed an issue with Stan code generation for models with time-varying covariates for categorical responses.
  * Fixed an issue with `formula.dynamitefit()` when the model contained a `splines` component.

# dynamite 1.3.1

  * Fixed an incorrect URL in the main vignette.
  * `"dynamitefit"` objects no longer contain the data used for Stan sampling by default. This data can still be retrieved via `get_data()`.
  * Added a new package data `gaussian_simulation_fit` that includes the model fit of the `dynamite_simulation` vignette for the example with time-varying effects.
  * The package data `latent_factor_example` and `latent_factor_example_fit` have been removed to accommodate CRAN package size requirements. The code to generate these data is still available in the `data_raw` directory.
  * Fixed an issue with `formula.dynamitefit()` when the model formula contained a `lags` component or a `lfactor` component.

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

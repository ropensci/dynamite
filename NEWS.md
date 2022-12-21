# dynamite 1.0.1
  * The package no longer depends on the `data.table` development version.

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

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
  * `categorical_logit` is now used instead of `categorical_logit_glm` on older
    `rstan` and `cmdstanr` versions.
  * Random intercepts with `random` now also support centered parameterization.
  * Added more comments to the generated Stan code.
  * `get_code` applied to fitted model now correctly returns only the model 
    code and not the stanmodel object.
  * Fixed the `.draw` column of the `as_data_frame` output.

# dynamite 0.0.2

  * Improved the memory usage of `predict` and `fitted` by separating the 
    simulated values from the predictors that are independent of the posterior 
    draws.
  * Added support for summarized predictions via a new argument `funs`, this
    can further significantly reduce memory usage when individual level 
    predictions are not of interest.

# dynamite 0.0.1

  * First version of `dynamite`

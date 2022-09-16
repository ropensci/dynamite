# dynamite 0.0.2

* Improved the memory usage of `predict` and `fitted` by separating the 
  simulated values from the predictors that are independent of the posterior 
  draws.
* Added support for summarized predictions via a new argument `funs`, this
  can further significantly reduce memory usage when individual level 
  predictions are not of interest.

# dynamite 0.0.1

* First version of `dynamite`

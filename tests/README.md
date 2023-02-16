### Running extended tests

The testthat directory contains several long-running tests which are not run by 
default. Running all of these tests can take several 
hours on an ordinary laptop. These tests can be switched on by defining an 
environmental variables `DYNAMITE_EXTENDED_TESTS = "true"` and 
`DYNAMITE_SCALING_TESTS = "true"` (e.g., 
`Sys.setenv("DYNAMITE_EXTENDED_TESTS" = "true")` in R), and on GitHub Actions by 
adding `run-extended` and `run-scaling` to the commit message respectively.

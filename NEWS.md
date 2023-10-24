# singleRcapture 0.2.1.1

* Bugfix for tests failing with `noLongDouble`

# singleRcapture 0.2.1

* Fixed bugs in `IRLS` fitting when providing `weights` argument when calling
  `estimatePopsize`
* The `weightsAsCounts` option in `controlModel` now works properly, 
  `dfbeta` and `dfpopsize` decrease weight of selected row in a model matrix 
  instead of deleting it when this is set to `TRUE`
* `simulate` method now works for both family object (like `ztpoisson()`) and
  for objects returned by `estimatePopsize`
* Introduced `singleRStaticCountData` sub class for `singleRclass` and made
  `estimatePopsize` a method so that a new package `singleRcaptureExtra` 
  (under development) can make all necessary calculations for pop size estimation
  when providing object fitted by `countreg::zerotrunc` or `VGAM::vglm`/`VGAM::vgam`
* Some bugfixes for multicore bootstrap
* Code was re-factored to make further development/maintenance for the package
  much easier
* Update will be uploaded to `CRAN`
* `semiparametric` bootstrap now has a much faster sampling algorithm (that does the same job)

Unit tests:
* Reduced computational burden of unit tests
* Multicore tests will only be performed after `TEST_SINGLERCAPTURE_MULTICORE_DEVELOPER`
  is set to `"true"` via `Sys.setenv` and `_R_CHECK_LIMIT_CORES_` to `false`

# singleRcapture 0.2.0.1

* Added `offset` argument to `estimatePopsize`
* Added options for parallel computing in `bootstrap` and in `dfbeta`
* Added deviance for all negative binomial based models. 
  (NOTE: They are very slow for now and I believe it may change after I verify
  one theoretical results that will lead to significant speed increase for 
  these computations)
* Overhaul of starting points (new methods and added linear predictors start in `IRLS`)
* Code for weights in `IRLS` fitting was speed up
* Minor bugfixes

# singleRcapture 0.2.0

### The package is now at [CRAN](https://cran.r-project.org/package=singleRcapture)

* features and improvements:
  * Added final `Hurdleztnegbin` model
  * Vastly improved `redoPopSize` which now handles bootstrap on a fitted model
  non standard covariance matrixes `newdata` argument user supplied `coef` etc.
  * Added `predict.singleR` method which offers standard error for both `link`, 
  `response` as well as `mean` predictions
  * No unexpected warnings should occur now in main function when using
  the package correctly
  * All control arguments are now verified before being passed
  * Fitting is now more reliable
  * Added information about `stats::optim` error codes
  * Added warnings for functions computing deviance
  

* bugfixes:
  * fixed bugs occurring when using mathematical functions as part of formulas
  i.e. when setting formula to something like: `y ~ log(x) + I(x ^ t) + I(t ^ 2)`

# singleRcapture 0.1.4

* features
  * Added `ztoinegbin`, `oiztnegbin` and `ztHurdlenegbin` models
  * Added an optional arguments to all family-functions to specify a link 
  function for distribution parameters
  * Updated and standardized documentation
  * Added more warnings
  * Added some more methods for `singleR` class in some commonly used `glm` 
  functions, in particular `texreg::screenreg` should work well now

* changes
  * Changed some default arguments
  * Added option to save logs from `IRLS` fitting

* bugfixes
  * Fixed some issues with intercept only models
  * Fixed some slight miscalculations in information matrixes for one inflated
  models making fitting them much more reliable

* github repository
  * More and better `Rcmd` tests

# singleRcapture 0.1.3.2 -- NTTS
* features:
  * Added function that implements population size estimates for stratas
  * More warnings in fitting
  * More options in control functions
  * Corrected/implemented deviance residuals for more models
  
* changes:
  * Now the whole package uses `cammelCase`
  * Performance upgrades
  * Corrected some miss calculated moments
  * Change exported data so that factors are actually factors not just characters
  * Removed unused dependency
  
* github repository
  * Added automated `R-cmd` check

# singleRcapture 0.1.3.1
* features:
  * Basically all of documentation was redone and now features most of important 
  theory on SSCR methods and some information on (v)glms
  * Added checks on positivity of working weights matrixes to stabilise `"IRLS"` algorithm
  * Added most of sandwich capabilities to the package, in particular:
    * S3 method for `vcovHC` was implemented
    * `vcovCL` should work on `singleR` class objects 
    should work with `"HC0"` and `"HC1"` `type` argument values
  * Basic version of function `redoPopEstimation` for updating the
  population size estimation after post-hoc procedures was implemented
  * `popSizeEst` function for extracting population size estimation
  results was implemented
  * Minor improvements to memory usage were made and 
  computation was speed up a little
  * Changed names of mle and robust fitting methods to optim and IRLS
  respectively
  * Some bugfixes
  * More warnings messages in `estimate_popsize.fit`

# singleRcapture 0.1.3

* features:
  * Multiple new models
  * `IRLS` generalised for distributions with multiple parameters
  * bugfixes
  * QOL improvements
  * extended bootstrap and most other methods for new models


# singleRcapture 0.1.2

* features:
  * control parameters for model
  * control parameters for regression in bootstrap sampling
  * leave one out diagnostics for popsize and regression parameters (`dfbetas` were corrected)
  * fixes for Goodness of fit tests in zero one truncated models
  * computational improvements in `IRLS`
  * other small bugfixes

# singleRcapture 0.1.1

* bug fixes and some of the promised features for 0.2.0 in particular
  * More tiny tests 
  * Some fixes for marginal frequencies
  * Deviance implemented
  * dfbetas and levarage matrix
  * Parametric bootstraps work correctly for the most part there is just some polishing left to do

# singleRcapture 0.1.0 

* first version of the package released

# singleRcapture 0.2.0

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

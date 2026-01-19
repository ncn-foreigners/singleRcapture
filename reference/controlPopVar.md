# Control parameters for population size estimation

Creating control parameters for population size estimation and
respective standard error and variance estimation.

## Usage

``` r
controlPopVar(
  alpha = 0.05,
  bootType = c("parametric", "semiparametric", "nonparametric"),
  B = 500,
  confType = c("percentilic", "normal", "basic"),
  keepbootStat = TRUE,
  traceBootstrapSize = FALSE,
  bootstrapVisualTrace = FALSE,
  fittingMethod = c("optim", "IRLS"),
  bootstrapFitcontrol = NULL,
  sd = c("sqrtVar", "normalMVUE"),
  covType = c("observedInform", "Fisher"),
  cores = 1L
)
```

## Arguments

- alpha:

  a significance level, 0.05 used by default.

- bootType:

  the bootstrap type to be used. Default is `"parametric"`, other
  possible values are: `"semiparametric"` and `"nonparametric"`.

- B:

  a number of bootstrap samples to be performed (default 500).

- confType:

  a type of confidence interval for bootstrap confidence interval,
  `"percentile"` by default. Other possibilities: `"studentized"` and
  `"basic"`.

- keepbootStat:

  a boolean value indicating whether to keep a vector of statistics
  produced by bootstrap.

- traceBootstrapSize:

  a boolean value indicating whether to print size of bootstrapped
  sample after truncation for semi- and fully parametric bootstraps.

- bootstrapVisualTrace:

  a boolean value indicating whether to plot bootstrap statistics in
  real time if `cores = 1` if `cores > 1` it instead indicates whether
  to make progress bar.

- fittingMethod:

  a method used for fitting models from bootstrap samples.

- bootstrapFitcontrol:

  control parameters for each regression works exactly like
  `controlMethod` but for fitting models from bootstrap samples.

- sd:

  a character indicating how to compute standard deviation of population
  size estimator either as:
  \\\hat{\sigma}=\sqrt{\hat{\text{var}}(\hat{N})}\\ for `sqrt` (which is
  slightly biased if \\\hat{N}\\ has a normal distribution) or for
  `normalMVUE` as the unbiased minimal variance estimator for normal
  distribution: \\\hat{\sigma}=\sqrt{\hat{\text{var}}(\hat{N})}
  \frac{\Gamma\left(\frac{N\_{obs}-1}{2}\right)}{\Gamma\left(\frac{N\_{obs}}{2}\right)}
  \sqrt{\frac{N\_{obs}}{2}}\\ where the ration involving gamma functions
  is computed by log gamma function.

- covType:

  a type of covariance matrix for regression parameters by default
  observed information matrix.

- cores:

  for bootstrap only, a number of processor cores to be used, any number
  greater than 1 activates code designed with `doParallel`, `foreach`
  and `parallel` packages. Note that for now using parallel computing
  makes tracing impossible so `traceBootstrapSize` and
  `bootstrapVisualTrace` parameters are ignored in this case.

## Value

A list with selected parameters, it is also possible to call list
directly.

## See also

[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)
[`controlModel()`](https://ncn-foreigners.github.io/singleRcapture/reference/controlModel.md)
[`controlMethod()`](https://ncn-foreigners.github.io/singleRcapture/reference/controlMethod.md)

## Author

Piotr Chlebicki, Maciej Beręsewicz

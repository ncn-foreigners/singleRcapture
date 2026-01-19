# Predict method for singleRStaticCountData class

A method for `predict` function, works analogous to `predict.glm` but
gives the possibility to get standard errors of mean/distribution
parameters and directly get pop size estimates for new data.

## Usage

``` r
# S3 method for class 'singleRStaticCountData'
predict(
  object,
  newdata,
  type = c("response", "link", "mean", "popSize", "contr"),
  se.fit = FALSE,
  na.action = NULL,
  weights,
  cov,
  ...
)
```

## Arguments

- object:

  an object of `singleRStaticCountData` class.

- newdata:

  an optional `data.frame` containing new data.

- type:

  the type of prediction required, possible values are:

  - `"response"`– For matrix containing estimated distributions
    parameters.

  - `"link"` – For matrix of linear predictors.

  - `"mean"` – For fitted values of both \\Y\\ and \\Y\|Y\>0\\.

  - `"contr"` – For inverse probability weights (here named for
    observation contribution to population size estimate).

  - `"popSize"` – For population size estimation. Note this results in a
    call to `redoPopEstimation` and it is usually better to call this
    function directly.

  by default set to `"response"`.

- se.fit:

  a logical value indicating whether standard errors should be computed.
  Only matters for `type` in `"response", "mean", "link"`.

- na.action:

  does nothing yet.

- weights:

  optional vector of weights for `type` in `"contr", "popSize"`.

- cov:

  optional matrix or function or character specifying either a
  covariance matrix or a function to compute that covariance matrix. By
  default `vcov.singleRStaticCountData` can be set to e.g. `vcovHC`.

- ...:

  arguments passed to other functions, for now this only affects
  `vcov.singleRStaticCountData` method and `cov` function.

## Value

Depending on `type` argument if one of `"response", "link", "mean"` a
matrix with fitted values and possibly standard errors if `se.fit`
argument was set to `TRUE`, if `type` was set to `"contr"` a vector with
inverses of probabilities, finally for `"popSize"` an object of class
`popSizeEstResults` with its own methods containing population size
estimation results.

## Details

Standard errors are computed with assumption of regression coefficients
being asymptotically normally distributed, if this assumption holds then
each of linear predictors i.e. each row of
\\\boldsymbol{\eta}=\boldsymbol{X}\_{vlm}\boldsymbol{\beta}\\ is
asymptotically normally distributed and their variances are expressed by
well known formula. The mean \\\mu\\ and distribution parameters are
then differentiable functions of asymptotically normally distributed
variables and therefore their variances can be computed using
(multivariate) delta method.

## See also

[`redoPopEstimation()`](https://ncn-foreigners.github.io/singleRcapture/reference/redoPopEstimation.md)
[`stats::summary.glm()`](https://rdrr.io/r/stats/summary.glm.html)
[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)

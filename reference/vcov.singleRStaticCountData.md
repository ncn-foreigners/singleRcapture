# Obtain covariance matrix estimation

A `vcov` method for `singleRStaticCountData` class.

## Usage

``` r
# S3 method for class 'singleRStaticCountData'
vcov(object, type = c("Fisher", "observedInform"), ...)
```

## Arguments

- object:

  object of singleRStaticCountData class.

- type:

  type of estimate for covariance matrix for now either expected
  (Fisher) information matrix or observed information matrix.

- ...:

  additional arguments for method functions

## Value

A covariance matrix for fitted coefficients, rows and columns of which
correspond to parameters returned by `coef` method.

## Details

Returns a estimated covariance matrix for model coefficients calculated
from analytic hessian or Fisher information matrix usually utilizing
asymptotic effectiveness of maximum likelihood estimates. Covariance
type is taken from control parameter that have been provided on call
that created `object` if arguments `type` was not specified.

## See also

[`vcovHC.singleRStaticCountData()`](https://ncn-foreigners.github.io/singleRcapture/reference/vcovHC.singleRStaticCountData.md)
[`sandwich::sandwich()`](https://sandwich.R-Forge.R-project.org/reference/sandwich.html)

# Extract population size estimation results

An extractor function with `singleRStaticCountData` method for
extracting important information regarding pop size estimate.

## Usage

``` r
# S3 method for class 'singleRRatioReg'
popSizeEst(object, estimator = c("primary", "HT", "SM"), ...)

popSizeEst(object, ...)

# S3 method for class 'singleRStaticCountData'
popSizeEst(object, ...)
```

## Arguments

- object:

  object with population size estimates.

- estimator:

  for `singleRRatioReg` objects, which population-size estimator to
  return. The default `"primary"` returns the estimator selected when
  calling
  [`ratioReg()`](https://ncn-foreigners.github.io/singleRcapture/reference/ratioReg.md).

- ...:

  additional optional arguments, currently not used in
  `singleRStaticCountData` class method.

## Value

An object of class `popSizeEstResults` containing population size
estimation results.

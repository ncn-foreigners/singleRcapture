# Confidence intervals for model parameters

A function that computes studentized confidence intervals for model
coefficients.

## Usage

``` r
# S3 method for class 'singleRStaticCountData'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  object of singleRStaticCountData class.

- parm:

  names of parameters for which confidence intervals are to be computed,
  if missing all parameters will be considered.

- level:

  confidence level for intervals.

- ...:

  currently does nothing.

## Value

An object with named columns that include upper and lower limit of
confidence intervals.

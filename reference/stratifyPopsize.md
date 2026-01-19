# Estimate size of sub-populations

A function that estimates sizes of specific sub populations based on a
capture-recapture model for the whole population.

## Usage

``` r
stratifyPopsize(object, strata, alpha, ...)

# S3 method for class 'singleRStaticCountData'
stratifyPopsize(object, strata, alpha, cov = NULL, ...)
```

## Arguments

- object:

  an object on which the population size estimates should be based in
  `singleRcapture` package this is a fitter `singleRStaticCountData`
  class object.

- strata:

  a specification of sub populations given by one of:

  - formula – a formula to be applied to `model.frame` extracted from
    the object.

  - Logical vector with number of entries equal to number of rows in the
    dataset.

  - A (named) list where each element is a logical vector, names of the
    list will be used to specify names variable in returned object.

  - Vector of names of explanatory variables. For
    `singleRStaticCountData` method for this function this specification
    of `strata` parameter will result in every level of explanatory
    variable having its own sub population for each variable specified.

  - If no value was provided the `singleRStaticCountData` method for
    this function will itself create sub populations based on levels of
    factor variables in `model.frame`.

- alpha:

  significance level for confidence intervals – Either a single numeric
  value or a vector of length equal to number of sub populations
  specified in `strata`. If missing it is set to `.05` in
  `singleRStaticCountData` method.

- ...:

  a vector of arguments to be passed to other functions. For
  `singleRStaticCountData` method for this functions arguments in `...`
  are passed to either `cov` if argument provided was a function or
  `vcov` if `cov` argument was missing at call.

- cov:

  for `singleRStaticCountData` method an estimate of variance-covariance
  matrix for estimate of regression parameters. It is possible to pass a
  function such as for example
  [`sandwich::vcovHC`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)
  which will be called as: `foo(object, ...)` and a user may specify
  additional arguments of a function in `...` argument. If not provided
  an estimate for covariance matrix will be set by calling appropriate
  `vcov` method.

## Value

A `data.frame` object with row names being the names of specified sub
populations either provided or inferred.

## Details

In single source capture-recapture models the most frequently used
estimate for population size is Horvitz-Thompson type estimate:

\\\hat{N} = \sum\_{k=1}^{N}\frac{I\_{k}}{\mathbb{P}(Y\_{k}\>0)} =
\sum\_{k=1}^{N\_{obs}}\frac{1}{1-\mathbb{P}(Y\_{k}=0)}\\

where \\I\_{k}=I\_{Y\_{k} \> 0}\\ are indicator variables, with value 1
if kth unit was observed at least once and 0 otherwise and the inverse
probabilistic weights weights for units observed in the data
\\\tfrac{1}{\mathbb{P}(Y\_{k}\>0)}\\ are estimated using fitted linear
predictors.

The estimates for different sub populations are made by changing the
\\I\_{k}=I\_{Y\_{k} \> 0}\\ indicator variables to refer not to the
population as a whole but to the sub populations that are being
considered i.e. by changing values from 1 to 0 if kth unit is not a
member of sub population that is being considered at the moment.

The estimation of variance for these estimates and estimation of
variance for estimate of population size for the whole population follow
the same relation as the one described above.

## See also

[`vcov.singleRStaticCountData()`](https://ncn-foreigners.github.io/singleRcapture/reference/vcov.singleRStaticCountData.md)
[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)

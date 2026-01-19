# Summary statistics for model of singleRStaticCountData class

A `summary` method for `singleRStaticCountData` class

## Usage

``` r
# S3 method for class 'singleRStaticCountData'
summary(
  object,
  test = c("t", "z"),
  resType = "pearson",
  correlation = FALSE,
  confint = FALSE,
  cov,
  popSizeEst,
  ...
)
```

## Arguments

- object:

  object of singleRStaticCountData class.

- test:

  type of test for significance of parameters `"t"` for t-test and `"z"`
  for normal approximation of students t distribution, by default `"z"`
  is used if there are more than 30 degrees of freedom and `"t"` is used
  in other cases.

- resType:

  type of residuals to summarize any value that is allowed in
  `residuals.singleRStaticCountData` except for `"all"` is allowed. By
  default pearson residuals are used.

- correlation:

  logical value indicating whether correlation matrix should be computed
  from covariance matrix by default `FALSE`.

- confint:

  logical value indicating whether confidence intervals for regression
  parameters should be constructed. By default `FALSE`.

- cov:

  covariance matrix corresponding to regression parameters. It is
  possible to give `cov` argument as a function of `object`. If not
  specified it will be constructed using `vcov.singleRStaticCountData`
  method. (i.e using Cramer-Rao lower bound)

- popSizeEst:

  a `popSizeEstResults` class object. If not specified population size
  estimation results will be drawn from `object`. If any post-hoc
  procedures, such as sandwich covariance matrix estimation or bias
  reduction, were taken it is possible to include them in population
  size estimation results by calling `redoPopEstimation`.

- ...:

  additional optional arguments passed to the following functions:

  - `vcov.singleRStaticCountData` – if no `cov` argument was provided.

  - `cov` – if `cov` parameter specified at call was a function.

  - `confint.singleRStaticCountData` – if `confint` parameter was set to
    `TRUE` at function call. In particular it is possible to set
    confidence level in `...`.

## Value

An object of `summarysingleRStaticCountData` class containing:

- `call` – A call which created `object`.

- `coefficients` – A dataframe with estimated regression coefficients
  and their summary statistics such as standard error Wald test
  statistic and p value for Wald test.

- `residuals` – A vector of residuals of type specified at call.

- `aic` – Akaike's information criterion.

- `bic` – Bayesian (Schwarz's) information criterion.

- `iter` – Number of iterations taken in fitting regression.

- `logL` – Logarithm of likelihood function evaluated at coefficients.

- `deviance` – Residual deviance.

- `populationSize` – Object with population size estimation results.

- `dfResidual` – Residual degrees of freedom.

- `sizeObserved` – Size of observed population.

- `correlation` – Correlation matrix if `correlation` parameter was set
  to `TRUE`

- `test` – Type of statistical test performed.

- `model` – Family class object specified in call for `object`.

- `skew` – If bootstrap sample was saved contains estimate of skewness.

## Details

Works analogically to `summary.glm` but includes population size
estimation results. If any additional statistics, such as confidence
intervals for coefficients or coefficient correlation, are specified
they will be printed.

## See also

[`redoPopEstimation()`](https://ncn-foreigners.github.io/singleRcapture/reference/redoPopEstimation.md)
[`stats::summary.glm()`](https://rdrr.io/r/stats/summary.glm.html)

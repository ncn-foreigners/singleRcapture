# Control parameters for regression

`controlMethod` constructs a list with all necessary control parameters
for regression fitting in `estimatePopsizeFit` and `estimatePopsize`.

## Usage

``` r
controlMethod(
  epsilon = 1e-08,
  maxiter = 1000,
  verbose = 0,
  printEveryN = 1L,
  coefStart = NULL,
  etaStart = NULL,
  optimMethod = "Nelder-Mead",
  silent = FALSE,
  optimPass = FALSE,
  stepsize = 1,
  checkDiagWeights = TRUE,
  weightsEpsilon = 1e-08,
  momentumFactor = 0,
  saveIRLSlogs = FALSE,
  momentumActivation = 5,
  criterion = c("coef", "abstol", "reltol")
)
```

## Arguments

- epsilon:

  a tolerance level for fitting algorithms by default `1e-8`.

- maxiter:

  a maximum number of iterations.

- verbose:

  a numeric value indicating whether to trace steps of fitting algorithm
  for `IRLS` fitting method different values of verbose give the
  following information:

  - 1 – Returns information on the number of current iteration and
    current log-likelihood.

  - 2 – Returns information on vector of regression parameters at
    current iteration (and all of the above).

  - 3 – Returns information on reduction of log-likelihood at current
    iteration (and all of the above).

  - 4 – Returns information on value of log-likelihood function gradient
    at current iteration (and all of the above).

  - 5 – Returns information on convergence criterion and values that are
    taken into account when considering convergence (and all of the
    above).

  if `optim` method was chosen verbose will be passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) as trace.

- printEveryN:

  an integer value indicating how often to print information specified
  in `verbose`, by default set to `1`.

- coefStart, etaStart:

  initial parameters for regression coefficients or linear predictors if
  `NULL`. For `IRLS` fitting only `etaStart` is needed so if `coefStart`
  is provided it will be converted to `etaStart`, for `optim` fitting
  `coefStart` is necessary and argument `etaStart` will be ignored.

- optimMethod:

  a method of [`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  used `"Nelder-Mead"` is the default .

- silent:

  a logical value, indicating whether warnings in `IRLS` method should
  be suppressed.

- optimPass:

  an optional list of parameters passed to
  `stats::optim(..., control = optimPass)` if FALSE then list of control
  parameters will be inferred from other parameters.

- stepsize:

  only for `IRLS`, scaling of updates to `beta` vector lower value means
  slower convergence but more accuracy by default 1. In general if
  fitting algorithm fails lowering this value tends to be most effective
  at correcting it.

- checkDiagWeights:

  a logical value indicating whether to check if diagonal elements of
  working weights matrixes in `IRLS` are sufficiently positive so that
  these matrixes are positive defined. By default `TRUE`.

- weightsEpsilon:

  a small number to ensure positive definedness of weights matrixes.
  Only matters if `checkDiagWeights` is set to `TRUE`. By default
  `1e-8`.

- momentumFactor:

  an experimental parameter in `IRLS` only allowing for taking previous
  step into account at current step, i.e instead of updating regression
  parameters as: \\\boldsymbol{\beta}\_{(a)} =
  \boldsymbol{\beta}\_{(a-1)} + \text{stepsize} \cdot
  \text{step}\_{(a)}\\ the update will be made as: \\
  \boldsymbol{\beta}\_{(a)} = \boldsymbol{\beta}\_{(a-1)} +
  \text{stepsize} \cdot (\text{step}\_{(a)} +
  \text{momentum}\cdot\text{step}\_{(a-1)})\\

- saveIRLSlogs:

  a logical value indicating if information specified in `verbose`
  should be saved to output object, by default `FALSE`.

- momentumActivation:

  the value of log-likelihood reduction bellow which momentum will
  apply.

- criterion:

  a criterion used to determine convergence in `IRLS`, multiple values
  may be provided. By default `c("coef", "abstol")`.

## Value

List with selected parameters, it is also possible to call list
directly.

## See also

[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)
[`estimatePopsizeFit()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsizeFit.md)
[`controlModel()`](https://ncn-foreigners.github.io/singleRcapture/reference/controlModel.md)
[`controlPopVar()`](https://ncn-foreigners.github.io/singleRcapture/reference/controlPopVar.md)

## Author

Piotr Chlebicki, Maciej Beręsewicz

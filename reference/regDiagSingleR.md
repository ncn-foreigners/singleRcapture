# Regression diagnostics in singleRcapture

List of some regression diagnostics implemented for
`singleRStaticCountData` class. Functions that either require no changes
from `glm` class or are not relevant to context of `singleRcapture` are
omitted.

## Usage

``` r
dfpopsize(model, ...)

# S3 method for class 'singleRStaticCountData'
dfpopsize(model, dfbeta = NULL, ...)

# S3 method for class 'singleRStaticCountData'
dfbeta(model, maxitNew = 1, trace = FALSE, cores = 1, ...)

# S3 method for class 'singleRStaticCountData'
hatvalues(model, ...)

# S3 method for class 'singleRStaticCountData'
residuals(
  object,
  type = c("pearson", "pearsonSTD", "response", "working", "deviance", "all"),
  ...
)

# S3 method for class 'singleRStaticCountData'
cooks.distance(model, ...)

# S3 method for class 'singleRStaticCountData'
sigma(object, ...)

# S3 method for class 'singleRStaticCountData'
influence(model, do.coef = FALSE, ...)

# S3 method for class 'singleRStaticCountData'
rstudent(model, ...)

# S3 method for class 'singleRStaticCountData'
rstandard(model, type = c("deviance", "pearson"), ...)
```

## Arguments

- model, object:

  an object of `singleRStaticCountData` class.

- ...:

  arguments passed to other methods. Notably
  `dfpopsize.singleRStaticCountData` calls
  `dfbeta.singleRStaticCountData` if no `dfbeta` argument was provided
  and `controlMethod` is called in `dfbeta` method.

- dfbeta:

  if `dfbeta` was already obtained it is possible to pass them into
  function so that they need not be computed for the second time.

- maxitNew:

  the maximal number of iterations for regressions with starting points
  \\\hat{\boldsymbol{\beta}}\\ on data specified at call for `model`
  after the removal of k'th row. By default 1.

- trace:

  a logical value specifying whether to tracking results when
  `cores > 1` it will result in a progress bar being created.

- cores:

  a number of processor cores to be used, any number greater than 1
  activates code designed with `doParallel`, `foreach` and `parallel`
  packages. Note that for now using parallel computing makes tracing
  impossible so `trace` parameter is ignored in this case.

- type:

  a type of residual to return.

- do.coef:

  logical indicating if `dfbeta` computation for influence should be
  done. `FALSE` by default.

## Value

- For `hatvalues` – A matrix with n rows and p columns where n is a
  number of observations in the data and p is number of regression
  parameters.

- For `dfpopsize` – A vector for which k'th element corresponds to the
  difference between point estimate of population size estimation on
  full data set and point estimate of population size estimation after
  the removal of k'th unit from the data set.

- For `dfbeta` – A matrix with n rows and p observations where p is a
  number of units in data and p is the number of regression parameters.
  K'th row of this matrix corresponds to
  \\\hat{\boldsymbol{\beta}}-\hat{\boldsymbol{\beta}}\_{-k}\\ where
  \\\hat{\boldsymbol{\beta}}\_{-k}\\ is a vector of estimates for
  regression parameters after the removal of k'th row from the data.

- `cooks.distance` – A matrix with a single columns with values of cooks
  distance for every unit in `model.matrix`

- `residuals.singleRStaticCountData` – A `data.frame` with chosen
  residuals.

## Details

`dfpopsize` and `dfbeta` are closely related. `dfbeta` fits a regression
after removing a specific row from the data and returns the difference
between regression coefficients estimated on full data set and data set
obtained after deletion of that row, and repeats procedure once for
every unit present in the data.`dfpopsize` does the same for population
size estimation utilizing coefficients computed by `dfbeta`.

`cooks.distance` is implemented (for now) only for models with a single
linear predictor and works exactly like the method for `glm` class.

`sigma` computes the standard errors of predicted means. Returns a
matrix with two columns first for truncated mean and the other for the
non-truncated mean.

`residuals.singleRStaticCountData` (can be abbreviated to `resid`) works
like `residuals.glm` with the exception that:

- `"pearson"` – returns non standardized residuals.

- `"pearsonSTD"` – is currently defined only for single predictors
  models but will be extended to all models in a near future, but for
  families with more than one distribution parameter it will be a
  multivariate residual.

- `"response"` – returns both residuals computed with truncated and non
  truncated fitted value.

- `"working"` – is possibly multivariate if more than one linear
  predictor is present.

- `"deviance"` – is not yet defined for all families in
  [`singleRmodels()`](https://ncn-foreigners.github.io/singleRcapture/reference/singleRmodels.md)
  e.g. negative binomial based methods.

- `"all"` – returns all available residual types.

`hatvalues.singleRStaticCountData` is method for
`singleRStaticCountData` class for extracting diagonal elements of
projection matrix.

Since `singleRcapture` supports not only regular glm's but also vglm's
the `hatvalues` returns a matrix with number of columns corresponding to
number of linear predictors in a model, where kth column corresponds to
elements of the diagonal of projection matrix associated with kth linear
predictor. For glm's \\\boldsymbol{W}^{\frac{1}{2}}\boldsymbol{X}
\left(\boldsymbol{X}^{T}\boldsymbol{W}\boldsymbol{X}\right)^{-1}
\boldsymbol{X}^{T}\boldsymbol{W}^{\frac{1}{2}}\\ where:
\\\boldsymbol{W}=\mathbb{E}\left(\text{Diag}
\left(\frac{\partial^{2}\ell}{\partial\boldsymbol{\eta}^{T}
\partial\boldsymbol{\eta}}\right)\right)\\ and \\\boldsymbol{X}\\ is a
model (lm) matrix. For vglm's present in the package it is instead :
\\\boldsymbol{X}\_{vlm}
\left(\boldsymbol{X}\_{vlm}^{T}\boldsymbol{W}\boldsymbol{X}\_{vlm}\right)^{-1}
\boldsymbol{X}\_{vlm}^{T}\boldsymbol{W}\\ where: \\ \boldsymbol{W} =
\mathbb{E}\left(\begin{bmatrix}
\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta\_{1}^{T}\partial\eta\_{1}}\right)
&
\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta\_{1}^{T}\partial\eta\_{2}}\right)
& \dotso &
\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta\_{1}^{T}\partial\eta\_{p}}\right)\cr
\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta\_{2}^{T}\partial\eta\_{1}}\right)
&
\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta\_{2}^{T}\partial\eta\_{2}}\right)
& \dotso &
\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta\_{2}^{T}\partial\eta\_{p}}\right)\cr
\vdots & \vdots & \ddots & \vdots\cr
\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta\_{p}^{T}\partial\eta\_{1}}\right)
&
\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta\_{p}^{T}\partial\eta\_{2}}\right)
& \dotso &
\text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta\_{p}^{T}\partial\eta\_{p}}\right)
\end{bmatrix}\right)\\ is a block matrix constructed by taking the
expected value from diagonal matrixes corresponding to second
derivatives with respect to each linear predictor (and mixed
derivatives) and \\\boldsymbol{X}\_{vlm}\\ is a model (vlm) matrix
constructed using specifications in `controlModel` and call to
`estimatePopsize`.

`influence` works like `glm` counterpart computing the most important
influence measures.

## See also

[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)
[`stats::hatvalues()`](https://rdrr.io/r/stats/influence.measures.html)
[`controlMethod()`](https://ncn-foreigners.github.io/singleRcapture/reference/controlMethod.md)
[`stats::dfbeta()`](https://rdrr.io/r/stats/influence.measures.html)
[`stats::cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)

## Author

Piotr Chlebicki, Maciej Beręsewicz

## Examples

``` r
# \donttest{
# For singleRStaticCountData class
# Get simple model
Model <- estimatePopsize(
  formula = capture ~ nation + age + gender, 
  data = netherlandsimmigrant, 
  model = ztpoisson, 
  method = "IRLS"
)
# Get dfbeta
dfb <- dfbeta(Model)
# The dfpopsize results are obtained via (It is also possible to not provide 
# dfbeta then they will be computed manually):
res <- dfpopsize(Model, dfbeta = dfb)
summary(res)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -4236.412     2.663     2.663     5.448    17.284   117.448 
plot(res)

# see vaious types of residuals:
head(resid(Model, "all"))
#>       lambda truncatedResponse nontruncatedResponse    pearson pearsonSTD
#> 1 -0.9328875       -0.25364898            0.5294668 -0.4864422 -0.4867283
#> 2 -0.9328875       -0.25364898            0.5294668 -0.4864422 -0.4867283
#> 3 -0.9328875       -0.25364898            0.5294668 -0.4864422 -0.4867283
#> 4 -0.9791760       -0.06666172            0.8695135 -0.2554869 -0.2559964
#> 5 -0.9791760       -0.06666172            0.8695135 -0.2554869 -0.2559964
#> 6  2.7449807        0.74635102            1.5294668  1.4313347  1.4321768
#>     deviance
#> 1 -0.6992491
#> 2 -0.6992491
#> 3 -0.6992491
#> 4 -0.3631875
#> 5 -0.3631875
#> 6  1.0619767
# }
```

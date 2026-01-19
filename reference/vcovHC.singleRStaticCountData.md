# Heteroscedasticity-consistent covariance matrix estimation for singleRStaticCountData class

S3 method for `vcovHC` to handle `singleRStaticCountData` class objects.
Works exactly like `vcovHC.default` the only difference being that this
method handles vector generalised linear models. Updating the covariance
matrix in variance/standard error estimation for population size
estimator can be done via
[`redoPopEstimation()`](https://ncn-foreigners.github.io/singleRcapture/reference/redoPopEstimation.md)

## Usage

``` r
# S3 method for class 'singleRStaticCountData'
estfun(x, ...)

# S3 method for class 'singleRStaticCountData'
bread(x, ...)

# S3 method for class 'singleRStaticCountData'
vcovHC(
  x,
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
  omega = NULL,
  sandwich = TRUE,
  ...
)
```

## Arguments

- x:

  a fitted `singleRStaticCountData` class object.

- ...:

  for `vcovHC` additional optional arguments passed to the following
  functions:

  - `estfun` – for empirical estimating functions.

  - `hatvalues` – for diagonal elements of projection matrix.

  - `sandwich` – only if `sandwich` argument in function call was set to
    `TRUE`.

  - `vcov` – when calling `bread` internally.

- type:

  a character string specifying the estimation type, same as in
  [`sandwich::vcovHC.default`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html).
  HC3 is the default value.

- omega:

  a vector or a function depending on the arguments residuals (i.e. the
  derivative of log-likelihood with respect to each linear predictor),
  diaghat (the diagonal of the corresponding hat matrix) and df (the
  residual degrees of freedom), same as in
  [`sandwich::vcovHC.default`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html).

- sandwich:

  logical. Should the sandwich estimator be computed? If set to FALSE
  only the meat matrix is returned. Same as in
  [`sandwich::vcovHC()`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)

## Value

Variance-covariance matrix estimation corrected for heteroscedasticity
of regression errors

## See also

[`sandwich::vcovHC()`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)
[`redoPopEstimation()`](https://ncn-foreigners.github.io/singleRcapture/reference/redoPopEstimation.md)

## Author

Piotr Chlebicki, Maciej Beręsewicz

## Examples

``` r
set.seed(1)
N <- 10000
gender <- rbinom(N, 1, 0.2)
eta <- -1 + 0.5*gender
counts <- rpois(N, lambda = exp(eta))
df <- data.frame(gender, eta, counts)
df2 <- subset(df, counts > 0)
mod1 <-  estimatePopsize(
  formula = counts ~ 1 + gender, 
  data = df2, 
  model = "ztpoisson", 
  method = "optim", 
  popVar = "analytic"
)
require(sandwich)
HC <- sandwich::vcovHC(mod1, type = "HC4")
Fisher <- vcov(mod1, "Fisher") # variance covariance matrix obtained from 
#Fisher (expected) information matrix
HC
#>             (Intercept)       gender
#> (Intercept)  0.00201216 -0.002012160
#> gender      -0.00201216  0.004790297
Fisher
#>              (Intercept)       gender
#> (Intercept)  0.002022145 -0.002022145
#> gender      -0.002022145  0.004881858
# usual results
summary(mod1)
#> 
#> Call:
#> estimatePopsize.default(formula = counts ~ 1 + gender, data = df2, 
#>     model = "ztpoisson", method = "optim", popVar = "analytic")
#> 
#> Pearson Residuals:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -0.5503 -0.4267 -0.4267  0.0000 -0.4267  6.2261 
#> 
#> Coefficients:
#> -----------------------
#> For linear predictors associated with: lambda 
#>             Estimate Std. Error z value  P(>|z|)    
#> (Intercept) -1.01346    0.04497 -22.537  < 2e-16 ***
#> gender       0.50288    0.06987   7.197 6.14e-13 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> AIC: 3990.538
#> BIC: 4002.804
#> Residual deviance: 2386.29
#> 
#> Log-likelihood: -1993.269 on 3403 Degrees of freedom 
#> Number of calls to log-likelihood function: 35
#> -----------------------
#> Population size estimation results: 
#> Point estimate 10146.02
#> Observed proportion: 33.6% (N obs = 3405)
#> Std. Error 341.7287
#> 95% CI for the population size:
#>           lowerBound upperBound
#> normal      9476.239   10815.79
#> logNormal   9508.827   10849.72
#> 95% CI for the share of observed population:
#>           lowerBound upperBound
#> normal      31.48175   35.93198
#> logNormal   31.38330   35.80883
# updated results
summary(mod1, cov = HC,
popSizeEst = redoPopEstimation(mod1, cov = HC))
#> 
#> Call:
#> estimatePopsize.default(formula = counts ~ 1 + gender, data = df2, 
#>     model = "ztpoisson", method = "optim", popVar = "analytic")
#> 
#> Pearson Residuals:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -0.5503 -0.4267 -0.4267  0.0000 -0.4267  6.2261 
#> 
#> Coefficients:
#> -----------------------
#> For linear predictors associated with: lambda 
#>             Estimate Std. Error z value  P(>|z|)    
#> (Intercept) -1.01346    0.04486 -22.593  < 2e-16 ***
#> gender       0.50288    0.06921   7.266 3.71e-13 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> AIC: 3990.538
#> BIC: 4002.804
#> Residual deviance: 2386.29
#> 
#> Log-likelihood: -1993.269 on 3403 Degrees of freedom 
#> Number of calls to log-likelihood function: 35
#> -----------------------
#> Population size estimation results: 
#> Point estimate 10146.02
#> Observed proportion: 33.6% (N obs = 3405)
#> Std. Error 340.7902
#> 95% CI for the population size:
#>           lowerBound upperBound
#> normal      9478.078   10813.95
#> logNormal   9510.489   10847.69
#> 95% CI for the share of observed population:
#>           lowerBound upperBound
#> normal      31.48710   35.92500
#> logNormal   31.38916   35.80257
# estimating equations
mod1_sims <- sandwich::estfun(mod1)
head(mod1_sims)
#>   (Intercept)    gender
#> 1  -0.1924342 0.0000000
#> 2   0.6700925 0.6700925
#> 3  -0.1924342 0.0000000
#> 4  -0.1924342 0.0000000
#> 5  -0.1924342 0.0000000
#> 6  -0.1924342 0.0000000
# bread method
all(vcov(mod1, "Fisher") * nrow(df2) == sandwich::bread(mod1, type = "Fisher"))
#> [1] TRUE
```

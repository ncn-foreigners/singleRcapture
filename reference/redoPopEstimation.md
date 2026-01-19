# Updating population size estimation results

A function that applies all post-hoc procedures that were taken (such as
heteroscedastic consistent covariance matrix estimation or bias
reduction) to population size estimation and standard error estimation.

## Usage

``` r
redoPopEstimation(object, newdata, ...)

# S3 method for class 'singleRStaticCountData'
redoPopEstimation(
  object,
  newdata,
  cov,
  weights,
  coef,
  control,
  popVar,
  offset,
  weightsAsCounts,
  ...
)
```

## Arguments

- object:

  object for which update of population size estimation results will be
  done.

- newdata:

  optional `data.frame` with new data for pop size estimation.

- ...:

  additional optional arguments, currently not used in
  `singleRStaticCountData` class method.

- cov:

  an updated covariance matrix estimate.

- weights:

  optional vector of weights to use in population size estimation.

- coef:

  optional vector of coefficients of regression on which to base
  population size estimation. If missing it is set to `coef(object)`.

- control:

  similar to `controlPopVar` in
  [`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md).
  If missing set to controls provided on call to `object`.

- popVar:

  similar to `popVar` in
  [`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md).
  If missing set to `"analytic"`.

- offset:

  offset argument for new data

- weightsAsCounts:

  for `singleRStaticCountData` method used to specify whether weights
  should be treated as number of occurrences for rows in data

## Value

An object of class `popSizeEstResults` containing updated population
size estimation results.

## Details

Any non specified arguments will be inferred from the `object`

## Examples

``` r
# Create simple model
Model <- estimatePopsize(
  formula = capture ~ nation + gender, 
  data = netherlandsimmigrant, 
  model = ztpoisson, 
  method = "IRLS"
)
# Apply heteroscedasticity consistent covariance matrix estimation
require(sandwich)
#> Loading required package: sandwich
cov <- vcovHC(Model, type = "HC3")
summary(Model, cov = cov,
popSizeEst = redoPopEstimation(Model, cov = cov))
#> 
#> Call:
#> estimatePopsize.default(formula = capture ~ nation + gender, 
#>     data = netherlandsimmigrant, model = ztpoisson, method = "IRLS")
#> 
#> Pearson Residuals:
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#> -0.4796685 -0.4796685 -0.3518333  0.0004493 -0.2257173 14.0588575 
#> 
#> Coefficients:
#> -----------------------
#> For linear predictors associated with: lambda 
#>                      Estimate Std. Error z value  P(>|z|)    
#> (Intercept)           -1.3977     0.2553  -5.475 4.37e-08 ***
#> nationAsia            -1.0560     0.3940  -2.680  0.00736 ** 
#> nationNorth Africa     0.2327     0.2399   0.970  0.33200    
#> nationRest of Africa  -0.8864     0.3514  -2.523  0.01164 *  
#> nationSurinam         -2.3519     1.0273  -2.289  0.02205 *  
#> nationTurkey          -1.6845     0.6110  -2.757  0.00583 ** 
#> gendermale             0.3833     0.2014   1.904  0.05695 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> AIC: 1718.993
#> BIC: 1757.766
#> Residual deviance: 1136.645
#> 
#> Log-likelihood: -852.4963 on 1873 Degrees of freedom 
#> Number of iterations: 8
#> -----------------------
#> Population size estimation results: 
#> Point estimate 11879.92
#> Observed proportion: 15.8% (N obs = 1880)
#> Std. Error 2531.196
#> 95% CI for the population size:
#>           lowerBound upperBound
#> normal      6918.872   16840.98
#> logNormal   8015.862   18177.38
#> 95% CI for the share of observed population:
#>           lowerBound upperBound
#> normal      11.16325   27.17206
#> logNormal   10.34252   23.45350
# Compare to results with usual covariance matrix estimation
summary(Model)
#> 
#> Call:
#> estimatePopsize.default(formula = capture ~ nation + gender, 
#>     data = netherlandsimmigrant, model = ztpoisson, method = "IRLS")
#> 
#> Pearson Residuals:
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#> -0.4796685 -0.4796685 -0.3518333  0.0004493 -0.2257173 14.0588575 
#> 
#> Coefficients:
#> -----------------------
#> For linear predictors associated with: lambda 
#>                      Estimate Std. Error z value  P(>|z|)    
#> (Intercept)           -1.3977     0.2146  -6.514 7.33e-11 ***
#> nationAsia            -1.0560     0.3016  -3.501 0.000464 ***
#> nationNorth Africa     0.2327     0.1939   1.200 0.230002    
#> nationRest of Africa  -0.8864     0.3009  -2.946 0.003224 ** 
#> nationSurinam         -2.3519     1.0137  -2.320 0.020337 *  
#> nationTurkey          -1.6845     0.6029  -2.794 0.005203 ** 
#> gendermale             0.3833     0.1630   2.352 0.018686 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> AIC: 1718.993
#> BIC: 1757.766
#> Residual deviance: 1136.645
#> 
#> Log-likelihood: -852.4963 on 1873 Degrees of freedom 
#> Number of iterations: 8
#> -----------------------
#> Population size estimation results: 
#> Point estimate 11879.92
#> Observed proportion: 15.8% (N obs = 1880)
#> Std. Error 2448.792
#> 95% CI for the population size:
#>           lowerBound upperBound
#> normal      7080.380   16679.47
#> logNormal   8111.333   17927.69
#> 95% CI for the share of observed population:
#>           lowerBound upperBound
#> normal      11.27134   26.55225
#> logNormal   10.48657   23.17745

## get confidence interval with larger significance level
redoPopEstimation(Model, control = controlPopVar(alpha = .000001))
#> Point estimate: 11879.92
#> Variance: 5996583
#> 99.9999% confidence intervals:
#>           lowerBound upperBound
#> normal      1880.000   23858.53
#> logNormal   4951.313   34438.87
```

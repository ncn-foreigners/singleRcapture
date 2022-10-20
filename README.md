
<!-- README.md is generated from README.Rmd. Please edit that file -->

# singleRcapture

<!-- badges: start -->
<!-- badges: end -->

The goal of singleRcapture is to automate single source
capture-recapture estimation of population size.

## Funding

Work on this package is supported by the the National Science Center,
OPUS 22 grant no. 2020/39/B/HS4/00941.

## Installation

You can install the development version of singleRcapture from
[GitHub](https://github.com/ncn-foreigners/singleRcapture) with:

``` r
# install.packages("devtools")
devtools::install_github("ncn-foreigners/singleRcapture")
```

### Examples

This is a basic example of zero truncated poisson model and zelterman
model with netherlands imigrant data with analytic variance:

``` r
library(singleRcapture)
ModelPo <- estimate_popsize(formula = capture ~ .,
                            data = netherlandsimmigrant,
                            pop.var = "analytic",
                            model = "ztpoisson",
                            method = "robust")
ModelZl <- estimate_popsize(formula = capture ~ .,
                            data = netherlandsimmigrant,
                            pop.var = "analytic",
                            model = "zelterman",
                            method = "robust")
summary(ModelPo)
#> estimate_popsize(formula = capture ~ ., data = netherlandsimmigrant, 
#>     model = "ztpoisson", method = "robust", pop.var = "analytic")
#> 
#> Pearson Residuals:
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -0.488779 -0.486085 -0.297859  0.002075 -0.210439 13.921578 
#> 
#> Coefficients:
#> -----------------------
#> For linear predictors associated with: lambda 
#>                      Estimate Std. Error z value P(>|z|)    
#> (Intercept)            -2.317      0.449   -5.16 2.5e-07 ***
#> gender                  0.397      0.163    2.44 1.5e-02   *
#> age                     0.975      0.408    2.39 1.7e-02   *
#> reason                  0.011      0.162    0.07 9.5e-01    
#> nationAsia             -1.092      0.302   -3.62 2.9e-04 ***
#> nationNorth Africa      0.190      0.194    0.98 3.3e-01    
#> nationRest of Africa   -0.911      0.301   -3.03 2.5e-03  **
#> nationSurinam          -2.337      1.014   -2.31 2.1e-02   *
#> nationTurkey           -1.675      0.603   -2.78 5.5e-03  **
#> -----------------------
#> Signif. codes:  0 '****' 0.001 '***' 0.01 '**' 0.05 '*' 0.1 '.' 1 ' '
#> 
#> AIC: 1714.896
#> BIC: 1764.747
#> Residual deviance: 1128.549
#> 
#> Log-likelihood: -848.4481 on 1871 Degrees of freedom 
#> Number of iterations: 8
#> -----------------------
#> Population size estimation results: 
#> Point estimate 12691.45
#> Observed proportion: 14.8% (N obs = 1880)
#> Std. Error 2809.508
#> 95% CI for the population size:
#>              lowerBound upperBound
#> Studentized    7184.917   18197.99
#> Logtransform   8430.749   19723.38
#> 95% CI for the share of observed population:
#>              lowerBound upperBound
#> Studentized   10.330814   26.16592
#> Logtransform   9.531836   22.29932
```

``` r
summary(ModelZl)
#> estimate_popsize(formula = capture ~ ., data = netherlandsimmigrant, 
#>     model = "zelterman", method = "robust", pop.var = "analytic")
#> 
#> Pearson Residuals:
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.47496 -0.43206 -0.32340  0.04906 -0.19071  9.54717 
#> 
#> Coefficients:
#> -----------------------
#> For linear predictors associated with: lambda 
#>                      Estimate Std. Error z value P(>|z|)    
#> (Intercept)            -3.359      0.528   -6.36 2.0e-10 ***
#> gender                  0.535      0.232    2.30 2.1e-02   *
#> age                     0.567      0.434    1.31 1.9e-01    
#> reason                  0.189      0.220    0.86 3.9e-01    
#> nationAsia             -1.056      0.448   -2.36 1.8e-02   *
#> nationNorth Africa      0.579      0.307    1.89 5.9e-02   .
#> nationRest of Africa   -0.664      0.425   -1.56 1.2e-01    
#> nationSurinam          -1.720      1.050   -1.64 1.0e-01    
#> nationTurkey           -1.030      0.657   -1.57 1.2e-01    
#> -----------------------
#> Signif. codes:  0 '****' 0.001 '***' 0.01 '**' 0.05 '*' 0.1 '.' 1 ' '
#> 
#> AIC: 1133.029
#> BIC: 1182.627
#> Residual deviance: 1115.029
#> 
#> Log-likelihood: -557.5143 on 1819 Degrees of freedom 
#> Number of iterations: 7
#> -----------------------
#> Population size estimation results: 
#> Point estimate 16188.3
#> Observed proportion: 11.6% (N obs = 1880)
#> Std. Error 3166.094
#> 95% CI for the population size:
#>              lowerBound upperBound
#> Studentized    9982.871   22393.73
#> Logtransform  11201.447   23843.06
#> 95% CI for the share of observed population:
#>              lowerBound upperBound
#> Studentized    8.395207   18.83226
#> Logtransform   7.884896   16.78355
```

Marginal frequencies and Goodness of fit test:

``` r
summary(marginalFreq(ModelPo), df = 2, dropl5 = "group")
#> Test for Goodness of fit of a regression model:
#> 
#>                  Test statistics df P(>X^2)
#> Chi-squared test           50.06  2 1.4e-11
#> G-test                     34.31  2 3.6e-08
#> 
#> -------------------------------------------------------------- 
#> Cells with fitted frequencies of < 5 have been grouped 
#> Names of cells used in calculating test(s) statistic: 1 2 3
```

Should goodness of fit tests prove insufficient in determining the best
model graphical comparisons may also be made. One such technique is a
rootogram and it is implemented as a part of singleR method for `plot`
function:

``` r
plot(ModelPo, plotType = "rootogram", main = "Zero Truncated Poisson model")
plot(ModelZl, plotType = "rootogram", main = "Logistic regression based Zelterman model")
```

<img src="man/figures/README-plot-1.png" width="45%" /><img src="man/figures/README-plot-2.png" width="45%" />

singleRcapture also includes bootstraps and models truncated at values 0
and 1 and non standard confidence levels

``` r
set.seed(123)
zotgeomBoot <- estimate_popsize(
    formula = TOTAL_SUB ~ .,
    data = farmsubmission,
    pop.var = "bootstrap",
    model = "zotgeom",
    method = "robust",
    control.pop.var = control.pop.var(B = 1000, alpha = .01)
)
summary(zotgeomBoot)
#> estimate_popsize(formula = TOTAL_SUB ~ ., data = farmsubmission, 
#>     model = "zotgeom", method = "robust", pop.var = "bootstrap", 
#>     control.pop.var = control.pop.var(B = 1000, alpha = 0.01))
#> 
#> Pearson Residuals:
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -0.952107 -0.727540 -0.426714  0.003655  0.322783 16.127909 
#> 
#> Coefficients:
#> -----------------------
#> For linear predictors associated with: lambda 
#>              Estimate Std. Error z value  P(>|z|)    
#> (Intercept)    -2.608      0.298   -8.76  2.0e-18 ***
#> log_size        0.585      0.022   26.47 2.2e-154 ***
#> log_distance   -0.068      0.025   -2.66  7.7e-03  **
#> C_TYPE          0.611      0.044   13.73  6.6e-43 ***
#> -----------------------
#> Signif. codes:  0 '****' 0.001 '***' 0.01 '**' 0.05 '*' 0.1 '.' 1 ' '
#> 
#> AIC: 19483.08
#> BIC: 19509.67
#> Residual deviance: 23179.43
#> 
#> Log-likelihood: -9737.539 on 5692 Degrees of freedom 
#> Number of iterations: 7
#> -----------------------
#> Population size estimation results: 
#> Point estimate 29087.96
#> Observed proportion: 41.4% (N obs = 12036)
#> Bootstrap Std. Error 1962.695
#> 99% CI for the population size:
#> lowerBound upperBound 
#>   25228.32   35576.04 
#> 99% CI for the share of observed population:
#> lowerBound upperBound 
#>   33.83176   47.70830
```

``` r
plot(zotgeomBoot, plotType = "bootHist", labels = TRUE, ylim = c(0, 425))
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

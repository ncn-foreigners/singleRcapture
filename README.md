
<!-- README.md is generated from README.Rmd. Please edit that file -->

# singleRcapture

<!-- badges: start -->
<!-- badges: end -->

The goal of singleRcapture is to automate single source
capture-recapture estimation of population size.

## Funding

Work on this package is supported by the the National Science Center, OPUS 22 grant no. 2020/39/B/HS4/00941.

## Installation

You can install the development version of singleRcapture from
[GitHub](https://github.com/ncn-foreigners/singleRcapture) with:

``` r
# install.packages("devtools")
devtools::install_github("ncn-foreigners/singleRcapture")
```

## Example

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
#> Response Residuals:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.5250  0.5302  0.8227  0.8578  0.9114  5.5302 
#> 
#> Coefficients:
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
#> Deviance: 3688.542
#> 
#> Log-likelihood: -848.4481 on 1871Degrees of freedom 
#> Number of iterations: 653
#> -----------------------
#> Population size estimation results: 
#> Point estimate 12691.36
#> Variance 7892651
#> Std. Error 2809.386
#> 95% CI:
#>              lowerBound upperBound
#> Studentized    7185.061   18197.65
#> Logtransform   8430.801   19722.92
```

``` r
summary(ModelZl)
#> estimate_popsize(formula = capture ~ ., data = netherlandsimmigrant, 
#>     model = "zelterman", method = "robust", pop.var = "analytic")
#> 
#> Response Residuals:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.5488  0.6266  0.8699  0.9286  0.9273  5.6266 
#> 
#> Coefficients:
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
#> BIC: 1182.88
#> Deviance: 0
#> 
#> Log-likelihood: -557.5143 on 1819Degrees of freedom 
#> Number of iterations: 8
#> -----------------------
#> Population size estimation results: 
#> Point estimate 16188.3
#> Variance 10024148
#> Std. Error 3166.094
#> 95% CI:
#>              lowerBound upperBound
#> Studentized    9982.871   22393.73
#> Logtransform  11201.447   23843.06
```

Marginal frequencies and Goodness of fit test:

``` r
summary(marginalFreq(ModelPo), df = 2)
#> Test for Goodness of fit of a regression model:
#> 
#>                  Test statistics df P(>X^2)
#> Chi-squared test           87.48  2 1.0e-19
#> G-test                     37.31  2 7.9e-09
#> 
#> --------------------------------------------------------
#> Cells with fitted frequencies of < 5 have been preserved
```

Here is a plot of marginal frequencies with matplot:

``` r
a1 <- marginalFreq(ModelPo)$table[-1]
a2 <- marginalFreq(ModelZl)$table[-1]
a3 <- table(netherlandsimmigrant$capture)
matplot(y = sqrt(cbind(a1, a2, a3)), x = 1:6, ylab = "Square root of Frequencies",
        xlab = "Counts", pch = 21:23, type = "o",
        lty = 1, col = 1:3, main = "Plot of observed and fitted marginal frequencies")
legend("topright",
       legend = c("Poisson model",
                  "Zelterman model",
                  "Observed data"),
       col = 1:3, pch = 21:23)
```

<img src="man/figures/README-plot-1.png" width="100%" />

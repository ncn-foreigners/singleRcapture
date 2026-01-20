# Regression fitting in single-source capture-recapture models

`estimatePopsizeFit` does for `estimatePopsize` what `glm.fit` does for
`glm`. It is internally called in `estimatePopsize`. Since
`estimatePopsize` does much more than just regression fitting
`estimatePopsizeFit` is much faster.

## Usage

``` r
estimatePopsizeFit(
  y,
  X,
  family,
  control,
  method,
  priorWeights,
  coefStart,
  etaStart,
  offset,
  ...
)
```

## Arguments

- y:

  vector of dependent variables.

- X:

  model matrix, the vglm one.

- family:

  same as model in `estimatePopsize`.

- control:

  control parameters created in `controlModel`.

- method:

  method of estimation same as in `estimatePopsize`.

- priorWeights:

  vector of prior weights its the same argument as weights in
  `estimatePopsize`.

- etaStart, coefStart:

  initial value of regression parameters or linear predictors.

- offset:

  offset passed from by default passed from
  [`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md).

- ...:

  arguments to pass to other methods.

## Value

List with regression parameters, working weights (if IRLS fitting
method) was chosen and number of iterations taken.

## Details

If `method` argument was set to `"optim"` the
[`stats::optim`](https://rdrr.io/r/stats/optim.html) function will be
used to fit regression with analytically computed gradient and (minus)
log likelihood functions as `gr` and `fn` arguments. Unfortunately
`optim` does not allow for hessian to be specified. More information
about how to modify `optim` fitting is included in
[`controlMethod()`](https://ncn-foreigners.github.io/singleRcapture/reference/controlMethod.md).

If `method` argument was set to `"IRLS"` the iteratively reweighted
least squares. The algorithm is well know in generalised linear models.
Thomas W. Yee later extended this algorithm to vector generalised linear
models and in more general terms it can roughly be described as (this is
Yee's description after changing some conventions):

1.  Initialize with:

    - `converged <- FALSE`

    - `iter <- 1`

    - \\\boldsymbol{\beta}\\` <- start`

    - \\\boldsymbol{W}\\` <- prior`

    - \\\ell\\` <- `\\\ell(\boldsymbol{\beta})\\

2.  If `converged` or `iter > Maxiter` move to step 7.

3.  Store values from previous algorithm step:

    - \\\boldsymbol{W}\_{-}\\` <- ` \\\boldsymbol{W}\\

    - \\\ell\_{-}\\` <- ` \\\ell\\

    - \\\boldsymbol{\beta}\_{-}\\` <- ` \\\boldsymbol{\beta}\\

    and assign values at current step:

    - \\\boldsymbol{\eta}\\` <- `
      \\\boldsymbol{X}\_{vlm}\boldsymbol{\beta}\\

    - \\Z\_{i}\\` <- ` \\
      \eta\_{i}+\frac{\partial\ell\_{i}}{\partial\eta\_{i}}
      \mathbb{E}\left(\frac{\partial^{2}\ell\_{i}}{
      \partial\eta\_{i}^{T}\partial\eta\_{i}}\right)^{-1}\\

    - \\\boldsymbol{W}\_{ij}\\` <- `
      \\\mathbb{E}\left(\frac{\partial^{2}\ell}{
      \partial\boldsymbol{\eta}\_{j}^{T}\partial\boldsymbol{\eta}\_{i}}\right)\\

    where \\\ell\_{i}\\ is the ith component of log likelihood function,
    \\\eta\_{i}\\ is the vector of linear predictors associated with ith
    row and \\\mathbb{E}\left(\frac{\partial^{2}\ell\_{i}}{
    \partial\eta\_{i}^{T}\partial\eta\_{i}}\right)\\ corresponds to
    weights associated with ith row and \\\boldsymbol{W}\\ is a block
    matrix, made of diagonal matrixes
    \\\mathbb{E}\left(\frac{\partial^{2}\ell}{
    \partial\boldsymbol{\eta}\_{j}^{T}\partial\boldsymbol{\eta}\_{i}}\right)\\

4.  Regress \\\boldsymbol{Z}\\ on \\\boldsymbol{X}\_{vlm}\\ to obtain
    \\\boldsymbol{\beta}\\ as: \\\boldsymbol{\beta}=
    \left(\boldsymbol{X}\_{vlm}^{T}\boldsymbol{W}\boldsymbol{X}\_{vlm}\right)^{-1}
    \boldsymbol{X}\_{vlm}^{T}\boldsymbol{W}\boldsymbol{Z}\\

5.  Assign:

    - `converged <- `\\ \ell(\boldsymbol{\beta})-\ell\_{-} \<
      \varepsilon\cdot\ell\_{-}\\ or \\
      \|\|\boldsymbol{\beta}-\boldsymbol{\beta}\_{-}\|\|\_{\infty} \<
      \varepsilon\\

    - `iter <- iter + 1`

    where \\\varepsilon\\ is the relative tolerance level, by default
    `1e-8`.

6.  Return to step 2.

7.  Return \\\boldsymbol{\beta}, \boldsymbol{W}\\, `iter`.

In this package we use different conventions for
\\\boldsymbol{X}\_{vlm}\\ matrix hence slight differences are present in
algorithm description but results are identical.

## References

Yee, T. W. (2015). Vector Generalized Linear and Additive Models: With
an Implementation in R. New York, USA: Springer. ISBN 978-1-4939-2817-0.

## See also

[`stats::glm()`](https://rdrr.io/r/stats/glm.html)
[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)
[`controlMethod()`](https://ncn-foreigners.github.io/singleRcapture/reference/controlMethod.md)
[`stats::optim()`](https://rdrr.io/r/stats/optim.html)

## Author

Piotr Chlebicki, Maciej Beresewicz

## Examples

``` r
# \donttest{
summary(farmsubmission)
#>    TOTAL_SUB        log_size       log_distance      C_TYPE    
#>  Min.   : 1.00   Min.   : 0.000   Min.   : 4.102   Beef :5336  
#>  1st Qu.: 1.00   1st Qu.: 4.673   1st Qu.:10.351   Dairy:6700  
#>  Median : 1.00   Median : 5.347   Median :10.778               
#>  Mean   : 2.34   Mean   : 5.259   Mean   :10.662               
#>  3rd Qu.: 3.00   3rd Qu.: 5.940   3rd Qu.:11.099               
#>  Max.   :47.00   Max.   :10.480   Max.   :12.097               

# construct vglm model matrix
X <- matrix(data = 0, nrow = 2 * NROW(farmsubmission), ncol = 7)
X[1:NROW(farmsubmission), 1:4] <- model.matrix(
~ 1 + log_size + log_distance + C_TYPE, 
farmsubmission
)


X[-(1:NROW(farmsubmission)), 5:7] <- X[1:NROW(farmsubmission), c(1, 3, 4)]

# this attribute tells the function which elements of the design matrix 
# correspond to which linear predictor 
attr(X, "hwm") <- c(4, 3)

# get starting points
start <- glm.fit(
y = farmsubmission$TOTAL_SUB, 
x = X[1:NROW(farmsubmission), 1:4], 
family = poisson()
)$coefficients

res <- estimatePopsizeFit(
y = farmsubmission$TOTAL_SUB, 
X = X, 
method = "IRLS", 
priorWeights = 1, 
family = ztoigeom(), 
control = controlMethod(verbose = 5), 
coefStart = c(start, 0, 0, 0),
etaStart = matrix(X %*% c(start, 0, 0, 0), ncol = 2),
offset = cbind(rep(0, NROW(farmsubmission)), rep(0, NROW(farmsubmission)))
)
#> Iteration number 1 log-likelihood: -17550.114
#> Parameter vector:  -2.376169746  0.526024958 -0.056369900  0.413872736 -1.952364227  0.029469787 -0.352737911
#> log-likelihood reduction:  Inf
#> Value of gradient at current step:
#>   1489.26417  8462.99644 15802.73557  1127.53995  -418.45144 -4455.57223  -253.63196
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 2.3761697
#> ----
#> Iteration number 2 log-likelihood: -17300.07
#> Parameter vector:  -2.633896686  0.582395291 -0.064073582  0.667740572 -3.095480026  0.047207764  0.046264308
#> log-likelihood reduction:  250.04378
#> Value of gradient at current step:
#>   -28.899411 -315.159403 -298.325366 -121.257418  -45.032685 -482.019954  -20.096746
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 1.1431158
#> ----
#> Iteration number 3 log-likelihood: -17287.188
#> Parameter vector:  -2.504876876  0.575295120 -0.073821708  0.612160832 -3.119640124  0.011888093  0.029514971
#> log-likelihood reduction:  12.882004
#> Value of gradient at current step:
#>   41.5353296 252.3498380 448.9046975  30.7496910  -8.7641682 -93.3640710  -5.3929815
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.12901981
#> ----
#> Iteration number 4 log-likelihood: -17286.557
#> Parameter vector:  -2.608145817  0.582912458 -0.069091913  0.622397984 -3.562325812  0.036715781  0.114980134
#> log-likelihood reduction:  0.63083188
#> Value of gradient at current step:
#>    1.65638471   9.43319152  16.29083146  -0.22561323  -1.22104226 -13.44032132  -0.49242892
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.44268569
#> ----
#> Iteration number 5 log-likelihood: -17286.527
#> Parameter vector:  -2.589615777  0.583968034 -0.071694931  0.624101247 -3.331535027  0.010372891  0.146052220
#> log-likelihood reduction:  0.030118595
#> Value of gradient at current step:
#>   0.0045622452  0.7724124445  0.4383942297  0.1410982952 -0.1743585631 -1.7157841593  0.0436948698
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.23079079
#> ----
#> Iteration number 6 log-likelihood: -17286.524
#> Parameter vector:  -2.603607106  0.584395570 -0.070803343  0.626306026 -3.469650528  0.021013822  0.172842019
#> log-likelihood reduction:  0.0033094037
#> Value of gradient at current step:
#>  -0.229287967 -1.254231939 -2.663094385 -0.136018786 -0.065501279 -0.875136888 -0.037830598
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.1381155
#> ----
#> Iteration number 7 log-likelihood: -17286.523
#> Parameter vector:  -2.595282607  0.584455592 -0.071643641  0.626377284 -3.357657769  0.010243233  0.173254080
#> log-likelihood reduction:  0.00058225262
#> Value of gradient at current step:
#>   0.01697688169  0.18152150469  0.30278294083  0.05491324450 -0.00954208534  0.00046176213  0.02942429986
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.11199276
#> ----
#> Iteration number 8 log-likelihood: -17286.523
#> Parameter vector:  -2.601354233  0.584513757 -0.071144386  0.626915261 -3.431719299  0.016721304  0.179806607
#> log-likelihood reduction:  0.00017410518
#> Value of gradient at current step:
#>  -0.064034760 -0.379810901 -0.778700450 -0.045639365 -0.011572568 -0.203677953 -0.019963628
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.07406153
#> ----
#> Iteration number 9 log-likelihood: -17286.523
#> Parameter vector:  -2.597126682  0.584506154 -0.071529657  0.626744005 -3.377281943  0.011722915  0.177519492
#> log-likelihood reduction:  0.000071507806
#> Value of gradient at current step:
#>  0.0284979151 0.1864404866 0.3666111521 0.0305300981 0.0016229304 0.0719860549 0.0142857723
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.054437356
#> ----
#> Iteration number 10 log-likelihood: -17286.523
#> Parameter vector:  -2.600150970  0.584522770 -0.071267073  0.626938070 -3.415354274  0.015139829  0.179959991
#> log-likelihood reduction:  0.000035001824
#> Value of gradient at current step:
#>  -0.0254323196 -0.1573738913 -0.3179746249 -0.0217129745 -0.0035760543 -0.0777211830 -0.0102459754
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.038072332
#> ----
#> Iteration number 11 log-likelihood: -17286.523
#> Parameter vector:  -2.598008714  0.584514838 -0.071457541  0.626825132 -3.388080609  0.012664531  0.178504506
#> log-likelihood reduction:  0.000017260521
#> Value of gradient at current step:
#>  0.0166844643 0.1053460426 0.2101040472 0.0156399950 0.0016132423 0.0451631907 0.0071883615
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.027273665
#> ----
#> Iteration number 12 log-likelihood: -17286.523
#> Parameter vector:  -2.599534056  0.584521813 -0.071323467  0.626914011 -3.407406275  0.014409036  0.179635638
#> log-likelihood reduction:  0.0000087810295
#> Value of gradient at current step:
#>  -0.0121554909 -0.0761209122 -0.1529982539 -0.0109236603 -0.0015139582 -0.0360631619 -0.0051776855
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.019325666
#> ----
#> Iteration number 13 log-likelihood: -17286.523
#> Parameter vector:  -2.598450181  0.584517321 -0.071419274  0.626853774 -3.393643298  0.013163381  0.178864482
#> log-likelihood reduction:  0.0000044078333
#> Value of gradient at current step:
#>  0.00863941235 0.05422850142 0.10851025481 0.00789351520 0.00092550566 0.02403409764 0.00364208947
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.013762977
#> ----
#> Iteration number 14 log-likelihood: -17286.523
#> Parameter vector:  -2.599221083  0.584520680 -0.071351319  0.626897628 -3.403423424  0.014047410  0.179424307
#> log-likelihood reduction:  0.0000022386557
#> Value of gradient at current step:
#>  -0.00609651399 -0.03826378502 -0.07677634455 -0.00553869224 -0.00072252915 -0.01776979990 -0.00261189222
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.0097801258
#> ----
#> Iteration number 15 log-likelihood: -17286.523
#> Parameter vector:  -2.598672961  0.584518350 -0.071399702  0.626866800 -3.396467246  0.013418234  0.179030220
#> log-likelihood reduction:  0.0000011281627
#> Value of gradient at current step:
#>  0.00437552328 0.02744244081 0.05496756154 0.00397917475 0.00048523685 0.01233742716 0.00184538567
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.0069561777
#> ----
#> Iteration number 16 log-likelihood: -17286.523
#> Parameter vector:  -2.599062724  0.584520028 -0.071365321  0.626888846 -3.401413268  0.013865450  0.179311848
#> log-likelihood reduction:  0.00000057176658
#> Value of gradient at current step:
#>  -0.00308520077 -0.01936683164 -0.03883603649 -0.00280772877 -0.00035808186 -0.00890811713 -0.00131846263
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.0049460217
#> ----
#> Iteration number 17 log-likelihood: -17286.523
#> Parameter vector:  -2.598785568  0.584518843 -0.071389778  0.626873213 -3.397896233  0.013547389  0.179112085
#> log-likelihood reduction:  0.00000028862269
#> Value of gradient at current step:
#>  0.00220900398 0.01385546526 0.02776312256 0.00200791788 0.00024850391 0.00627132383 0.00093429219
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.0035170353
#> ----
#> Iteration number 18 log-likelihood: -17286.523
#> Parameter vector:  -2.598982644  0.584519689 -0.071372391  0.626884345 -3.400397179  0.013773541  0.179254311
#> log-likelihood reduction:  0.00000014611214
#> Value of gradient at current step:
#>  -0.00156243814 -0.00980648384 -0.01966002003 -0.00142189173 -0.00017963917 -0.00448948244 -0.00066602383
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.0025009461
#> ----
#> Iteration number 19 log-likelihood: -17286.523
#> Parameter vector:  -2.598842504  0.584519089 -0.071384756  0.626876435 -3.398618850  0.013612726  0.179153243
#> log-likelihood reduction:  0.00000007381459
#> Value of gradient at current step:
#>  0.00111546600 0.00699753811 0.01402367467 0.00101410800 0.00012631018 0.00317780227 0.00047274423
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.0017783286
#> ----
#> Iteration number 20 log-likelihood: -17286.523
#> Parameter vector:  -2.598942153  0.584519516 -0.071375964  0.626882062 -3.399883416  0.013727078  0.179225136
#> log-likelihood reduction:  0.00000003734749
#> Value of gradient at current step:
#>  -0.000790860150 -0.004963117716 -0.009948972121 -0.000719565560 -0.000090519732 -0.002266834445 -0.000336593619
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.0012645652
#> ----
#> Iteration number 21 log-likelihood: -17286.523
#> Parameter vector:  -2.598871294  0.584519212 -0.071382216  0.626878062 -3.398984225  0.013645765  0.179174024
#> log-likelihood reduction:  0.000000018873834
#> Value of gradient at current step:
#>  0.000563568332 0.003535731060 0.007086435089 0.000512456788 0.000064017826 0.001608354191 0.000239124512
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.0008991904
#> ----
#> Iteration number 22 log-likelihood: -17286.523
#> Parameter vector:  -2.598921680  0.584519428 -0.071377771  0.626880906 -3.399623632  0.013703585  0.179210372
#> log-likelihood reduction:  0.0000000095496944
#> Value of gradient at current step:
#>  -0.000400125264 -0.002510830909 -0.005032906088 -0.000363998081 -0.000045696107 -0.001145438950 -0.000170148829
#> Algorithm will terminate if one of following conditions will be met:
#> The increase to minus log-likelihood will be bellow chosen value of epsilon 1e-08 
#> Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.
#> At current step the highest change was: 0.00063940667
#> ----
#> Value of analytically computed hessian at fitted regression coefficients:
#>            [,1]         [,2]         [,3]        [,4]       [,5]       [,6]
#> [1,]  -5941.804  -33237.6696  -63184.3939  -4042.8998  424.63598 23751.3293
#> [2,] -33237.670 -191051.4616 -353329.4001 -23640.8438 4526.54681  1389.4971
#> [3,] -63184.394 -353329.4001 -674494.4908 -42806.7283  247.79024  4526.5468
#> [4,]  -4042.900  -23640.8438  -42806.7283  -4042.8998 2229.15341 48421.1063
#> [5,]    424.636     247.7902   48421.1063   1389.4971  -69.94232  -744.9399
#> [6,]   2229.153    4526.5468    2624.5806   2624.5806 -744.93988 -7973.6057
#> [7,]   4526.547   23751.3293     247.7902    247.7902  -45.99667  -490.0571
#>            [,7]
#> [1,] 2624.58059
#> [2,]  247.79024
#> [3,] 2624.58059
#> [4,]  247.79024
#> [5,]  -45.99667
#> [6,] -490.05713
#> [7,]  -45.99667
#> The matrix above has the following eigen values:
#>  -869785.7+0i -17240.18+0i -2472.809+9708.042i -2472.809-9708.042i 8149.489+0i 1615.552+0i -1413.766+0i 

# extract results

# regression coefficient vector
res$beta
#> [1] -2.59892168  0.58451943 -0.07137777  0.62688091 -3.39962363  0.01370359
#> [7]  0.17921037

# check likelihood
ll <- ztoigeom()$makeMinusLogLike(y = farmsubmission$TOTAL_SUB, X = X)

-ll(res$beta)
#> [1] -17286.52

# number of iterations
res$iter
#> [1] 22

# working weights
head(res$weights)
#>         lambda       mixed       mixed       omega
#> [1,] 0.2733879 -0.03187827 -0.03187827 0.004860059
#> [2,] 0.6783958 -0.03738217 -0.03738217 0.007747769
#> [3,] 0.1145154 -0.02534924 -0.02534924 0.006028454
#> [4,] 0.1790770 -0.02931302 -0.02931302 0.005523789
#> [5,] 0.3987180 -0.03319312 -0.03319312 0.004479829
#> [6,] 0.4920636 -0.03350614 -0.03350614 0.004531546

# Compare with optim call

res2 <- estimatePopsizeFit(
  y = farmsubmission$TOTAL_SUB, 
  X = X, 
  method = "optim", 
  priorWeights = 1, 
  family = ztoigeom(), 
  coefStart = c(start, 0, 0, 0),
  control = controlMethod(verbose = 1, silent = TRUE),
  offset = cbind(rep(0, NROW(farmsubmission)), rep(0, NROW(farmsubmission)))
)
#>   Nelder-Mead direct search function minimizer
#> function value for initial parameters = 19135.970356
#>   Scaled convergence tolerance is 0.00019136
#> Stepsize computed as 0.082584
#> BUILD              8 21628.437230 19133.539224
#> REFLECTION        10 20756.644804 18940.732447
#> LO-REDUCTION      12 19606.916214 18926.129600
#> LO-REDUCTION      14 19309.305936 18847.773068
#> EXTENSION         16 19269.501942 18303.096276
#> LO-REDUCTION      18 19153.670565 18303.096276
#> EXTENSION         20 19135.970356 17912.380725
#> LO-REDUCTION      22 19133.539224 17912.380725
#> LO-REDUCTION      24 18940.732447 17912.380725
#> LO-REDUCTION      26 18926.129600 17912.380725
#> LO-REDUCTION      28 18847.773068 17912.380725
#> LO-REDUCTION      30 18330.888204 17912.380725
#> EXTENSION         32 18303.096276 17699.648860
#> EXTENSION         34 18212.293899 17558.677603
#> HI-REDUCTION      36 18078.266152 17558.677603
#> LO-REDUCTION      38 18077.244658 17558.677603
#> EXTENSION         40 18046.635595 17434.246747
#> LO-REDUCTION      42 17966.526721 17434.246747
#> LO-REDUCTION      44 17912.380725 17434.246747
#> LO-REDUCTION      46 17871.858218 17434.246747
#> LO-REDUCTION      48 17699.648860 17434.246747
#> LO-REDUCTION      50 17650.163809 17434.246747
#> LO-REDUCTION      52 17558.677603 17434.246747
#> HI-REDUCTION      54 17516.785947 17434.246747
#> REFLECTION        56 17467.768303 17420.151936
#> HI-REDUCTION      58 17467.699076 17420.151936
#> EXTENSION         60 17465.840145 17381.144605
#> LO-REDUCTION      62 17465.538855 17381.144605
#> LO-REDUCTION      64 17457.288417 17381.144605
#> LO-REDUCTION      66 17442.024150 17381.144605
#> LO-REDUCTION      68 17440.392477 17381.144605
#> EXTENSION         70 17434.246747 17334.895518
#> LO-REDUCTION      72 17420.151936 17334.895518
#> LO-REDUCTION      74 17416.961045 17334.895518
#> LO-REDUCTION      76 17402.660537 17334.895518
#> LO-REDUCTION      78 17394.959214 17334.895518
#> LO-REDUCTION      80 17392.532797 17334.895518
#> EXTENSION         82 17381.144605 17316.504861
#> LO-REDUCTION      84 17375.832799 17316.504861
#> LO-REDUCTION      86 17363.774766 17316.504861
#> REFLECTION        88 17343.751331 17305.644043
#> LO-REDUCTION      90 17342.366351 17305.644043
#> EXTENSION         92 17338.983203 17299.308930
#> LO-REDUCTION      94 17334.895518 17299.308930
#> LO-REDUCTION      96 17324.646822 17299.308930
#> LO-REDUCTION      98 17323.590625 17299.308930
#> LO-REDUCTION     100 17316.627163 17299.308930
#> LO-REDUCTION     102 17316.504861 17299.308930
#> LO-REDUCTION     104 17309.426844 17299.308930
#> LO-REDUCTION     106 17306.896574 17299.213022
#> LO-REDUCTION     108 17305.644043 17299.213022
#> HI-REDUCTION     110 17301.855713 17299.213022
#> REFLECTION       112 17301.803640 17298.374210
#> HI-REDUCTION     114 17301.796394 17298.374210
#> LO-REDUCTION     116 17301.489802 17298.374210
#> LO-REDUCTION     118 17300.957660 17298.374210
#> LO-REDUCTION     120 17300.016126 17298.374210
#> HI-REDUCTION     122 17299.611047 17298.374210
#> REFLECTION       124 17299.308930 17298.093107
#> EXTENSION        126 17299.213022 17297.231097
#> REFLECTION       128 17298.925386 17297.171892
#> REFLECTION       130 17298.732554 17297.077851
#> LO-REDUCTION     132 17298.452859 17297.077851
#> LO-REDUCTION     134 17298.382263 17297.077851
#> REFLECTION       136 17298.374210 17296.980993
#> REFLECTION       138 17298.093107 17296.468183
#> EXTENSION        140 17297.268148 17295.433620
#> LO-REDUCTION     142 17297.231097 17295.433620
#> LO-REDUCTION     144 17297.201571 17295.433620
#> REFLECTION       146 17297.171892 17295.341847
#> LO-REDUCTION     148 17297.077851 17295.341847
#> EXTENSION        150 17296.980993 17294.284096
#> LO-REDUCTION     152 17296.468183 17294.284096
#> EXTENSION        154 17295.993173 17293.316067
#> LO-REDUCTION     156 17295.517735 17293.316067
#> LO-REDUCTION     158 17295.494544 17293.316067
#> LO-REDUCTION     160 17295.433620 17293.316067
#> EXTENSION        162 17295.341847 17292.072213
#> LO-REDUCTION     164 17295.283117 17292.072213
#> LO-REDUCTION     166 17294.357739 17292.072213
#> LO-REDUCTION     168 17294.284096 17292.072213
#> LO-REDUCTION     170 17294.231320 17292.072213
#> EXTENSION        172 17294.104223 17291.423540
#> LO-REDUCTION     174 17293.669383 17291.423540
#> REFLECTION       176 17293.316067 17291.060821
#> EXTENSION        178 17292.290481 17290.565091
#> LO-REDUCTION     180 17292.217296 17290.565091
#> LO-REDUCTION     182 17292.188283 17290.565091
#> LO-REDUCTION     184 17292.072213 17290.565091
#> LO-REDUCTION     186 17291.567818 17290.565091
#> LO-REDUCTION     188 17291.423540 17290.517025
#> EXTENSION        190 17291.060821 17289.934833
#> LO-REDUCTION     192 17291.007425 17289.934833
#> LO-REDUCTION     194 17290.866709 17289.934833
#> LO-REDUCTION     196 17290.817445 17289.934833
#> LO-REDUCTION     198 17290.758168 17289.934833
#> LO-REDUCTION     200 17290.565091 17289.934833
#> LO-REDUCTION     202 17290.530576 17289.934833
#> EXTENSION        204 17290.517025 17289.755906
#> LO-REDUCTION     206 17290.282088 17289.755906
#> LO-REDUCTION     208 17290.248001 17289.755906
#> LO-REDUCTION     210 17290.223189 17289.755906
#> REFLECTION       212 17290.127309 17289.725940
#> LO-REDUCTION     214 17290.098882 17289.725940
#> LO-REDUCTION     216 17290.036552 17289.725940
#> LO-REDUCTION     218 17289.934833 17289.725940
#> LO-REDUCTION     220 17289.902270 17289.725940
#> LO-REDUCTION     222 17289.890661 17289.725940
#> REFLECTION       224 17289.820937 17289.702686
#> LO-REDUCTION     226 17289.777081 17289.702312
#> REFLECTION       228 17289.758951 17289.686667
#> EXTENSION        230 17289.755906 17289.666388
#> EXTENSION        232 17289.747982 17289.610366
#> EXTENSION        234 17289.745697 17289.562699
#> HI-REDUCTION     236 17289.725940 17289.562699
#> LO-REDUCTION     238 17289.702686 17289.562699
#> LO-REDUCTION     240 17289.702312 17289.562699
#> LO-REDUCTION     242 17289.686667 17289.562699
#> REFLECTION       244 17289.666388 17289.553928
#> EXTENSION        246 17289.659099 17289.430319
#> LO-REDUCTION     248 17289.625783 17289.430319
#> LO-REDUCTION     250 17289.610366 17289.430319
#> LO-REDUCTION     252 17289.589639 17289.430319
#> EXTENSION        254 17289.581544 17289.366344
#> LO-REDUCTION     256 17289.562699 17289.366344
#> EXTENSION        258 17289.553928 17289.355186
#> REFLECTION       260 17289.524309 17289.347390
#> EXTENSION        262 17289.463415 17289.256311
#> LO-REDUCTION     264 17289.438892 17289.256311
#> EXTENSION        266 17289.430319 17289.194525
#> LO-REDUCTION     268 17289.422342 17289.194525
#> LO-REDUCTION     270 17289.366344 17289.194525
#> REFLECTION       272 17289.355186 17289.169478
#> HI-REDUCTION     274 17289.347390 17289.169478
#> LO-REDUCTION     276 17289.258763 17289.169478
#> LO-REDUCTION     278 17289.256311 17289.169478
#> LO-REDUCTION     280 17289.233544 17289.169478
#> HI-REDUCTION     282 17289.220868 17289.169478
#> LO-REDUCTION     284 17289.211407 17289.169478
#> LO-REDUCTION     286 17289.194525 17289.169478
#> LO-REDUCTION     288 17289.192199 17289.169478
#> LO-REDUCTION     290 17289.186075 17289.169478
#> REFLECTION       292 17289.177779 17289.164394
#> LO-REDUCTION     294 17289.177465 17289.164394
#> LO-REDUCTION     296 17289.176353 17289.164394
#> LO-REDUCTION     298 17289.173142 17289.164394
#> LO-REDUCTION     300 17289.171672 17289.164394
#> REFLECTION       302 17289.170351 17289.162643
#> LO-REDUCTION     304 17289.169550 17289.162479
#> REFLECTION       306 17289.169478 17289.160797
#> REFLECTION       308 17289.165906 17289.159411
#> LO-REDUCTION     310 17289.165109 17289.159411
#> LO-REDUCTION     312 17289.164417 17289.159411
#> LO-REDUCTION     314 17289.164394 17289.159411
#> REFLECTION       316 17289.162643 17289.158035
#> LO-REDUCTION     318 17289.162479 17289.158035
#> LO-REDUCTION     320 17289.161336 17289.158035
#> LO-REDUCTION     322 17289.160972 17289.158035
#> REFLECTION       324 17289.160797 17289.157923
#> LO-REDUCTION     326 17289.160192 17289.157923
#> REFLECTION       328 17289.159411 17289.157174
#> LO-REDUCTION     330 17289.159135 17289.157174
#> LO-REDUCTION     332 17289.159042 17289.157174
#> REFLECTION       334 17289.158322 17289.156947
#> REFLECTION       336 17289.158042 17289.156521
#> LO-REDUCTION     338 17289.158035 17289.156521
#> LO-REDUCTION     340 17289.157923 17289.156521
#> REFLECTION       342 17289.157547 17289.156342
#> LO-REDUCTION     344 17289.157424 17289.156342
#> REFLECTION       346 17289.157174 17289.156151
#> LO-REDUCTION     348 17289.157122 17289.156151
#> LO-REDUCTION     350 17289.157110 17289.156151
#> REFLECTION       352 17289.156947 17289.155928
#> LO-REDUCTION     354 17289.156532 17289.155928
#> EXTENSION        356 17289.156521 17289.155659
#> EXTENSION        358 17289.156501 17289.155347
#> EXTENSION        360 17289.156342 17289.155165
#> LO-REDUCTION     362 17289.156239 17289.155165
#> EXTENSION        364 17289.156151 17289.154673
#> LO-REDUCTION     366 17289.156024 17289.154673
#> EXTENSION        368 17289.155928 17289.153645
#> LO-REDUCTION     370 17289.155659 17289.153645
#> LO-REDUCTION     372 17289.155347 17289.153645
#> LO-REDUCTION     374 17289.155336 17289.153645
#> LO-REDUCTION     376 17289.155165 17289.153645
#> REFLECTION       378 17289.155017 17289.153558
#> EXTENSION        380 17289.154673 17289.153090
#> LO-REDUCTION     382 17289.154524 17289.153090
#> EXTENSION        384 17289.154080 17289.152577
#> EXTENSION        386 17289.154065 17289.152162
#> REFLECTION       388 17289.153701 17289.152077
#> LO-REDUCTION     390 17289.153645 17289.152077
#> REFLECTION       392 17289.153558 17289.151837
#> LO-REDUCTION     394 17289.153156 17289.151837
#> LO-REDUCTION     396 17289.153090 17289.151837
#> LO-REDUCTION     398 17289.152577 17289.151837
#> HI-REDUCTION     400 17289.152418 17289.151837
#> EXTENSION        402 17289.152162 17289.150909
#> LO-REDUCTION     404 17289.152109 17289.150909
#> LO-REDUCTION     406 17289.152077 17289.150909
#> LO-REDUCTION     408 17289.152055 17289.150909
#> LO-REDUCTION     410 17289.151913 17289.150909
#> LO-REDUCTION     412 17289.151868 17289.150909
#> LO-REDUCTION     414 17289.151837 17289.150909
#> EXTENSION        416 17289.151677 17289.150529
#> LO-REDUCTION     418 17289.151620 17289.150529
#> EXTENSION        420 17289.151536 17289.150258
#> LO-REDUCTION     422 17289.151371 17289.150258
#> EXTENSION        424 17289.151226 17289.149282
#> LO-REDUCTION     426 17289.151019 17289.149282
#> LO-REDUCTION     428 17289.150909 17289.149282
#> LO-REDUCTION     430 17289.150778 17289.149282
#> LO-REDUCTION     432 17289.150529 17289.149282
#> LO-REDUCTION     434 17289.150459 17289.149282
#> LO-REDUCTION     436 17289.150258 17289.149282
#> EXTENSION        438 17289.150196 17289.149030
#> LO-REDUCTION     440 17289.150015 17289.149030
#> EXTENSION        442 17289.149732 17289.148215
#> LO-REDUCTION     444 17289.149724 17289.148215
#> LO-REDUCTION     446 17289.149661 17289.148215
#> EXTENSION        448 17289.149490 17289.148082
#> REFLECTION       450 17289.149282 17289.147882
#> HI-REDUCTION     452 17289.149154 17289.147882
#> EXTENSION        454 17289.149030 17289.147689
#> EXTENSION        456 17289.148634 17289.146503
#> LO-REDUCTION     458 17289.148548 17289.146503
#> LO-REDUCTION     460 17289.148356 17289.146503
#> LO-REDUCTION     462 17289.148215 17289.146503
#> LO-REDUCTION     464 17289.148082 17289.146503
#> LO-REDUCTION     466 17289.147882 17289.146503
#> LO-REDUCTION     468 17289.147689 17289.146503
#> LO-REDUCTION     470 17289.147390 17289.146503
#> EXTENSION        472 17289.147372 17289.146209
#> EXTENSION        474 17289.147325 17289.145509
#> LO-REDUCTION     476 17289.147206 17289.145509
#> LO-REDUCTION     478 17289.147159 17289.145509
#> EXTENSION        480 17289.146556 17289.144353
#> LO-REDUCTION     482 17289.146550 17289.144353
#> LO-REDUCTION     484 17289.146503 17289.144353
#> EXTENSION        486 17289.146209 17289.143600
#> LO-REDUCTION     488 17289.145788 17289.143600
#> LO-REDUCTION     490 17289.145629 17289.143600
#> EXTENSION        492 17289.145509 17289.142801
#> LO-REDUCTION     494 17289.145288 17289.142801
#> EXTENSION        496 17289.144749 17289.142580
#> LO-REDUCTION     498 17289.144353 17289.142580
#> REFLECTION       500 17289.143924 17289.142388
#> LO-REDUCTION     502 17289.143813 17289.142388
#> LO-REDUCTION     504 17289.143710 17289.142388
#> EXTENSION        506 17289.143600 17289.141784
#> LO-REDUCTION     508 17289.142801 17289.141784
#> HI-REDUCTION     510 17289.142748 17289.141784
#> LO-REDUCTION     512 17289.142616 17289.141784
#> LO-REDUCTION     514 17289.142580 17289.141784
#> REFLECTION       516 17289.142500 17289.141716
#> LO-REDUCTION     518 17289.142388 17289.141716
#> REFLECTION       520 17289.142334 17289.141585
#> LO-REDUCTION     522 17289.142234 17289.141585
#> EXTENSION        524 17289.141991 17289.140765
#> LO-REDUCTION     526 17289.141796 17289.140765
#> LO-REDUCTION     528 17289.141786 17289.140765
#> LO-REDUCTION     530 17289.141784 17289.140765
#> LO-REDUCTION     532 17289.141772 17289.140765
#> EXTENSION        534 17289.141716 17289.140229
#> LO-REDUCTION     536 17289.141585 17289.140229
#> EXTENSION        538 17289.141362 17289.139857
#> LO-REDUCTION     540 17289.141111 17289.139857
#> EXTENSION        542 17289.140987 17289.139368
#> LO-REDUCTION     544 17289.140866 17289.139368
#> LO-REDUCTION     546 17289.140765 17289.139368
#> REFLECTION       548 17289.140709 17289.139248
#> LO-REDUCTION     550 17289.140229 17289.139248
#> REFLECTION       552 17289.139904 17289.139012
#> REFLECTION       554 17289.139857 17289.138897
#> LO-REDUCTION     556 17289.139487 17289.138897
#> LO-REDUCTION     558 17289.139371 17289.138897
#> LO-REDUCTION     560 17289.139368 17289.138897
#> LO-REDUCTION     562 17289.139253 17289.138897
#> EXTENSION        564 17289.139248 17289.138535
#> LO-REDUCTION     566 17289.139012 17289.138535
#> LO-REDUCTION     568 17289.139007 17289.138535
#> LO-REDUCTION     570 17289.138999 17289.138535
#> LO-REDUCTION     572 17289.138990 17289.138535
#> LO-REDUCTION     574 17289.138975 17289.138535
#> LO-REDUCTION     576 17289.138897 17289.138535
#> REFLECTION       578 17289.138817 17289.138485
#> LO-REDUCTION     580 17289.138745 17289.138485
#> EXTENSION        582 17289.138723 17289.138282
#> LO-REDUCTION     584 17289.138641 17289.138282
#> LO-REDUCTION     586 17289.138571 17289.138282
#> LO-REDUCTION     588 17289.138543 17289.138282
#> LO-REDUCTION     590 17289.138535 17289.138282
#> LO-REDUCTION     592 17289.138526 17289.138282
#> EXTENSION        594 17289.138485 17289.138060
#> EXTENSION        596 17289.138356 17289.137752
#> LO-REDUCTION     598 17289.138332 17289.137752
#> LO-REDUCTION     600 17289.138307 17289.137752
#> LO-REDUCTION     602 17289.138302 17289.137752
#> EXTENSION        604 17289.138289 17289.137475
#> LO-REDUCTION     606 17289.138282 17289.137475
#> LO-REDUCTION     608 17289.138060 17289.137475
#> EXTENSION        610 17289.138056 17289.137172
#> LO-REDUCTION     612 17289.137937 17289.137172
#> EXTENSION        614 17289.137833 17289.136763
#> LO-REDUCTION     616 17289.137808 17289.136763
#> EXTENSION        618 17289.137752 17289.136209
#> LO-REDUCTION     620 17289.137512 17289.136209
#> EXTENSION        622 17289.137475 17289.135715
#> EXTENSION        624 17289.137183 17289.134720
#> LO-REDUCTION     626 17289.137172 17289.134720
#> LO-REDUCTION     628 17289.136910 17289.134720
#> EXTENSION        630 17289.136763 17289.133898
#> LO-REDUCTION     632 17289.136395 17289.133898
#> EXTENSION        634 17289.136209 17289.132040
#> LO-REDUCTION     636 17289.135715 17289.132040
#> LO-REDUCTION     638 17289.135220 17289.132040
#> EXTENSION        640 17289.134996 17289.130704
#> LO-REDUCTION     642 17289.134720 17289.130704
#> EXTENSION        644 17289.133962 17289.130357
#> EXTENSION        646 17289.133898 17289.128792
#> EXTENSION        648 17289.132444 17289.125823
#> LO-REDUCTION     650 17289.132232 17289.125823
#> LO-REDUCTION     652 17289.132040 17289.125823
#> LO-REDUCTION     654 17289.130832 17289.125823
#> LO-REDUCTION     656 17289.130704 17289.125823
#> LO-REDUCTION     658 17289.130357 17289.125823
#> LO-REDUCTION     660 17289.128792 17289.125823
#> LO-REDUCTION     662 17289.128276 17289.125823
#> LO-REDUCTION     664 17289.128208 17289.125823
#> EXTENSION        666 17289.128166 17289.124162
#> LO-REDUCTION     668 17289.128165 17289.124162
#> EXTENSION        670 17289.127639 17289.122521
#> LO-REDUCTION     672 17289.126750 17289.122521
#> LO-REDUCTION     674 17289.126352 17289.122521
#> EXTENSION        676 17289.125995 17289.120555
#> LO-REDUCTION     678 17289.125823 17289.120555
#> EXTENSION        680 17289.125122 17289.120009
#> EXTENSION        682 17289.124162 17289.118228
#> LO-REDUCTION     684 17289.122960 17289.118228
#> LO-REDUCTION     686 17289.122937 17289.118228
#> REFLECTION       688 17289.122521 17289.118049
#> LO-REDUCTION     690 17289.121998 17289.118049
#> EXTENSION        692 17289.120555 17289.114665
#> HI-REDUCTION     694 17289.120009 17289.114665
#> LO-REDUCTION     696 17289.119590 17289.114665
#> LO-REDUCTION     698 17289.119356 17289.114665
#> LO-REDUCTION     700 17289.118291 17289.114665
#> LO-REDUCTION     702 17289.118288 17289.114665
#> EXTENSION        704 17289.118228 17289.113381
#> LO-REDUCTION     706 17289.118049 17289.113381
#> LO-REDUCTION     708 17289.117504 17289.113381
#> EXTENSION        710 17289.116483 17289.112471
#> EXTENSION        712 17289.115678 17289.109345
#> LO-REDUCTION     714 17289.115541 17289.109345
#> LO-REDUCTION     716 17289.114665 17289.109345
#> LO-REDUCTION     718 17289.114477 17289.109345
#> LO-REDUCTION     720 17289.113634 17289.109345
#> LO-REDUCTION     722 17289.113381 17289.109345
#> LO-REDUCTION     724 17289.112766 17289.109345
#> LO-REDUCTION     726 17289.112471 17289.109345
#> EXTENSION        728 17289.111818 17289.107773
#> LO-REDUCTION     730 17289.111746 17289.107773
#> LO-REDUCTION     732 17289.111245 17289.107773
#> EXTENSION        734 17289.111091 17289.106063
#> LO-REDUCTION     736 17289.110663 17289.106063
#> LO-REDUCTION     738 17289.109600 17289.106063
#> EXTENSION        740 17289.109345 17289.105024
#> LO-REDUCTION     742 17289.108880 17289.105024
#> EXTENSION        744 17289.107989 17289.103245
#> LO-REDUCTION     746 17289.107773 17289.103245
#> EXTENSION        748 17289.106271 17289.102524
#> EXTENSION        750 17289.106140 17289.100785
#> HI-REDUCTION     752 17289.106063 17289.100785
#> LO-REDUCTION     754 17289.105467 17289.100785
#> LO-REDUCTION     756 17289.105024 17289.100785
#> EXTENSION        758 17289.104800 17289.099737
#> LO-REDUCTION     760 17289.103902 17289.099737
#> HI-REDUCTION     762 17289.103245 17289.099737
#> EXTENSION        764 17289.102760 17289.098971
#> EXTENSION        766 17289.102723 17289.097388
#> LO-REDUCTION     768 17289.102524 17289.097388
#> EXTENSION        770 17289.101955 17289.095649
#> LO-REDUCTION     772 17289.100785 17289.095649
#> LO-REDUCTION     774 17289.100237 17289.095649
#> EXTENSION        776 17289.099737 17289.094176
#> LO-REDUCTION     778 17289.098971 17289.094176
#> LO-REDUCTION     780 17289.098143 17289.094176
#> LO-REDUCTION     782 17289.097880 17289.094176
#> EXTENSION        784 17289.097591 17289.091226
#> EXTENSION        786 17289.097388 17289.089069
#> LO-REDUCTION     788 17289.095649 17289.089069
#> LO-REDUCTION     790 17289.094674 17289.089069
#> LO-REDUCTION     792 17289.094652 17289.089069
#> LO-REDUCTION     794 17289.094566 17289.089069
#> EXTENSION        796 17289.094176 17289.087095
#> LO-REDUCTION     798 17289.092354 17289.087095
#> LO-REDUCTION     800 17289.092156 17289.087095
#> LO-REDUCTION     802 17289.091715 17289.087095
#> EXTENSION        804 17289.091226 17289.083939
#> LO-REDUCTION     806 17289.089665 17289.083939
#> LO-REDUCTION     808 17289.089506 17289.083939
#> EXTENSION        810 17289.089069 17289.080477
#> LO-REDUCTION     812 17289.087649 17289.080477
#> LO-REDUCTION     814 17289.087188 17289.080477
#> LO-REDUCTION     816 17289.087095 17289.080477
#> LO-REDUCTION     818 17289.085216 17289.080477
#> LO-REDUCTION     820 17289.084941 17289.080477
#> LO-REDUCTION     822 17289.084029 17289.080477
#> EXTENSION        824 17289.083939 17289.078527
#> LO-REDUCTION     826 17289.082956 17289.078527
#> LO-REDUCTION     828 17289.082866 17289.078527
#> EXTENSION        830 17289.081381 17289.076568
#> LO-REDUCTION     832 17289.081227 17289.076568
#> LO-REDUCTION     834 17289.080496 17289.076568
#> EXTENSION        836 17289.080477 17289.074508
#> HI-REDUCTION     838 17289.079939 17289.074508
#> LO-REDUCTION     840 17289.079200 17289.074508
#> EXTENSION        842 17289.078527 17289.073373
#> EXTENSION        844 17289.077815 17289.072709
#> LO-REDUCTION     846 17289.077222 17289.072709
#> LO-REDUCTION     848 17289.076735 17289.072709
#> LO-REDUCTION     850 17289.076568 17289.072709
#> EXTENSION        852 17289.074832 17289.070756
#> LO-REDUCTION     854 17289.074508 17289.070756
#> LO-REDUCTION     856 17289.074125 17289.070756
#> LO-REDUCTION     858 17289.073623 17289.070756
#> EXTENSION        860 17289.073373 17289.068622
#> LO-REDUCTION     862 17289.073250 17289.068622
#> LO-REDUCTION     864 17289.072709 17289.068622
#> LO-REDUCTION     866 17289.072460 17289.068622
#> LO-REDUCTION     868 17289.072372 17289.068622
#> LO-REDUCTION     870 17289.071455 17289.068622
#> LO-REDUCTION     872 17289.071020 17289.068622
#> LO-REDUCTION     874 17289.070882 17289.068622
#> EXTENSION        876 17289.070756 17289.068081
#> EXTENSION        878 17289.070498 17289.065876
#> LO-REDUCTION     880 17289.069611 17289.065876
#> LO-REDUCTION     882 17289.069527 17289.065876
#> EXTENSION        884 17289.069417 17289.064672
#> EXTENSION        886 17289.068794 17289.061734
#> LO-REDUCTION     888 17289.068622 17289.061734
#> LO-REDUCTION     890 17289.068081 17289.061734
#> EXTENSION        892 17289.066584 17289.057845
#> LO-REDUCTION     894 17289.066230 17289.057845
#> LO-REDUCTION     896 17289.065876 17289.057845
#> EXTENSION        898 17289.064672 17289.055947
#> EXTENSION        900 17289.063358 17289.053844
#> EXTENSION        902 17289.063141 17289.052345
#> EXTENSION        904 17289.061734 17289.050642
#> EXTENSION        906 17289.059974 17289.047222
#> EXTENSION        908 17289.058292 17289.043892
#> LO-REDUCTION     910 17289.057845 17289.043892
#> LO-REDUCTION     912 17289.055947 17289.043892
#> LO-REDUCTION     914 17289.053844 17289.043892
#> HI-REDUCTION     916 17289.052345 17289.043892
#> HI-REDUCTION     918 17289.050642 17289.043892
#> LO-REDUCTION     920 17289.048103 17289.043892
#> HI-REDUCTION     922 17289.047899 17289.043892
#> LO-REDUCTION     924 17289.047866 17289.043892
#> LO-REDUCTION     926 17289.047222 17289.043892
#> LO-REDUCTION     928 17289.047131 17289.043892
#> EXTENSION        930 17289.046353 17289.042209
#> LO-REDUCTION     932 17289.046255 17289.042209
#> LO-REDUCTION     934 17289.045471 17289.042209
#> EXTENSION        936 17289.044714 17289.040269
#> EXTENSION        938 17289.044649 17289.037362
#> LO-REDUCTION     940 17289.044112 17289.037362
#> LO-REDUCTION     942 17289.043892 17289.037362
#> LO-REDUCTION     944 17289.042584 17289.037362
#> LO-REDUCTION     946 17289.042324 17289.037362
#> LO-REDUCTION     948 17289.042209 17289.037362
#> EXTENSION        950 17289.040382 17289.033960
#> LO-REDUCTION     952 17289.040269 17289.033960
#> LO-REDUCTION     954 17289.039777 17289.033960
#> EXTENSION        956 17289.039477 17289.031085
#> LO-REDUCTION     958 17289.039384 17289.031085
#> EXTENSION        960 17289.037477 17289.029757
#> LO-REDUCTION     962 17289.037362 17289.029757
#> LO-REDUCTION     964 17289.036963 17289.029757
#> EXTENSION        966 17289.034508 17289.025230
#> LO-REDUCTION     968 17289.033960 17289.025230
#> LO-REDUCTION     970 17289.033351 17289.025230
#> EXTENSION        972 17289.031085 17289.023573
#> LO-REDUCTION     974 17289.030785 17289.023573
#> LO-REDUCTION     976 17289.030431 17289.023573
#> LO-REDUCTION     978 17289.029757 17289.023573
#> LO-REDUCTION     980 17289.029215 17289.023573
#> LO-REDUCTION     982 17289.027198 17289.023573
#> EXTENSION        984 17289.026203 17289.022041
#> LO-REDUCTION     986 17289.025551 17289.022041
#> LO-REDUCTION     988 17289.025295 17289.022041
#> LO-REDUCTION     990 17289.025230 17289.022041
#> HI-REDUCTION     992 17289.024774 17289.022041
#> REFLECTION       994 17289.023776 17289.021364
#> LO-REDUCTION     996 17289.023573 17289.021364
#> LO-REDUCTION     998 17289.023375 17289.021364
#> EXTENSION       1000 17289.023258 17289.020647
#> Exiting from Nelder Mead minimizer
#>     1002 function evaluations used
# extract results

# regression coefficient vector
res2$beta
#> [1] -2.32798372  0.58760533 -0.09963433  0.63521793  0.41916222 -0.36081430
#> [7]  0.27574700


# check likelihood
-ll(res2$beta)
#> [1] -17289.02

# number of calls to log lik function
# since optim does not return the number of
# iterations
res2$iter
#> function gradient 
#>     1002       NA 

# optim does not calculated working weights
head(res2$weights)
#> [1] 1
# }
```

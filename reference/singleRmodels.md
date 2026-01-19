# Family functions in singleRcapture package

Package `singleRcapture` utilizes various family type functions that
specify variable parts of population size estimation, regression,
diagnostics and other necessary information that depends on the model.
These functions are used as `model` argument in `estimatePopsize`
function.

## Usage

``` r
chao(lambdaLink = "loghalf", ...)

Hurdleztgeom(
  lambdaLink = c("log", "neglog"),
  piLink = c("logit", "cloglog", "probit"),
  ...
)

Hurdleztnegbin(
  nSim = 1000,
  epsSim = 1e-08,
  eimStep = 6,
  lambdaLink = c("log", "neglog"),
  alphaLink = c("log", "neglog"),
  piLink = c("logit", "cloglog", "probit"),
  ...
)

Hurdleztpoisson(
  lambdaLink = c("log", "neglog"),
  piLink = c("logit", "cloglog", "probit"),
  ...
)

oiztgeom(
  lambdaLink = c("log", "neglog"),
  omegaLink = c("logit", "cloglog", "probit"),
  ...
)

oiztnegbin(
  nSim = 1000,
  epsSim = 1e-08,
  eimStep = 6,
  lambdaLink = c("log", "neglog"),
  alphaLink = c("log", "neglog"),
  omegaLink = c("logit", "cloglog", "probit"),
  ...
)

oiztpoisson(
  lambdaLink = c("log", "neglog"),
  omegaLink = c("logit", "cloglog", "probit"),
  ...
)

zelterman(lambdaLink = "loghalf", ...)

zotgeom(lambdaLink = c("log", "neglog"), ...)

zotnegbin(
  nSim = 1000,
  epsSim = 1e-08,
  eimStep = 6,
  lambdaLink = c("log", "neglog"),
  alphaLink = c("log", "neglog"),
  ...
)

zotpoisson(lambdaLink = c("log", "neglog"), ...)

ztHurdlegeom(
  lambdaLink = c("log", "neglog"),
  piLink = c("logit", "cloglog", "probit"),
  ...
)

ztHurdlenegbin(
  nSim = 1000,
  epsSim = 1e-08,
  eimStep = 6,
  lambdaLink = c("log", "neglog"),
  alphaLink = c("log", "neglog"),
  piLink = c("logit", "cloglog", "probit"),
  ...
)

ztHurdlepoisson(
  lambdaLink = c("log", "neglog"),
  piLink = c("logit", "cloglog", "probit"),
  ...
)

ztgeom(lambdaLink = c("log", "neglog"), ...)

ztnegbin(
  nSim = 1000,
  epsSim = 1e-08,
  eimStep = 6,
  lambdaLink = c("log", "neglog"),
  alphaLink = c("log", "neglog"),
  ...
)

ztoigeom(
  lambdaLink = c("log", "neglog"),
  omegaLink = c("logit", "cloglog", "probit"),
  ...
)

ztoinegbin(
  nSim = 1000,
  epsSim = 1e-08,
  eimStep = 6,
  lambdaLink = c("log", "neglog"),
  alphaLink = c("log", "neglog"),
  omegaLink = c("logit", "cloglog", "probit"),
  ...
)

ztoipoisson(
  lambdaLink = c("log", "neglog"),
  omegaLink = c("logit", "cloglog", "probit"),
  ...
)

ztpoisson(lambdaLink = c("log", "neglog"), ...)
```

## Arguments

- lambdaLink:

  a link for Poisson parameter, `"log"` by default except for
  zelterman's and chao's models where only
  \\\ln\left(\frac{x}{2}\right)\\ is possible.

- ...:

  Additional arguments, not used for now.

- piLink:

  a link for probability parameter, `"logit"` by default

- nSim, epsSim:

  if working weights cannot be computed analytically these arguments
  specify maximum number of simulations allowed and precision level for
  finding them numerically respectively.

- eimStep:

  a non negative integer describing how many values should be used at
  each step of approximation of information matrixes when no analytic
  solution is available (e.g. `"ztnegbin"`), default varies depending on
  a function. Higher value usually means faster convergence but may
  potentially cause issues with convergence.

- alphaLink:

  a link for dispersion parameter, `"log"` by default

- omegaLink:

  a link for inflation parameter, `"logit"` by default

## Value

A object of class `family` containing objects:

- `makeMinusLogLike` – A factory function for creating the following
  functions: \\\ell(\boldsymbol{\beta}),
  \frac{\partial\ell}{\partial\boldsymbol{\beta}},
  \frac{\partial^{2}\ell}{\partial\boldsymbol{\beta}^{T}\partial\boldsymbol{\beta}}
  \\ functions from the \\\boldsymbol{y}\\ vector and the
  \\\boldsymbol{X}\_{vlm}\\ matrix (or just \\\boldsymbol{X}\\ if
  applied to model with single linear predictor)which has the `deriv`
  argument with possible values in `c(0, 1, 2)` that determine which
  derivative to return; the default value is `0`, which represents the
  minus log-likelihood.

- `links` – A list with link functions.

- `mu.eta, variance` – Functions of linear predictors that return
  expected value and variance. The `type` argument with 2 possible
  values (`"trunc"` and `"nontrunc"`) that specifies whether to return
  \\\mathbb{E}(Y\|Y\>0), \text{var}(Y\|Y\>0)\\ or \\\mathbb{E}(Y),
  \text{var}(Y)\\ respectively; the `deriv` argument with values in
  `c(0, 1, 2)` is used for indicating the derivative with respect to the
  linear predictors, which is used for providing standard errors in the
  `predict` method.

- `family` – A string that specifies name of the model.

- `valideta, validmu` – For now it only returns `TRUE`. In the near
  future, it will be used to check whether applied linear predictors are
  valid (i.e. are transformed into some elements of the parameter space
  subjected to the inverse link function).

- `funcZ, Wfun` – Functions that create pseudo residuals and working
  weights used in IRLS algorithm.

- `devResids` – A function that given the linear predictors prior
  weights vector and response vector returns deviance residuals. Not all
  family functions have these functions implemented yet.

- `pointEst, popVar` – Functions that given prior weights linear
  predictors and in the latter case also estimate of
  \\\text{cov}(\hat{\boldsymbol{\beta}})\\ and \\\boldsymbol{X\_{vlm}}\\
  matrix return point estimate for population size and analytic
  estimation of its variance.There is a additional boolean parameter
  `contr` in the former function that if set to true returns
  contribution of each unit.

- `etaNames` – Names of linear predictors.

- `densityFunction` – A function that given linear predictors returns
  value of PMF at values `x`. Additional argument `type` specifies
  whether to return \\\mathbb{P}(Y\|Y\>0)\\ or \\\mathbb{P}(Y)\\.

- `simulate` – A function that generates values of a dependent vector
  given linear predictors.

- `getStart` – An expression for generating starting points.

## Details

Most of these functions are based on some "base" distribution with
support \\\mathbb{N}\_{0}=\mathbb{N}\cup\lbrace 0\rbrace\\ that describe
distribution of \\Y\\ before truncation. Currently they include:
\\\mathbb{P}(Y=y\|\lambda,\alpha)=\left\lbrace \begin{array}{cc}
\frac{\lambda^{y}e^{-\lambda}}{y!} & \text{Poisson distribution} \cr
\frac{\Gamma(y+\alpha^{-1})}{\Gamma(\alpha^{-1})y!}
\left(\frac{\alpha^{-1}}{\alpha^{-1}+\lambda}\right)^{\alpha^{-1}}
\left(\frac{\lambda}{\alpha^{-1}+\lambda}\right)^{y} & \text{negative
binomial distribution} \cr \frac{\lambda^{y}}{(1+\lambda)^{y+1}} &
\text{geometric distribution} \end{array} \right.\\ where \\\lambda\\ is
the Poisson parameter and \\\alpha\\ is the dispersion parameter.
Geometric distribution is a special case of negative binomial
distribution when \\\alpha=1\\ it is included because negative binomial
distribution is quite troublesome numerical regression in fitting. It is
important to know that PMF of negative binomial distribution approaches
the PMF of Poisson distribution when \\\alpha\rightarrow 0^{+}\\.

**Note** in literature on single source capture recapture models the
dispersion parameter which introduces greater variability in negative
binomial distribution compared to Poisson distribution is generally
interpreted as explaining the *unobserved* heterogeneity i.e. presence
of important unobserved independent variables. All these methods for
estimating population size are tied to Poisson processes hence we use
\\\lambda\\ as parameter symbol instead of \\\mu\\ to emphasize this
connection. Also will not be hard to see that **all** estimators derived
from modifying the "base" distribution are unbiased if assumptions made
by respective models are not violated.

The **zero-truncated** models corresponding to "base" distributions are
characterized by relation: \\\mathbb{P}(Y=y\|Y\>0)=\left\lbrace
\begin{array}{cc} \frac{\mathbb{P}(Y=y)}{1-\mathbb{P}(Y=0)} & \text{when
}y\neq 0 \cr 0 & \text{when }y=0 \end{array}\right.\\ which allows us to
estimate parameter values using only observed part of population. These
models lead to the following estimates, respectively: \\ \begin{aligned}
\hat{N} &= \sum\_{k=1}^{N\_{obs}}\frac{1}{1-\exp(-\lambda\_{k})} &
\text{ For Poisson distribution} \cr \hat{N} &=
\sum\_{k=1}^{N\_{obs}}\frac{1}{1-(1+\alpha\_{k}\lambda\_{k})^{-\alpha\_{k}^{-1}}}
& \text{ For negative binomial distribution} \cr \hat{N} &=
\sum\_{k=1}^{N\_{obs}}\frac{1+\lambda\_{k}}{\lambda\_{k}} & \text{ For
geometric distribution} \end{aligned} \\

One common way in which assumptions of zero-truncated models are
violated is presence of **one-inflation** the presence of which is
somewhat similar in single source capture-recapture models to
zero-inflation in usual count data analysis. There are two ways in which
one-inflation may be understood, they relate to whether
\\\mathbb{P}(Y=0)\\ is modified by inflation. The first approach is
inflate (\\\omega\\ parameter) zero-truncated distribution as: \\
\mathbb{P}\_{new}(Y=y\|Y\>0) = \left\lbrace\begin{array}{cc} \omega +
(1 - \omega)\mathbb{P}\_{old}(Y=1\|Y\>0)& \text{when: } y = 1 \cr (1 -
\omega) \mathbb{P}\_{old}(Y=y\|Y\>0) & \text{when: } y \neq 1
\end{array}\right.\\ which corresponds to: \\ \mathbb{P}\_{new}(Y=y) =
\left\lbrace\begin{array}{cc} \mathbb{P}\_{old}(Y=0) & \text{when: } y =
0 \cr \omega(1 - \mathbb{P}(Y=0)) + (1 - \omega)\mathbb{P}\_{old}(Y=1) &
\text{when: } y = 1 \cr (1 - \omega) \mathbb{P}\_{old}(Y=y) &
\text{when: } y \> 1 \end{array}\right. \\ before zero-truncation.
Models that utilize this approach are commonly referred to as
*one-inflated zero-truncated models* (`oizt*()` family). Another way of
accommodating one-inflation in SSCR is by putting inflation parameter on
base distribution as: \\ \mathbb{P}\_{new}(Y=y) =
\left\lbrace\begin{array}{cc} \omega + (1 -
\omega)\mathbb{P}\_{old}(Y=1)& \text{when: } y = 1 \cr (1 - \omega)
\mathbb{P}\_{old}(Y=y) & \text{when: } y \neq 1 \end{array}\right. \\
which then becomes: \\ \mathbb{P}\_{new}(Y=y\|Y\>0) =
\left\lbrace\begin{array}{cc} \frac{\omega}{1 -
(1-\omega)\mathbb{P}\_{old}(Y=0)} + \frac{(1 - \omega)}{1 -
(1-\omega)\mathbb{P}\_{old}(Y=0)}\mathbb{P}\_{old}(Y=1)& \text{when: } y
= 1 \cr \frac{(1 - \omega)}{1 -
(1-\omega)\mathbb{P}\_{old}(Y=0)}\mathbb{P}\_{old}(Y=y) & \text{when: }
y \> 1 \end{array}\right. \\ after truncation and are commonly referred
to as *zero-truncated one-inflated models* (`ztoi*()` family). It was
shown by Böhning in 2022 paper that these approaches are equivalent in
terms of maximizing likelihoods if we do not put formula on \\\omega\\.
They can however lead to different population size estimates.

For *one-inflated zero-truncated models* the formula for population size
estimate \\\hat{N}\\ does not change since \\\mathbb{P}(y=0)\\ remains
the same but estimation of parameters changes all calculations.

For *zero-truncated one-inflated models* population size estimates are
expressed, respectively by: \\ \begin{aligned} \hat{N} &=
\sum\_{k=1}^{N\_{obs}}\frac{1}{1-(1-\omega\_{k})\exp(-\lambda\_{k})}
&\text{ For base Poisson distribution} \cr \hat{N} &=
\sum\_{k=1}^{N\_{obs}}\frac{1}{1-(1-\omega\_{k})(1+\alpha\_{k}\lambda\_{k})^{-\alpha\_{k}^{-1}}}
&\text{ For base negative binomial distribution} \cr \hat{N} &=
\sum\_{k=1}^{N\_{obs}}\frac{1+\lambda\_{k}}{\lambda\_{k} + \omega\_{k}}
&\text{ For base geometric distribution} \end{aligned} \\

**Zero-one-truncated** models ignore one counts instead of accommodating
one-inflation by utilizing the identity \\
\ell\_{\text{ztoi}}=\boldsymbol{f}\_{1}\ln{\frac{\boldsymbol{f}\_{1}}{N\_{obs}}}
+(N\_{obs}-\boldsymbol{f}\_{1})\ln{\left(1-\frac{\boldsymbol{f}\_{1}}{N\_{obs}}
\right)} + \ell\_{\text{zot}} \\ where \\\ell\_{\text{zot}}\\ is the log
likelihood of zero-one-truncated distribution characterized by
probability mass function: \\\mathbb{P}(Y=y\|Y\>1)=\left\lbrace
\begin{array}{cc}
\frac{\mathbb{P}(Y=y)}{1-\mathbb{P}(Y=0)-\mathbb{P}(Y=1)} & \text{when
}y \> 1 \cr 0 & \text{when }y\in\lbrace 0, 1\rbrace \end{array}\right.\\
where \\\mathbb{P}(Y)\\ is the probability mass function of the "base"
distribution. The identity above justifies use of zero-one-truncated,
unfortunately it was only proven for intercept only models, however
numerical simulations seem to indicate that even if the theorem cannot
be extended for (non trivial) regression population size estimation is
still possible.

For *zero-one-truncated models* population size estimates are expressed
by: \\ \begin{aligned} \hat{N} &= \boldsymbol{f}\_{1} +
\sum\_{k=1}^{N\_{obs}}
\frac{1-\lambda\_{k}\exp(-\lambda\_{k})}{1-\exp(-\lambda\_{k})-\lambda\_{k}\exp(-\lambda\_{k})}
&\text{ For base Poisson distribution} \cr \hat{N} &=
\boldsymbol{f}\_{1} + \sum\_{k=1}^{N\_{obs}}
\frac{1-\lambda\_{k}(1+\alpha\_{k}\lambda\_{k})^{-1-\alpha\_{k}^{-1}}}{
1-(1+\alpha\_{k}\lambda\_{k})^{-\alpha\_{k}^{-1}}-\lambda\_{k}(1+\alpha\_{k}\lambda\_{k})^{-1-\alpha\_{k}^{-1}}}
&\text{ For base negative binomial distribution} \cr \hat{N} &=
\boldsymbol{f}\_{1} + \sum\_{k=1}^{N\_{obs}}
\frac{\lambda\_{k}^{2}+\lambda\_{k}+1}{\lambda\_{k}^{2}} &\text{ For
base geometric distribution} \end{aligned} \\

Pseudo hurdle models are experimental and not yet described in
literature.

Lastly there are **chao** and **zelterman** models which are based on
logistic regression on the dummy variable \\ Z =
\left\lbrace\begin{array}{cc} 0 & \text{if }Y = 1 \cr 1 & \text{if }Y =
2 \end{array}\right.\\ based on the equation: \\ \text{logit}(p\_{k})=
\ln\left(\frac{\lambda\_{k}}{2}\right)=
\boldsymbol{\beta}\mathbf{x}\_{k}=\eta\_{k}\\ where \\\lambda\_{k}\\ is
the Poisson parameter.

The *zelterman* estimator of population size is expressed as:
\\\hat{N}=\sum\_{k=1}^{N\_{obs}}{1-\exp\left(-\lambda\_{k}\right)}\\ and
*chao* estimator has the form: \\
\hat{N}=N\_{obs}+\sum\_{k=1}^{\boldsymbol{f}\_{1}+\boldsymbol{f}\_{2}}
\frac{1}{\lambda\_{k}+ \frac{\lambda\_{k}^{2}}{2}} \\

## See also

[`estimatePopsize()`](https://ncn-foreigners.github.io/singleRcapture/reference/estimatePopsize.md)

## Author

Piotr Chlebicki, Maciej Beręsewicz

#' @import mathjaxr
NULL
#' \loadmathjax
#' @title Family functions in singleRcapture package
#' @author Piotr Chlebicki, Maciej Beręsewicz
#'
#' @description Package \code{singleRcapture} utilizes various family type
#' functions that specify variable parts of population size estimation,
#' regression, diagnostics and other necessary information that depends
#' on the model. These functions are used as \code{model} argument in
#' \code{estimatePopsize} function.
#' 
#' @param nSim,epsSim if working weights cannot be computed analytically these 
#' arguments specify maximum number of simulations allowed and
#' precision level for finding them numerically respectively.
#' @param lambdaLink a link for Poisson parameter, \code{"log"} 
#' by default except for zelterman's and chao's models where only 
#' \mjseqn{\ln\left(\frac{x}{2}\right)} is possible.
#' @param alphaLink a link for dispersion parameter, \code{"log"} by default
#' @param omegaLink a link for inflation parameter, \code{"logit"} by default 
#' @param piLink a link for probability parameter,  \code{"logit"} by default
#' @param eimStep a non negative integer describing 
#' how many values should be used at each step of approximation
#' of information matrixes when no analytic solution is available 
#' (e.g. \code{"ztnegbin"}), default varies depending on a function.
#' Higher value usually means faster convergence but may potentially cause
#' issues with convergence.
#' @param ... Additional arguments, not used for now.
#' 
#' @details Most of these functions are based on some "base" distribution with
#' support \mjseqn{\mathbb{N}_{0}=\mathbb{N}\cup\lbrace 0\rbrace} that describe
#' distribution of \mjseqn{Y} before truncation. Currently they include:
#' \mjsdeqn{\mathbb{P}(Y=y|\lambda,\alpha)=\left\lbrace
#' \begin{array}{cc}
#' \frac{\lambda^{y}e^{-\lambda}}{y!}    & \text{Poisson distribution}  \cr
#' \frac{\Gamma(y+\alpha^{-1})}{\Gamma(\alpha^{-1})y!} 
#' \left(\frac{\alpha^{-1}}{\alpha^{-1}+\lambda}\right)^{\alpha^{-1}}
#' \left(\frac{\lambda}{\alpha^{-1}+\lambda}\right)^{y} & 
#' \text{negative binomial distribution} \cr
#' \frac{\lambda^{y}}{(1+\lambda)^{y+1}} & 
#' \text{geometric distribution}
#' \end{array}
#' \right.}
#' where \mjseqn{\lambda} is the Poisson parameter and 
#' \mjseqn{\alpha} is the dispersion parameter. Geometric distribution
#' is a special case of negative binomial distribution when 
#' \mjseqn{\alpha=1} it is included because negative binomial 
#' distribution is quite troublesome numerical regression in fitting.
#' It is important to know that PMF of negative binomial distribution 
#' approaches the PMF of Poisson distribution when 
#' \mjseqn{\alpha\rightarrow 0^{+}}.
#' 
#' 
#' **Note** in literature on single source capture recapture models 
#' the dispersion parameter which introduces greater variability 
#' in negative binomial distribution compared to Poisson distribution is 
#' generally interpreted as explaining the *unobserved* heterogeneity
#' i.e. presence of important unobserved independent variables.
#' All these methods for estimating population size are tied to Poisson 
#' processes hence we use \mjseqn{\lambda} as parameter symbol
#' instead of \mjseqn{\mu} to emphasize this connection.
#' Also will not be hard to see that **all** estimators derived from 
#' modifying the "base" distribution are unbiased if assumptions 
#' made by respective models are not violated.
#' 
#' 
#' The **zero truncated** models corresponding to "base" distributions are
#' characterized by relation:
#' \mjsdeqn{\mathbb{P}(Y=y|Y>0)=\left\lbrace
#' \begin{array}{cc}
#' \frac{\mathbb{P}(Y=y)}{1-\mathbb{P}(Y=0)} & \text{when }y\neq 0 \cr
#' 0 & \text{when }y=0
#' \end{array}\right.}
#' which allows us to estimate parameter values using only observed part of
#' population. These models lead to the following estimates, respectively:
#' \mjsdeqn{
#' \begin{aligned}
#' \hat{N} &= \sum_{k=1}^{N_{obs}}\frac{1}{1-\exp(-\lambda_{k})} &
#' \text{ For Poisson distribution} \cr
#' \hat{N} &= \sum_{k=1}^{N_{obs}}\frac{1}{1-(1+\\alpha_{k}\lambda_{k})^{-\alpha_{k}^{-1}}} &
#' \text{ For negative binomial distribution} \cr
#' \hat{N} &= \sum_{k=1}^{N_{obs}}\frac{1+\lambda_{k}}{\lambda_{k}} &
#' \text{ For geometric distribution}
#' \end{aligned}
#' }
#' 
#' One common way in which assumptions of zero truncated models are violated is
#' presence of **one inflation** the presence of which is somewhat similar in
#' single source capture-recapture models to zero inflation in usual count data
#' analysis. There are two ways in which one inflation may be understood,
#' they relate to whether \mjseqn{\mathbb{P}(Y=0)} is 
#' modified by inflation. The first approach is inflate 
#' (\mjseqn{\omega} parameter) zero truncated distribution as:
#' \mjsdeqn{
#' \mathbb{P}_{new}(Y=y|Y>0) = \left\lbrace\begin{array}{cc}
#' \omega + (1 - \omega)\mathbb{P}_{old}(Y=1|Y>0)& \text{when: } y = 1 \cr
#' (1 - \omega) \mathbb{P}_{old}(Y=y|Y>0) & \text{when: } y \neq 1
#' \end{array}\right.}
#' which corresponds to:
#' \mjsdeqn{
#' \mathbb{P}_{new}(Y=y) = \left\lbrace\begin{array}{cc}
#' \mathbb{P}_{old}(Y=0) & \text{when: } y = 0 \cr
#' \omega(1 - \mathbb{P}(Y=0)) + (1 - \omega)\mathbb{P}_{old}(Y=1) & \text{when: } y = 1 \cr
#' (1 - \omega) \mathbb{P}_{old}(Y=y) & \text{when: } y > 1
#' \end{array}\right.
#' }
#' before zero truncation. Models that utilize this
#' approach are commonly referred to as *zero truncated one inflated models*.
#' Another way of accommodating one inflation in SSCR is by putting inflation
#' parameter on base distribution as:
#' \mjsdeqn{
#' \mathbb{P}_{new}(Y=y) = \left\lbrace\begin{array}{cc}
#' \omega + (1 - \omega)\mathbb{P}_{old}(Y=1)& \text{when: } y = 1 \cr
#' (1 - \omega) \mathbb{P}_{old}(Y=y) & \text{when: } y \neq 1
#' \end{array}\right.
#' }
#' which then becomes:
#' \mjsdeqn{
#' \mathbb{P}_{new}(Y=y|Y>0) = \left\lbrace\begin{array}{cc}
#' \frac{\omega}{1 - (1-\omega)\mathbb{P}_{old}(Y=0)} + \frac{(1 - \omega)}{1 - (1-\omega)\mathbb{P}_{old}(Y=0)}\mathbb{P}_{old}(Y=1)& \text{when: } y = 1 \cr
#' \frac{(1 - \omega)}{1 - (1-\omega)\mathbb{P}_{old}(Y=0)}\mathbb{P}_{old}(Y=y) & \text{when: } y > 1
#' \end{array}\right.
#' }
#' after truncation.
#' It was shown by Böhning in 2022 paper that these approaches are equivalent 
#' in terms of maximizing likelihoods if we do not put formula on 
#' \mjseqn{\omega}. They can however lead to different 
#' population size estimates.
#' 
#' For *zero truncated one inflated models* the formula for population size
#' estimate \mjseqn{\hat{N}} does not change since 
#' \mjseqn{\mathbb{P}(y=0)} remains the same but estimation of parameters 
#' changes all calculations.
#' 
#' For *one inflated zero truncated models* population size estimates are 
#' expressed, respectively by:
#' \mjsdeqn{
#' \begin{aligned}
#' \hat{N} &= \sum_{k=1}^{N_{obs}}\frac{1}{1-(1-\omega_{k})\exp(-\lambda_{k})} &\text{ For base Poisson distribution} \cr
#' \hat{N} &= \sum_{k=1}^{N_{obs}}\frac{1}{1-(1-\omega_{k})(1+\\alpha_{k}\lambda_{k})^{-\alpha_{k}^{-1}}} &\text{ For base negative binomial distribution} \cr
#' \hat{N} &= \sum_{k=1}^{N_{obs}}\frac{1+\lambda_{k}}{\lambda_{k} + \omega_{k}} &\text{ For base geometric distribution}
#' \end{aligned}
#' }
#' 
#' **Zero one truncated** models ignore one counts instead of accommodating
#' one inflation by utilizing the identity
#' \mjsdeqn{
#' \ell_{\text{ztoi}}=\boldsymbol{f}_{1}\ln{\frac{\boldsymbol{f}_{1}}{N_{obs}}}
#' +(N_{obs}-\boldsymbol{f}_{1})\ln{\left(1-\frac{\boldsymbol{f}_{1}}{N_{obs}}
#' \right)} + \ell_{\text{zot}}
#' }
#' where \mjseqn{\ell_{\text{zot}}} is the log likelihood 
#' of zero one truncated distribution characterized by probability mass function:
#' \mjsdeqn{\mathbb{P}(Y=y|Y>1)=\left\lbrace
#' \begin{array}{cc}
#' \frac{\mathbb{P}(Y=y)}{1-\mathbb{P}(Y=0)-\mathbb{P}(Y=1)} & \text{when }y > 1 \cr
#' 0 & \text{when }y\in\lbrace 0, 1\rbrace
#' \end{array}\right.}
#' where \mjseqn{\mathbb{P}(Y)} is the probability mass function of 
#' the "base" distribution. The identity above justifies use of zero one truncated,
#' unfortunately it was only proven for intercept only models, however
#' numerical simulations seem to indicate that even if the theorem cannot be
#' extended for (non trivial) regression population size estimation is still
#' possible. 
#' 
#' For *zero one truncated models* population size estimates are expressed by:
#' \mjsdeqn{
#' \begin{aligned}
#' \hat{N} &= \boldsymbol{f}_{1} + \sum_{k=1}^{N_{obs}}
#' \frac{1-\lambda_{k}\exp(-\lambda_{k})}{1-\exp(-\lambda_{k})-\lambda_{k}\exp(-\lambda_{k})} 
#' &\text{ For base Poisson distribution} \cr
#' \hat{N} &= \boldsymbol{f}_{1} + \sum_{k=1}^{N_{obs}}
#' \frac{1-\lambda_{k}(1+\alpha_{k}\lambda_{k})^{-1-\alpha_{k}^{-1}}}{
#' 1-(1+\alpha_{k}\lambda_{k})^{-\alpha_{k}^{-1}}-\lambda_{k}(1+\alpha_{k}\lambda_{k})^{-1-\alpha_{k}^{-1}}} 
#' &\text{ For base negative binomial distribution} \cr
#' \hat{N} &= \boldsymbol{f}_{1} + \sum_{k=1}^{N_{obs}}
#' \frac{\lambda_{k}^{2}+\lambda_{k}+1}{\lambda_{k}^{2}} 
#' &\text{ For base geometric distribution}
#' \end{aligned}
#' }
#' 
#' Pseudo hurdle models are experimental and not yet described in literature.
#' 
#' Lastly there are **chao** and **zelterman** models which are based on 
#' logistic regression on the dummy variable
#' \mjsdeqn{
#' Z = \left\lbrace\begin{array}{cc}
#' 0     & \text{if }Y = 1  \cr
#' 1     & \text{if }Y = 2
#' \end{array}\right.}
#' based on the equation:
#' \mjsdeqn{
#' \text{logit}(p_{k})=
#' \ln\left(\frac{\lambda_{k}}{2}\right)=
#' \boldsymbol{\beta}\mathbf{x}_{k}=\eta_{k}}
#' where \mjseqn{\lambda_{k}} is the Poisson parameter.
#' 
#' The *zelterman* estimator of population size is expressed as:
#' \mjsdeqn{\hat{N}=\sum_{k=1}^{N_{obs}}{1-\exp\left(-\lambda_{k}\right)}}
#' and *chao* estimator has the form:
#' \mjsdeqn{
#' \hat{N}=N_{obs}+\sum_{k=1}^{\boldsymbol{f}_{1}+\boldsymbol{f}_{2}}
#' \frac{1}{\lambda_{k}+ \frac{\lambda_{k}^{2}}{2}}
#' }
#' 
#' @seealso [estimatePopsize()]
#'
#' @return A object of class \code{family} containing objects:
#' \itemize{
#' \item \code{makeMinusLogLike} -- A factory function for creating the following functions:
#' \mjseqn{\ell(\boldsymbol{\beta}), \frac{\partial\ell}{\partial\boldsymbol{\beta}},
#' \frac{\partial^{2}\ell}{\partial\boldsymbol{\beta}^{T}\partial\boldsymbol{\beta}}
#' } functions from the
#' \mjseqn{\boldsymbol{y}} vector and the
#' \mjseqn{\boldsymbol{X}_{vlm}} matrix
#' (or just \mjseqn{\boldsymbol{X}} if applied to model 
#' with single linear predictor)which has the \code{deriv} argument with possible 
#' values in \code{c(0, 1, 2)} that determine which derivative to return; 
#' the default value is \code{0}, which represents the minus log-likelihood.
#' \item \code{links} -- A list with link functions.
#' \item \code{mu.eta, variance} -- Functions of linear predictors that
#' return expected value and variance. The \code{type} argument with 2 possible 
#' values (\code{"trunc"} and \code{"nontrunc"}) that specifies whether
#' to return \mjseqn{\mathbb{E}(Y|Y>0), \text{var}(Y|Y>0)} or 
#' \mjseqn{\mathbb{E}(Y), \text{var}(Y)} respectively; the \code{deriv} argument 
#' with values in \code{c(0, 1, 2)} is used for indicating the derivative with 
#' respect to the linear predictors, which is used for providing 
#' standard errors in the \code{predict} method.
#' \item \code{family} -- A string that specifies name of the model.
#' \item \code{valideta, validmu} -- For now it only returns \code{TRUE}. 
#' In the near future, it will be used to check whether applied linear 
#' predictors are valid (i.e. are transformed into some elements of the 
#' parameter space subjected to the inverse link function).
#' \item \code{funcZ, Wfun} -- Functions that create pseudo residuals and
#' working weights used in IRLS algorithm.
#' \item \code{devResids} -- A function that given the linear predictors
#' prior weights vector and response vector returns deviance residuals.
#' Not all family functions have these functions implemented yet.
#' \item \code{pointEst, popVar} -- Functions that given prior weights
#' linear predictors and in the latter case also estimate of 
#' \mjseqn{\text{cov}(\hat{\boldsymbol{\beta}})} and \mjseqn{\boldsymbol{X_{vlm}}}
#' matrix return point estimate for population size and analytic estimation 
#' of its variance.There is a additional boolean parameter \code{contr} in the 
#' former function that if set to true returns contribution of each unit.
#' \item \code{etaNames} -- Names of linear predictors.
#' \item \code{densityFunction} -- A function that given linear predictors 
#' returns value of PMF at values \code{x}. Additional argument \code{type}
#' specifies whether to return \mjseqn{\mathbb{P}(Y|Y>0)} or
#' \mjseqn{\mathbb{P}(Y)}.
#' \item \code{simulate} -- A function that generates values of a dependent 
#' vector given linear predictors.
#' \item \code{getStart} -- An expression for generating starting points.
#' }
#' @name singleRmodels
NULL

#' \loadmathjax
#' @title Regression diagnostics in singleRcapture
#' @author Piotr Chlebicki, Maciej Beręsewicz
#' 
#'
#' @description List of some regression diagnostics implemented for
#' \code{singleRStaticCountData} class. Functions that either require no changes from 
#' \code{glm} class or are not relevant to context of \code{singleRcapture}
#' are omitted.
#' 
#' @param model,object an object of \code{singleRStaticCountData} class.
#' @param dfbeta if \code{dfbeta} was already obtained it is possible to pass 
#' them into function so that they need not be computed for the second time.
#' @param do.coef logical indicating if \code{dfbeta} computation for influence 
#' should be done. \code{FALSE} by default.
#' @param cores a number of processor cores to be used,
#' any number greater than 1 activates code designed with \code{doParallel}, 
#' \code{foreach} and \code{parallel} packages. Note that for now using parallel 
#' computing makes tracing impossible so \code{trace} parameter is ignored in this case.
#' @param type a type of residual to return.
#' @param trace a logical value specifying whether to tracking results when
#' \code{cores > 1} it will result in a progress bar being created.
#' @param maxitNew the maximal number of iterations for regressions with starting 
#' points \mjseqn{\hat{\boldsymbol{\beta}}} on data 
#' specified at call for \code{model} after the removal of k'th row. By default 1.
#' @param ... arguments passed to other methods. 
#' Notably \code{dfpopsize.singleRStaticCountData} calls 
#' \code{dfbeta.singleRStaticCountData} if no \code{dfbeta} argument was 
#' provided and \code{controlMethod} is called in \code{dfbeta} method.
#' 
#' @details 
#' 
#' \code{dfpopsize} and \code{dfbeta} are closely related. \code{dfbeta} 
#' fits a regression after removing a specific row from the data and returns the 
#' difference between regression coefficients estimated on full data set and
#' data set obtained after deletion of that row, and repeats procedure once
#' for every unit present in the data.\code{dfpopsize} does the same for 
#' population size estimation utilizing coefficients computed by \code{dfbeta}.
#' 
#' \code{cooks.distance} is implemented (for now) only for models with a single
#' linear predictor and works exactly like the method for \code{glm} class.
#' 
#' \code{sigma} computes the standard errors of predicted means. Returns a matrix
#' with two columns first for truncated mean and the other for the non-truncated mean.
#' 
#' \code{residuals.singleRStaticCountData} (can be abbreviated to \code{resid})
#'  works like \code{residuals.glm} with the exception that:
#' \itemize{
#' \item \code{"pearson"} -- returns non standardized residuals.
#' \item \code{"pearsonSTD"} -- is currently defined only for single predictors
#' models but will be extended to all models in a near future, but for families 
#' with more than one distribution parameter it will be a multivariate residual.
#' \item \code{"response"} -- returns both residuals computed with truncated
#' and non truncated fitted value.
#' \item \code{"working"} -- is possibly multivariate if more than one linear 
#' predictor is present.
#' \item \code{"deviance"} -- is not yet defined for all families in 
#' [singleRmodels()] e.g. negative binomial based methods.
#' \item \code{"all"} -- returns all available residual types.
#' }
#' 
#' 
#' \code{hatvalues.singleRStaticCountData} is method for \code{singleRStaticCountData} 
#' class for extracting diagonal elements of projection matrix. 
#' 
#' Since \code{singleRcapture} supports not only regular glm's but also vglm's the 
#' \code{hatvalues} returns a matrix with number of columns corresponding to number 
#' of linear predictors in a model, where kth column corresponds to elements of 
#' the diagonal of projection matrix associated with kth linear predictor. 
#' For glm's  
#' \mjsdeqn{\boldsymbol{W}^{\frac{1}{2}}\boldsymbol{X}
#' \left(\boldsymbol{X}^{T}\boldsymbol{W}\boldsymbol{X}\right)^{-1}
#' \boldsymbol{X}^{T}\boldsymbol{W}^{\frac{1}{2}}}
#' where: \mjseqn{\boldsymbol{W}=\mathbb{E}\left(\text{Diag}
#' \left(\frac{\partial^{2}\ell}{\partial\boldsymbol{\eta}^{T}
#' \partial\boldsymbol{\eta}}\right)\right)}
#' and \mjseqn{\boldsymbol{X}} is a model (lm) matrix. 
#' For vglm's present in the package it is instead :
#' \mjsdeqn{\boldsymbol{X}_{vlm}
#' \left(\boldsymbol{X}_{vlm}^{T}\boldsymbol{W}\boldsymbol{X}_{vlm}\right)^{-1}
#' \boldsymbol{X}_{vlm}^{T}\boldsymbol{W}}
#' where:
#' \mjsdeqn{
#' \boldsymbol{W} = \mathbb{E}\left(\begin{bmatrix}
#' \text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta_{1}^{T}\partial\eta_{1}}\right) &
#' \text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta_{1}^{T}\partial\eta_{2}}\right) &
#' \dotso & \text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta_{1}^{T}\partial\eta_{p}}\right)\cr
#' \text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta_{2}^{T}\partial\eta_{1}}\right) &
#' \text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta_{2}^{T}\partial\eta_{2}}\right) &
#' \dotso & \text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta_{2}^{T}\partial\eta_{p}}\right)\cr
#' \vdots & \vdots & \ddots & \vdots\cr
#' \text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta_{p}^{T}\partial\eta_{1}}\right) &
#' \text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta_{p}^{T}\partial\eta_{2}}\right) &
#' \dotso & \text{Diag}\left(\frac{\partial^{2}\ell}{\partial\eta_{p}^{T}\partial\eta_{p}}\right)
#' \end{bmatrix}\right)}
#' is a block matrix constructed by taking the expected  value from diagonal 
#' matrixes corresponding to second derivatives with respect to each linear 
#' predictor (and mixed derivatives) and 
#' \mjseqn{\boldsymbol{X}_{vlm}} is a model (vlm) 
#' matrix constructed using specifications in \code{controlModel} and 
#' call to \code{estimatePopsize}.
#'
#' \code{influence} works like \code{glm} counterpart computing the most important
#' influence measures.
#'
#' @seealso [estimatePopsize()] [stats::hatvalues()] [controlMethod()] [stats::dfbeta()]
#' [stats::cooks.distance()]
#' @return \itemize{
#' \item For \code{hatvalues} -- A matrix with n rows and p columns where n is a 
#' number of observations in the data and p is number of regression parameters.
#' \item For \code{dfpopsize} -- A vector for which k'th element corresponds
#' to the difference between point estimate of population size estimation on 
#' full data set and point estimate of population size estimation after the 
#' removal of k'th unit from the data set.
#' \item For \code{dfbeta} -- A matrix with n rows and p observations where p 
#' is a number of units in data and p is the number of regression parameters. 
#' K'th row of this matrix corresponds to 
#' \mjseqn{\hat{\boldsymbol{\beta}}-\hat{\boldsymbol{\beta}}_{-k}}
#' where \mjseqn{\hat{\boldsymbol{\beta}}_{-k}} is a vector of estimates for 
#' regression parameters after the removal of k'th row from the data.
#' \item \code{cooks.distance} -- A matrix with a single columns with
#' values of cooks distance for every unit in \code{model.matrix}
#' \item \code{residuals.singleRStaticCountData} -- A \code{data.frame} 
#' with chosen residuals.
#' }
#' 
#' @examples
#' \donttest{
#' # For singleRStaticCountData class
#' # Get simple model
#' Model <- estimatePopsize(
#'   formula = capture ~ nation + age + gender, 
#'   data = netherlandsimmigrant, 
#'   model = ztpoisson, 
#'   method = "IRLS"
#' )
#' # Get dfbeta
#' dfb <- dfbeta(Model)
#' # The dfpopsize results are obtained via (It is also possible to not provide 
#' # dfbeta then they will be computed manually):
#' res <- dfpopsize(Model, dfbeta = dfb)
#' summary(res)
#' plot(res)
#' # see vaious types of residuals:
#' head(resid(Model, "all"))
#' }
#' @name regDiagSingleR
NULL
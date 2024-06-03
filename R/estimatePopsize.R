#' @import mathjaxr
NULL
#' @title Single source capture-recapture models
#' @author Piotr Chlebicki, Maciej Beręsewicz
#'
#' \loadmathjax
#' @description \code{estimatePopsize} first fits appropriate (v)glm model and 
#' then estimates full (observed and unobserved) population size.
#' In this types of models it is assumed that the response vector 
#' (i.e. the dependent variable) corresponds to the number of times a given unit 
#' was observed in the source.
#' Population size is then usually estimated by Horvitz-Thompson type estimator:
#' 
#' \mjsdeqn{\hat{N} = \sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)} = 
#' \sum_{k=1}^{N_{obs}}\frac{1}{1-\mathbb{P}(Y_{k}=0)}}
#'
#' where \mjseqn{I_{k}=I_{Y_{k} > 0}} are indicator 
#' variables, with value 1 if kth unit was observed at least once and 0 otherwise.
#'
#' @param data data frame or object coercible to data.frame class containing 
#' data for the regression and population size estimation.
#' @param formula formula for the model to be fitted, only applied to the "main" 
#' linear predictor. Only single response models are available.
#' @param ratioReg Not yet implemented
#' @param model model for regression and population estimate full description in [singleRmodels()]. 
#' @param weights optional object of prior weights used in fitting the model. 
#' Can be used to specify number of occurrences of rows in data see [controlModel()]
#' @param subset a logical vector indicating which observations should be used 
#' in regression and population size estimation. It will be evaluated on \code{data} argument provided on call.
#' @param naAction Not yet implemented.
#' @param method method for fitting values currently supported: iteratively 
#' reweighted least squares (\code{IRLS}) and maximum likelihood (\code{optim}).
#' @param popVar a method of constructing confidence interval either analytic 
#' or bootstrap. Bootstrap confidence interval type may be specified in 
#' \code{controlPopVar.} There is also the third possible value of \code{noEst} 
#' which skips the population size estimate all together.
#' @param controlMethod a list indicating parameters to use in fitting the model 
#' may be constructed with \code{singleRcapture::controlMethod} function. 
#' More information included in [controlMethod()].
#' @param controlModel a list indicating additional formulas for regression 
#' (like formula for inflation parameter/dispersion parameter) may be 
#' constructed with \code{singleRcapture::controlModel} function. 
#' More information will eventually be included in [controlModel()].
#' @param controlPopVar a list indicating parameters to use in estimating variance 
#' of population size estimation may be constructed with 
#' \code{singleRcapture::controlPopVar} function. 
#' More information included in [controlPopVar()].
#' @param modelFrame,x,y logical value indicating whether to return model matrix, 
#' dependent vector and model matrix as a part of output.
#' @param contrasts not yet implemented.
#' @param offset a matrix of offset values with number of columns matching the
#' number of distribution parameters providing offset values to each of 
#' linear predictors.
#' @param ... additional optional arguments passed to other methods eg. 
#' \code{estimatePopsizeFit}.
#' 
#' @details The generalized linear model is characterized by equation
#' \mjsdeqn{\boldsymbol{\eta}=\boldsymbol{X}\boldsymbol{\beta}}
#' where \mjseqn{\boldsymbol{X}} is the (lm) model matrix. 
#' The vector generalized linear model is similarly characterized by equations
#' \mjsdeqn{\boldsymbol{\eta}_{k}=\boldsymbol{X}_{k}\boldsymbol{\beta}_{k}}
#' where \mjseqn{\boldsymbol{X}_{k}} is a (lm) model 
#' matrix constructed from appropriate formula 
#' (specified in \code{controlModel} parameter).
#' 
#' The \mjseqn{\boldsymbol{\eta}} is then a vector constructed as:
#' 
#' \mjsdeqn{\boldsymbol{\eta}=\begin{pmatrix}
#' \boldsymbol{\eta}_{1} \cr 
#' \boldsymbol{\eta}_{2} \cr
#' \dotso \cr
#' \boldsymbol{\eta}_{p}
#' \end{pmatrix}^{T}}
#' 
#' and in cases of models in our package the (vlm) model matrix 
#' is constructed as a block matrix:
#' 
#' \mjsdeqn{\boldsymbol{X}_{vlm}=
#' \begin{pmatrix}
#' \boldsymbol{X}_{1} & \boldsymbol{0} &\dotso &\boldsymbol{0}\cr
#' \boldsymbol{0} & \boldsymbol{X}_{2} &\dotso &\boldsymbol{0}\cr
#' \vdots & \vdots & \ddots & \vdots\cr
#' \boldsymbol{0} & \boldsymbol{0} &\dotso &\boldsymbol{X}_{p}
#' \end{pmatrix}}
#' 
#' this differs from convention in \code{VGAM} package (if we only consider our 
#' special cases of vglm models) but this is just a convention and does not 
#' affect the model, this convention is taken because it makes fitting with 
#' IRLS (explanation of algorithm in [estimatePopsizeFit()]) algorithm easier.
#' (If \code{constraints} matrixes in \code{vglm} match the ones we implicitly
#' use the \code{vglm} model matrix differs with respect to order of 
#' \code{kronecker} multiplication of \code{X} and \code{constraints}.)
#' In this package we use observed likelihood to fit regression models.
#' 
#' As mentioned above usually the population size estimation is done via:
#' \mjsdeqn{\hat{N} = \sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)} = 
#' \sum_{k=1}^{N_{obs}}\frac{1}{1-\mathbb{P}(Y_{k}=0)}}
#'
#' where \mjseqn{I_{k}=I_{Y_{k} > 0}} are indicator variables, 
#' with value 1 if kth unit was observed at least once and 0 otherwise.
#' The \mjseqn{\mathbb{P}(Y_{k}>0)} are estimated by maximum likelihood.
#' 
#' The following assumptions are usually present when using 
#' the method of estimation described above:
#' 1. The specified regression model is correct. This entails linear
#' relationship between independent variables and dependent ones and
#' dependent variable being generated by appropriate distribution.
#' 2. No unobserved heterogeneity. If this assumption is broken there
#' are some possible (admittedly imperfect) workarounds see details in 
#' [singleRmodels()].
#' 3. The population size is constant in relevant time frame.
#' 4. Depending on confidence interval construction (asymptotic) normality
#' of \mjseqn{\hat{N}} statistic is assumed.
#' 
#' There are two ways of estimating variance of estimate \mjseqn{\hat{N}},
#' the first being \code{"analytic"} usually done by application of 
#' law of total variance to \mjseqn{\hat{N}}:
#'
#' \mjsdeqn{\text{var}(\hat{N})=\mathbb{E}\left(\text{var}
#' \left(\hat{N}|I_{1},\dots,I_{n}\right)\right)+
#' \text{var}\left(\mathbb{E}(\hat{N}|I_{1},\dots,I_{n})\right)}
#' 
#' and then by \mjseqn{\delta} method to 
#' \mjseqn{\hat{N}|I_{1},\dots I_{N}}:
#' 
#' \mjsdeqn{\mathbb{E}\left(\text{var}
#' \left(\hat{N}|I_{1},\dots,I_{n}\right)\right)=
#' \left.\left(\frac{\partial(N|I_1,\dots,I_N)}{\partial\boldsymbol{\beta}}\right)^{T}
#' \text{cov}\left(\boldsymbol{\beta}\right)
#' \left(\frac{\partial(N|I_1,\dots,I_N)}{\partial\boldsymbol{\beta}}\right)
#' \right|_{\boldsymbol{\beta}=\hat{\boldsymbol{\beta}}}}
#' 
#' and the \mjseqn{\text{var}\left(\mathbb{E}(\hat{N}|I_{1},\dots,I_{n})\right)}
#' term may be derived analytically (if we assume independence of
#' observations) since \mjseqn{\hat{N}|I_{1},\dots,I_{n}}
#' is just a constant. 
#' 
#' In general this gives us:
#' \mjsdeqn{
#' \begin{aligned}
#' \text{var}\left(\mathbb{E}(\hat{N}|I_{1},\dots,I_{n})\right)&=
#' \text{var}\left(\sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)}\right)\cr
#' &=\sum_{k=1}^{N}\text{var}\left(\frac{I_{k}}{\mathbb{P}(Y_{k}>0)}\right)\cr
#' &=\sum_{k=1}^{N}\frac{1}{\mathbb{P}(Y_{k}>0)^{2}}\text{var}(I_{k})\cr
#' &=\sum_{k=1}^{N}\frac{1}{\mathbb{P}(Y_{k}>0)^{2}}\mathbb{P}(Y_{k}>0)(1-\mathbb{P}(Y_{k}>0))\cr
#' &=\sum_{k=1}^{N}\frac{1}{\mathbb{P}(Y_{k}>0)}(1-\mathbb{P}(Y_{k}>0))\cr
#' &\approx\sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)^{2}}(1-\mathbb{P}(Y_{k}>0))\cr
#' &=\sum_{k=1}^{N_{obs}}\frac{1-\mathbb{P}(Y_{k}>0)}{\mathbb{P}(Y_{k}>0)^{2}}
#' \end{aligned}
#' }
#' 
#' Where the approximation on 6th line appears because in 5th line we sum over 
#' all units, that includes unobserved units, since \mjseqn{I_{k}} are 
#' independent and \mjseqn{I_{k}\sim b(\mathbb{P}(Y_{k}>0))} the 6th line
#' is an unbiased estimator of the 5th line.
#' 
#' The other method for estimating variance is \code{"bootstrap"}, but since
#' \mjseqn{N_{obs}=\sum_{k=1}^{N}I_{k}} is also a random variable bootstrap 
#' will not be as simple as just drawing \mjseqn{N_{obs}} units from data 
#' with replacement and just computing \mjseqn{\hat{N}}.
#' 
#' Method described above is referred to in literature as \code{"nonparametric"}
#' bootstrap (see [controlPopVar()]), due to ignoring variability in observed
#' sample size it is likely to underestimate variance.
#' 
#' A more sophisticated bootstrap procedure may be described as follows:
#' 1. Compute the probability distribution as: 
#' \mjsdeqn{
#' \frac{\hat{\boldsymbol{f}}_{0}}{\hat{N}}, 
#' \frac{\boldsymbol{f}_{1}}{\hat{N}}, 
#' \dots,
#' \frac{\boldsymbol{f}_{\max{y}}}{\hat{N}}}
#' where \mjseqn{\boldsymbol{f}_{n}} denotes observed 
#' marginal frequency of units being observed exactly n times.
#' 2. Draw \mjseqn{\hat{N}} units from the distribution above 
#' (if \mjseqn{\hat{N}} is not an integer than draw 
#' \mjseqn{\lfloor\hat{N}\rfloor + b(\hat{N}-\lfloor\hat{N}\rfloor)}),
#' where \mjseqn{\lfloor\cdot\rfloor} is the floor function.
#' 3. Truncated units with \mjseqn{y=0}.
#' 4. If there are covariates draw them from original data with replacement from 
#' uniform distribution. For example if unit drawn to new data has 
#' \mjseqn{y=2} choose one of covariate vectors from original data that 
#' was associated with unit for which was observed 2 times.
#' 5. Regress \mjseqn{\boldsymbol{y}_{new}} on \mjseqn{\boldsymbol{X}_{vlm new}} 
#' and obtain \mjseqn{\hat{\boldsymbol{\beta}}_{new}}, with starting point 
#' \mjseqn{\hat{\boldsymbol{\beta}}} to make it slightly faster, use them 
#' to compute \mjseqn{\hat{N}_{new}}.
#' 6. Repeat 2-5 unit there are at least \code{B} statistics are obtained.
#' 7. Compute confidence intervals based on \code{alpha} and \code{confType} 
#' specified in [controlPopVar()].
#' 
#' To do step 1 in procedure above it is convenient to first draw binary vector of length
#' \mjseqn{\lfloor\hat{N}\rfloor + b(\hat{N}-\lfloor\hat{N}\rfloor)} with probability 
#' \mjseqn{1-\frac{\hat{\boldsymbol{f}}_{0}}{\hat{N}}}, sum elements in that vector to
#' determine the sample size and then draw sample of this size uniformly from the data.
#' 
#' This procedure is known in literature as \code{"semiparametric"} bootstrap
#' it is necessary to assume that the have a correct estimate 
#' \mjseqn{\hat{N}} in order to use this type of bootstrap.
#' 
#' Lastly there is \code{"paramteric"} bootstrap where we assume that the 
#' probabilistic model used to obtain \mjseqn{\hat{N}} is correct the 
#' bootstrap procedure may then be described as:
#' 
#' 1. Draw \mjseqn{\lfloor\hat{N}\rfloor + b(\hat{N}-\lfloor\hat{N}\rfloor)} 
#' covariate information vectors with replacement from data according to 
#' probability distribution that is proportional to: \mjseqn{N_{k}},
#' where \mjseqn{N_{k}} is the contribution of kth unit i.e. 
#' \mjseqn{\dfrac{1}{\mathbb{P}(Y_{k}>0)}}.
#' 2. Determine \mjseqn{\boldsymbol{\eta}} matrix using estimate 
#' \mjseqn{\hat{\boldsymbol{\beta}}}.
#' 3. Generate \mjseqn{\boldsymbol{y}} (dependent variable) 
#' vector using \mjseqn{\boldsymbol{\eta}} and 
#' probability mass function associated with chosen model.
#' 4. Truncated units with \mjseqn{y=0} and construct 
#' \mjseqn{\boldsymbol{y}_{new}} and \mjseqn{\boldsymbol{X}_{vlm new}}.
#' 5. Regress \mjseqn{\boldsymbol{y}_{new}} on \mjseqn{\boldsymbol{X}_{vlm new}} 
#' and obtain \mjseqn{\hat{\boldsymbol{\beta}}_{new}}
#' use them to compute \mjseqn{\hat{N}_{new}}.
#' 6. Repeat 1-5 unit there are at least \code{B} statistics are obtained.
#' 7. Compute confidence intervals based on \code{alpha} and \code{confType}
#' specified in [controlPopVar()]
#' 
#' It is also worth noting that in the \code{"analytic"} method \code{estimatePopsize}
#' only uses "standard" covariance matrix estimation. It is possible that improper
#' covariance matrix estimate is the only part of estimation that has its assumptions
#' violated. In such cases post-hoc procedures are implemented in this package
#' to address this issue.
#' 
#' Lastly confidence intervals for \mjseqn{\hat{N}} are computed 
#' (in analytic case) either by assuming that it follows a normal distribution 
#' or that variable \mjseqn{\ln(N-\hat{N})}
#' follows a normal distribution. 
#' 
#' These estimates may be found using either \code{summary.singleRStaticCountData}
#' method or \code{popSizeEst.singleRStaticCountData} function. They're labelled as 
#' \code{normal} and \code{logNormal} respectively.
#'
#' @references General single source capture recapture literature:
#' 
#' Zelterman, Daniel (1988). ‘Robust estimation in truncated discrete distributions
#' with application to capture-recapture experiments’. In: Journal of statistical
#' planning and inference 18.2, pp. 225–237.
#' 
#' Heijden, Peter GM van der et al. (2003). ‘Point and interval estimation of the
#' population size using the truncated Poisson regression model’. 
#' In: Statistical Modelling 3.4, pp. 305–322. doi: 10.1191/1471082X03st057oa.
#' 
#' Cruyff, Maarten J. L. F. and Peter G. M. van der Heijden (2008). ‘Point and 
#' Interval Estimation of the Population Size Using a Zero-Truncated Negative 
#' Binomial Regression Model’. In: Biometrical Journal 50.6, pp. 1035–1050. 
#' doi: 10.1002/bimj.200810455
#' 
#' Böhning, Dankmar and Peter G. M. van der Heijden (2009). ‘A covariate adjustment 
#' for zero-truncated approaches to estimating the size of hidden and 
#' elusive populations’. In: The Annals of Applied Statistics 3.2, pp. 595–610. 
#' doi: 10.1214/08-AOAS214.
#' 
#' Böhning, Dankmar, Alberto Vidal-Diez et al. (2013). ‘A Generalization of 
#' Chao’s Estimator for Covariate Information’. In: Biometrics 69.4, pp. 1033– 
#' 1042. doi: 10.1111/biom.12082
#' 
#' Böhning, Dankmar and Peter G. M. van der Heijden (2019). ‘The identity of the 
#' zero-truncated, one-inflated likelihood and the zero-one-truncated likelihood 
#' for general count densities with an application to drink-driving in Britain’. 
#' In: The Annals of Applied Statistics 13.2, pp. 1198–1211. 
#' doi: 10.1214/18-AOAS1232.
#' 
#' Navaratna WC, Del Rio Vilas VJ, Böhning D. Extending Zelterman's approach for 
#' robust estimation of population size to zero-truncated clustered Data. 
#' Biom J. 2008 Aug;50(4):584-96. doi: 10.1002/bimj.200710441.
#' 
#' Böhning D. On the equivalence of one-inflated zero-truncated and zero-truncated 
#' one-inflated count data likelihoods. Biom J. 2022 Aug 15. doi: 10.1002/bimj.202100343.
#' 
#' Böhning, D., Friedl, H. Population size estimation based upon zero-truncated, 
#' one-inflated and sparse count data. Stat Methods Appl 30, 1197–1217 (2021). 
#' doi: 10.1007/s10260-021-00556-8
#' 
#' Bootstrap:
#' 
#' 
#' Zwane, PGM EN and Van der Heijden, Implementing the parametric bootstrap in capture-recapture 
#' models with continuous covariates 2003 Statistics & probability letters 65.2 pp 121-125
#' 
#' Norris, James L and Pollock, Kenneth H Including model uncertainty in estimating variances 
#' in multiple capture studies 1996 in Environmental and Ecological Statistics 3.3 pp 235-244
#' 
#' Vector generalized linear models: 
#' 
#' 
#' Yee, T. W. (2015). Vector Generalized Linear and Additive Models: 
#' With an Implementation in R. New York, USA: Springer. ISBN 978-1-4939-2817-0.
#'
#' @returns Returns an object of class \code{c("singleRStaticCountData", "singleR", "glm", "lm")}
#' with type \code{list} containing:\cr
#' \itemize{
#'  \item{\code{y} -- Vector of dependent variable if specified at function call.}
#'  \item{\code{X} -- Model matrix if specified at function call.}
#'  \item{\code{formula} -- A list with formula provided on call and additional formulas specified in \code{controlModel}.}
#'  \item{\code{call} -- Call matching original input.}
#'  \item{\code{coefficients} -- A vector of fitted coefficients of regression.}
#'  \item{\code{control} -- A list of control parameters for \code{controlMethod} and \code{controlModel}, \code{controlPopVar} is included in populationSize.}
#'  \item{\code{model} -- Model which estimation of population size and regression was built, object of class family.}
#'  \item{\code{deviance} -- Deviance for the model.}
#'  \item{\code{priorWeights} -- Prior weight provided on call.}
#'  \item{\code{weights} -- If \code{IRLS} method of estimation was chosen weights returned by \code{IRLS}, otherwise same as \code{priorWeights}.}
#'  \item{\code{residuals} -- Vector of raw residuals.}
#'  \item{\code{logL} -- Logarithm likelihood obtained at final iteration.}
#'  \item{\code{iter} -- Numbers of iterations performed in fitting or if \code{stats::optim} was used number of call to loglikelihood function.}
#'  \item{\code{dfResiduals} -- Residual degrees of freedom.}
#'  \item{\code{dfNull} -- Null degrees of freedom.}
#'  \item{\code{fittValues} -- Data frame of fitted values for both mu (the expected value) and lambda (Poisson parameter).}
#'  \item{\code{populationSize} -- A list containing information of population size estimate.}
#'  \item{\code{modelFrame} -- Model frame if specified at call.}
#'  \item{\code{linearPredictors} -- Vector of fitted linear predictors.}
#'  \item{\code{sizeObserved} -- Number of observations in original model frame.}
#'  \item{\code{terms} -- terms attribute of model frame used.}
#'  \item{\code{contrasts} -- contrasts specified in function call.}
#'  \item{\code{naAction} -- naAction used.}
#'  \item{\code{which} -- list indicating which observations were used in regression/population size estimation.}
#'  \item{\code{fittingLog} -- log of fitting information for \code{"IRLS"} fitting if specified in \code{controlMethod}.}
#' }
#' 
#' @seealso 
#' [stats::glm()] -- For more information on generalized linear models.
#' 
#' [stats::optim()] -- For more information on \code{optim} function used in 
#' \code{optim} method of fitting regression.
#' 
#' [controlMethod()] -- For control parameters related to regression.
#' 
#' [controlPopVar()] -- For control parameters related to population size estimation.
#' 
#' [controlModel()] -- For control parameters related to model specification.
#' 
#' [estimatePopsizeFit()] -- For more information on fitting procedure in
#' \code{esitmate_popsize}.
#' 
#' [popSizeEst()] [redoPopEstimation()] -- For extracting population size 
#' estimation results are applying post-hoc procedures.
#' 
#' [summary.singleRStaticCountData()] -- For summarizing important information about the
#' model and population size estimation results.
#' 
#' [marginalFreq()] -- For information on marginal frequencies and comparison
#' between observed and fitted quantities.
#' 
#' [VGAM::vglm()] -- For more information on vector generalized linear models.
#' 
#' [singleRmodels()] -- For description of various models.

#' @examples
#' \donttest{
#' # Model from 2003 publication 
#' # Point and interval estimation of the
#' # population size using the truncated Poisson regression mode
#' # Heijden, Peter GM van der et al. (2003)
#' model <- estimatePopsize(
#'    formula = capture ~ gender + age + nation, 
#'    data = netherlandsimmigrant, 
#'    model = ztpoisson
#' )
#' summary(model)
#' # Graphical presentation of model fit
#' plot(model, "rootogram")
#' # Statistical test
#' # see documentation for summary.singleRmargin
#' summary(marginalFreq(model), df = 1, "group")
#' 
#' # We currently support 2 methods of numerical fitting
#' # (generalized) IRLS algorithm and via stats::optim
#' # the latter one is faster when fitting negative binomial models 
#' # (and only then) due to IRLS having to numerically compute
#' # (expected) information matrixes, optim is also less reliable when
#' # using alphaFormula argument in controlModel
#' modelNegBin <- estimatePopsize(
#'     formula = TOTAL_SUB ~ ., 
#'     data = farmsubmission, 
#'     model = ztnegbin, 
#'     method = "optim"
#' )
#' summary(modelNegBin)
#' summary(marginalFreq(modelNegBin))
#' 
#' # More advanced call that specifies additional formula and shows
#' # in depth information about fitting procedure
#' pseudoHurdleModel <- estimatePopsize(
#'     formula = capture ~ nation + age,
#'     data = netherlandsimmigrant,
#'     model = Hurdleztgeom,
#'     method = "IRLS",
#'     controlMethod = controlMethod(verbose = 5),
#'     controlModel = controlModel(piFormula = ~ gender)
#' )
#' summary(pseudoHurdleModel)
#' # Assessing model fit
#' plot(pseudoHurdleModel, "rootogram")
#' summary(marginalFreq(pseudoHurdleModel), "group", df = 1)
#' 
#' 
#' # A advanced input with additional information for fitting procedure and
#' # additional formula specification and different link for inflation parameter.
#' Model <- estimatePopsize(
#'  formula = TOTAL_SUB ~ ., 
#'  data = farmsubmission, 
#'  model = oiztgeom(omegaLink = "cloglog"), 
#'  method = "IRLS", 
#'  controlMethod = controlMethod(
#'    stepsize = .85, 
#'    momentumFactor = 1.2,
#'    epsilon = 1e-10, 
#'    silent = TRUE
#'  ),
#'  controlModel = controlModel(omegaFormula = ~ C_TYPE + log_size)
#' )
#' summary(marginalFreq(Model), df = 18 - length(Model$coefficients))
#' summary(Model)
#' }
#' @export
estimatePopsize <- function(formula, 
                            ...) {
  UseMethod("estimatePopsize")
}
#' @rdname estimatePopsize
#' @importFrom stats model.frame model.matrix model.response lm predict
#' @export
estimatePopsize.default <- function(formula,
                                    data,
                                    model = c(
                                      "ztpoisson", "ztnegbin", "ztgeom",
                                      "zotpoisson", "ztoipoisson", "oiztpoisson",
                                      "ztHurdlepoisson", "Hurdleztpoisson", "zotnegbin",
                                      "ztoinegbin", "oiztnegbin", "ztHurdlenegbin",
                                      "Hurdleztnegbin", "zotgeom", "ztoigeom",
                                      "oiztgeom", "ztHurdlegeom", "ztHurdlegeom",
                                      "zelterman", "chao"
                                    ),
                                    ratioReg = FALSE,
                                    weights  = NULL,
                                    subset   = NULL,
                                    naAction = NULL,
                                    method   = c("optim", 
                                                 "IRLS"),
                                    popVar   = c("analytic",
                                                 "bootstrap",
                                                 "noEst"),
                                    controlMethod = NULL,
                                    controlModel  = NULL,
                                    controlPopVar = NULL,
                                    modelFrame    = TRUE,
                                    x             = FALSE,
                                    y             = TRUE,
                                    contrasts     = NULL,
                                    offset,
                                    ...) {
  if (missing(method)) method <- "IRLS"
  if (missing(popVar)) popVar <- "analytic"
  
  subset <- parse(text = deparse(substitute(subset)))
  
  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }
  
  if (!is.logical(subset)) {subset <- eval(subset, data)}
  if (is.null(subset)) {subset <- TRUE}
  # subset is often in conflict with some common packages hence explicit call
  data <- base::subset(data, subset = subset)
  
  family <- model
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  
  returnElements <- list(y, x, modelFrame)
  
  if (!ratioReg) {
    # adding control parameters that may possibly be missing
    # since passing simple lists as control arguments is allowed
    m1 <- controlPopVar
    m1 <- m1[sapply(m1, is.null) == FALSE]
    m2 <- controlPopVar(
      fittingMethod = match.arg(method), 
      bootstrapFitcontrol = controlMethod(
        epsilon = 1e-3, 
        maxiter = 20, 
        optimMethod = "Nelder-Mead", 
        silent = TRUE
      )
    )
    m2 <- m2[names(m2) %in% names(m1) == FALSE]
    controlPopVar <- append(m1, m2)
    
    m1 <- controlMethod
    m2 <- controlMethod(
      optimMethod = if (grepl(x = family$family, pattern = "negbin") || 
                        grepl(x = family$family, pattern = "^ztoi")) 
        "Nelder-Mead" else "BFGS",
      maxiter = if (isTRUE(method == "IRLS")) 100
                else 1000
    )
    m2 <- m2[names(m2) %in% names(m1) == FALSE]
    controlMethod <- append(m1, m2)
    
    m1 <- controlModel
    m2 <- controlModel()
    m2 <- m2[names(m2) %in% names(m1) == FALSE]
    controlModel <- append(m1, m2)
    
    #this may be changed for formula to be a list later
    formulas <- list(formula)
    
    if ("alpha" %in% family$etaNames) {
      formulas <- append(x = formulas, controlModel$alphaFormula)
    }
    if ("omega" %in% family$etaNames) {
      formulas <- append(x = formulas, controlModel$omegaFormula)
    }
    if ("pi" %in% family$etaNames) {
      formulas <- append(x = formulas, controlModel$piFormula)
    }
  
    combinedFromula <- singleRinternalMergeFormulas(formulas)
    
    modelFrame <- stats::model.frame(combinedFromula, 
                                     data,  
                                     ...)
    variables  <- stats::model.matrix(combinedFromula, 
                                      modelFrame, 
                                      contrasts = contrasts, 
                                      ...)
    
    terms     <- attr(modelFrame, "terms")
    contrasts <- attr(variables, "contrasts")
    observed  <- model.response(modelFrame)
    
    # It is both necessary and sufficient to check X_lm matrix to know whether
    # X_vlm matrix is full rank
    ## TODO:: Either start using Matrix package everywhere or use QR
    ## in fitting
    
    if (NCOL(observed) > 1) 
      stop("Single source capture-recapture models support only single dependent variable.")
  
    if (!is.null(weights)) {
      priorWeights <- as.numeric(weights)
    } else {
      priorWeights <- rep(1, nrow(modelFrame))
    }
    weights <- 1
    
    if (controlModel$weightsAsCounts) {
      sizeObserved <- sum(priorWeights)
    } else {
      sizeObserved <- nrow(data)
    }
    
    if (!all(observed > 0)) {
      stop("Error in function estimatePopsize, data contains zero-counts.")
    }
    # wch was here and is now deleted
    Xvlm <- singleRinternalGetXvlmMatrix(
      # this preserves terms attribute
      X = modelFrame,
      formulas = formulas, 
      parNames = family$etaNames
    )
    
    if (missing(offset)) {
      offset <- matrix(0, nrow = NROW(modelFrame), ncol = length(family$etaNames))
    } else if (!is.matrix(offset)) {
      offset <- matrix(offset, nrow = NROW(modelFrame), ncol = length(family$etaNames))
    }
    
    colnames(offset) <- family$etaNames
    
    if (!is.null(controlMethod$coefStart) || !is.null(controlMethod$etaStart)) {
      if (method == "IRLS") {
        etaStart  <- controlMethod$etaStart
        coefStart <- controlMethod$coefStart
        if (is.null(etaStart)) {
          etaStart <- matrix(Xvlm %*% coefStart, ncol = length(family$etaNames))
        }
      } else if (method == "optim") {
        etaStart  <- controlMethod$etaStart
        coefStart <- controlMethod$coefStart
        if (is.null(coefStart)) {
          eval(family$getStart)
        }
      }
    }
    else {
      eval(family$getStart)
    }
    
    
    FITT <- estimatePopsizeFit(
      y            = observed,
      X            = Xvlm,
      family       = family,
      control      = controlMethod,
      method       = method,
      priorWeights = priorWeights,
      coefStart    = coefStart,
      etaStart     = etaStart,
      offset       = offset
    )
    
    coefficients        <- FITT$beta
    names(coefficients) <- colnames(Xvlm)
    iter                <- FITT$iter
    dfReduced           <- sizeObserved * length(family$etaNames) - length(coefficients)
    IRLSlog             <- FITT$logg
    
    
    logLike <- family$makeMinusLogLike(
      y      = observed, 
      X      = Xvlm,
      weight = priorWeights,
      offset = offset
    )
    
    grad <- family$makeMinusLogLike(
      y      = observed, 
      X      = Xvlm, 
      weight = priorWeights, 
      deriv  = 1,
      offset = offset
    )
    
    hessian <- family$makeMinusLogLike(
      y      = observed,
      X      = Xvlm,
      weight = priorWeights,
      deriv  = 2,
      offset = offset
    )
  
    eta           <- matrix(as.matrix(Xvlm) %*% coefficients, 
                            ncol = length(family$etaNames)) + offset
    
    colnames(eta) <- family$etaNames
    rownames(eta) <- rownames(variables)
    weights       <- FITT$weights
    
    if (family$family == "zelterman") {
      eta <- matrix(as.matrix(variables) %*% coefficients, ncol = 1) + offset
      colnames(eta) <- family$etaNames
      rownames(eta) <- rownames(variables)
    }
  
    fitt <- data.frame(
      family$mu.eta(eta = eta),
      family$mu.eta(eta = eta, type = "nontrunc")
    )
    colnames(fitt) <- c("truncated", "nontruncated")
    
    # (Real square) Matrix is negative define iff all eigen values have 
    # negative sign. This is very fast. 
    # We only use eigen values so only.values is set to true.
    eig <- eigen(hessian(coefficients), symmetric = TRUE, only.values = TRUE)
    if (!all(sign(eig$values) == -1)) {
      warningMessage <- paste0(
        "The (analytically computed) hessian of the score function ",
        "is not negative define.\nNOTE: Second derivative test failing does not 
        necessarily mean that the maximum of score function that was found 
        numericaly is invalid since R^k is not a bounded space.\n",
        "Additionally in one inflated and hurdle models second ",
        "derivative test often fails even on valid arguments."
      )
      if (!isTRUE(controlMethod$silent)) warning(warningMessage)
      #cat("The eigen values were: ", eig$values) 
      # Add some option that will give much more 
      # information everywhere including here.
      if (controlPopVar$covType == "observedInform") {
        if (!isTRUE(controlMethod$silent)) 
          warning(paste0(
            "Switching from observed information matrix to Fisher information",
            " matrix because hessian of log-likelihood is not negative define."
          ))
        controlPopVar$covType <- "Fisher"
      }
    }
    
    nullDeviance <- as.numeric(NULL)
    LOG          <- -logLike(coefficients)
    resRes       <- priorWeights * (observed - fitt)
    
    if (family$family %in% c("zelterman", "chao")) {resRes <- resRes - 1}
  
    deviance <- sum(family$devResids(y   = observed, 
                                     wt  = priorWeights,
                                     eta = eta + offset) ^ 2)
    
    if (popVar == "noEst") {
      Pop <- NULL #TODO:: make sure methods are prepared for this
    } else {
      POP <- singleRcaptureinternalpopulationEstimate(
        y = observed,
        formulas = formulas,
        X = variables,
        grad = grad,
        hessian = hessian,
        popVar = popVar,
        weights = priorWeights,
        eta = eta,
        family = family,
        beta = coefficients,
        control = controlPopVar,
        Xvlm = Xvlm,
        W = if (method == "IRLS")
              weights 
            else 
              family$Wfun(prior = priorWeights, eta = eta),
        sizeObserved = sizeObserved,
        modelFrame = modelFrame,
        cov = NULL,
        offset = offset,
        weightsFlag = controlModel$weightsAsCounts
      )
    }
    
    structure(
      list(
        y                = if(isTRUE(returnElements[[1]])) as.numeric(observed) else NULL, # drop names
        X                = if(isTRUE(returnElements[[2]])) variables else NULL,
        modelFrame       = if (isTRUE(returnElements[[3]])) modelFrame else NULL,
        formula          = formulas,
        call             = match.call(),
        coefficients     = coefficients,
        control          = list(controlModel  = controlModel,
                                controlMethod = controlMethod),
        nullDeviance     = nullDeviance,
        model            = family,
        deviance         = deviance,
        priorWeights     = priorWeights,
        weights          = weights,
        residuals        = resRes,
        logL             = LOG,
        iter             = iter,
        dfResidual       = dfReduced,
        dfNull           = length(observed) - 1,
        fittValues       = fitt,
        populationSize   = POP,
        linearPredictors = eta,
        offset           = offset,
        sizeObserved     = sizeObserved,
        terms            = terms,
        contrasts        = contrasts,
        naAction         = naAction,
        fittingLog       = if (is.null(IRLSlog)) "IRLS logs were not saved." else IRLSlog
      ),
      class = c("singleRStaticCountData", "singleR", "glm", "lm")
    )
  } else {
    stop("Ratio regression is not yet implemented")
    # ff <- formula
    # # TODO make this inherit from family
    # a <- function(x) {x+1}
    # if (length(ff) == 3) {ff[[3]] <- 1}
    # modelFrame <- stats::model.frame(ff, data, ...)
    # 
    # observed <- modelFrame |>
    #   model.response() |>
    #   as.vector() # dropping names, won't be needed
    # 
    # delete <- modelFrame |>
    #   attr("names")
    # 
    # ff <- log(r) ~ 1
    # ff[[3]] <- formula[[3]]
    # 
    # counts <- table(observed)
    # 
    # if (TRUE) {
    #   r <- sapply(1:(max(observed)-1), function(x) {
    #     family$ratioFunc(x) * counts[as.character(x+1)] / counts[as.character(x)]
    #   }) |> as.vector()
    #   
    #   if (!is.null(weights)) {
    #     priorWeights <- as.numeric(weights)
    #   } else {
    #     priorWeights <- rep(1, nrow(modelFrame))
    #   }
    #   
    #   if (TRUE) {
    #     weights <- (1/counts[1:(max(observed)-1)] + 1/counts[2:max(observed)]) ^ -1
    #     #weights <- priorWeights * weights
    #   } else {
    #     # weighting outgh to be optional
    #   }
    #   
    #   # maybe include data here
    #   linearModel <- lm(ff, data = data.frame(
    #                     r = r, x = 1:(max(observed) - 1)
    #                     ), weights = weights, ...)
    #   
    #   fitRatio <- predict(linearModel, data.frame(x = 0:(max(observed)-1)))
    #   fitRatio <- exp(fitRatio)
    #   
    #   # first est for N
    #   N <- sum(counts) / (1 - 1 / sum(c(1, cumprod(fitRatio))))
    #   # second est for N
    #   N <- c(N, sum(counts) + family$ratioFunc(0)*counts["1"] / fitRatio[1])
    #   names(N) <- c("ht", "reg")
    #   f0 <- N - sum(counts)
    #   
    #   # TODO:: idk if this applies to HT estimate
    #   variation <- model.matrix(ff, model.frame(ff, data.frame(x = 0, r = 0)))
    #   variation <- f0^2 * as.vector(variation %*% vcov(linearModel) %*% t(variation))
    #   ## TODO -- this is an approximation
    #   # we're assuming var(counts["1"]) ~~ counts["1"]
    #   # this can be made better
    #   variation <- variation + counts["1"] * (a(0) ^ 2) * 
    #     exp(-predict(linearModel, data.frame(x = 0))) ^ 2
    #   
    #   sd <- sqrt(variation)
    #   sc <- qnorm(p = 1 - .05 / 2)
    #   
    #   
    #   confidenceInterval <- 
    #     data.frame(lowerBound = pmax(N - sc * sd, sum(counts)), 
    #                upperBound = N + sc * sd)
    #   print(N)
    #   print(confidenceInterval)
    #   
    # } else {
    #   ### TODO
    #   # put observed likelihood method here
    # }
    # 
    # structure(
    #   list(
    #     y                = if(isTRUE(returnElements[[1]])) as.numeric(observed) else NULL, # drop names
    #     X                = if(isTRUE(returnElements[[2]])) variables else NULL,
    #     modelFrame       = if (isTRUE(returnElements[[3]])) modelFrame else NULL,
    #     formula          = formula,
    #     call             = match.call(),
    #     coefficients     = coefficients
    #   ),
    #   class = c("singleRRatioReg", "singleR", "glm", "lm")
    # )
  }
}

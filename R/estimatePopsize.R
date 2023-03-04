#' @import mathjaxr
NULL
#' @title Single source capture-recapture models
#' @author Piotr Chlebicki, Maciej Beręsewicz
#'
#' @description \code{estimatePopsize} first fits appropriate (v)glm model and 
#' then estimates full (observed and unobserved) population size.
#' In this types of models it is assumed that the response vector 
#' (i.e. the dependent variable) corresponds to the number of times a given unit 
#' was observed in the source.
#' Population size is then usually estimated by Horvitz-Thompson type estimator:
#' 
#' \loadmathjax
#' \mjdeqn{\hat{N} = \sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)} = \sum_{k=1}^{N_{obs}}\frac{1}{1-\mathbb{P}(Y_{k}=0)}}{N = Sum_k=1^N I_k/P(Y_k > 0) = Sum_k=1^N_obs 1/(1-P(Y_k = 0))}
#'
#' where \mjeqn{I_{k}=I_{Y_{k} > 0}}{I_k=I_(Y_k > 0)} are indicator variables, 
#' with value 1 if kth unit was observed at least once and 0 otherwise.
#'
#' @param data Data frame or object coercible to data.frame class containing 
#' data for the regression and population size estimation.
#' @param formula Formula for the model to be fitted, only applied to the "main" 
#' linear predictor. Only single response models are available.
#' @param model Model for regression and population estimate full description in [singleRmodels()]. 
#' @param weights Optional object of a priori weights used in fitting the model.
#' @param subset A logical vector indicating which observations should be used 
#' in regression and population size estimation. It will be evaluated on \code{data} argument provided on call.
#' @param naAction Not yet implemented.
#' @param method Method for fitting values currently supported: iteratively 
#' reweighted least squares (\code{IRLS}) and maximum likelihood (\code{optim}).
#' @param popVar A method of constructing confidence interval either analytic or bootstrap.
#' Bootstrap confidence interval type may be specified in \code{controlPopVar.} 
#' There is also the third possible value of \code{noEst} which skips the 
#' population size estimate all together.
#' @param controlMethod A list indicating parameters to use in fitting the model 
#' may be constructed with \code{singleRcapture::controlMethod} function. 
#' More information included in [controlMethod()].
#' @param controlModel A list indicating additional formulas for regression 
#' (like formula for inflation parameter/dispersion parameter) may be 
#' constructed with \code{singleRcapture::controlModel} function. 
#' More information will eventually be included in [controlModel()].
#' @param controlPopVar A list indicating parameters to use in estimating variance 
#' of population size estimation may be constructed with 
#' \code{singleRcapture::controlPopVar} function. 
#' More information included in [controlPopVar()].
#' @param modelFrame,x,y Logical value indicating whether to return model matrix, 
#' dependent vector and model matrix as a part of output.
#' @param contrasts Not yet implemented.
#' @param ... Additional optional arguments passed to the following functions:
#' \itemize{
#'   \item \code{stats::model.frame} -- for creating data frame with all information about model specified with "main" formula.
#'   \item \code{stats::model.matrix} -- for creating model matrix (the lm matrix).
#'   \item \code{estimatePopsize.fit} -- possibly for picking starting points from zero truncated poisson regression.
#' } 
#' 
#' @details The generalized linear model is characterised by equation
#' \mjdeqn{\boldsymbol{\eta}=\boldsymbol{X}\boldsymbol{\beta}}{eta=X*beta}
#' where \mjeqn{\boldsymbol{X}}{X} is the (lm) model matrix. The vector 
#' generalized linear model is similarly characterised by equations
#' \mjdeqn{\boldsymbol{\eta}_{k}=\boldsymbol{X}_{k}\boldsymbol{\beta}_{k}}{eta_k=X_k*beta_k}
#' where \mjeqn{\boldsymbol{X}_{k}}{X_k} is a (lm) model matrix constructed
#' from appropriate formula (specified in \code{controlModel} parameter).
#' The \mjeqn{\boldsymbol{\eta}}{eta} is then a vector constructed as:
#' \mjdeqn{\boldsymbol{\eta}=\begin{pmatrix}\boldsymbol{\eta}_{1}^{T} & \boldsymbol{\eta}_{2}^{T} & \dotso & \boldsymbol{\eta}_{p}^{T}\end{pmatrix}^{T}}{eta = (eta_1', eta_2', ..., eta_p')'}
#' and in cases of models in our package the (vlm) model matrix 
#' is constructed as a block matrix:
#' \mjdeqn{\boldsymbol{X}_{vlm}=
#' \begin{pmatrix}
#' \boldsymbol{X}_{1} & \boldsymbol{0} &\dotso &\boldsymbol{0}\cr
#' \boldsymbol{0} & \boldsymbol{X}_{2} &\dotso &\boldsymbol{0}\cr
#' \vdots & \vdots & \ddots & \vdots\cr
#' \boldsymbol{0} & \boldsymbol{0} &\dotso &\boldsymbol{X}_{p}
#' \end{pmatrix}}{X_vlm = matrix(
#' X_1, 0, ..., 0
#' 0, X_2, ..., 0
#' ...........
#' 0, 0, ..., X_p)}
#' this differs from convention in \code{VGAM} package (if we only consider our 
#' special cases of vglm models) but this is just a convention and does not 
#' affect the model, this convention is taken because it makes fitting with 
#' IRLS (explanation of algorithm in [estimatePopsize.fit()]) algorithm easier.
#' (If \code{constraints} matrixes in \code{vglm} match the ones we implicitly
#' use the \code{vglm} model matrix differs with respect to order of 
#' \code{kronecker} multiplication of \code{X} and \code{constraints}.)
#' In this package we use observed likelihood to fit regression models.
#' 
#' As mentioned aboce usually the population size estimation is done via:
#' \mjdeqn{\hat{N} = \sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)} = \sum_{k=1}^{N_{obs}}\frac{1}{1-\mathbb{P}(Y_{k}=0)}}{N = Sum_k=1^N I_k/P(Y_k > 0) = Sum_k=1^N_obs 1/(1-P(Y_k = 0))}
#'
#' where \mjeqn{I_{k}=I_{Y_{k} > 0}}{I_k=I_(Y_k > 0)} are indicator variables, 
#' with value 1 if k'th unit was observed at least once and 0 otherwise.
#' The \mjeqn{\mathbb{P}(Y_{k}>0)}{P(Y_k > 0)} are estimated by maximum likelihood.
#' 
#' The following assumptions are usually present when using 
#' the method of estimation described above:
#' 1. The specified regression model is correct. This entails linear
#' relationship between independent variables and dependent ones and
#' dependent variable being generated by appropriate distribution.
#' 2. No unobserved heterogeneity. If this assumption is broken there
#' are some possible (admittedly imperfect) workarounds see details in [singleRmodels()].
#' 3. The population size is constant in relevant time frame.
#' 4. Depending on confidence interval construction (asymptotic) normality
#' of \mjeqn{\hat{N}}{N} statistic is assumed.
#' 
#' There are two ways of estimating variance of estimate \mjeqn{\hat{N}}{N},
#' the first being \code{"analytic"} usually done by application of 
#' law of total variance to \mjeqn{\hat{N}}{N}:
#' \mjdeqn{\text{var}(\hat{N})=\mathbb{E}\left(\text{var}
#' \left(\hat{N}|I_{1},\dots,I_{n}\right)\right)+
#' \text{var}\left(\mathbb{E}(\hat{N}|I_{1},\dots,I_{n})\right)}{
#' var(N)=E(var(N|I_1,...,I_N))+var(E(N|I_1,...,I_N))}
#' and then by \mjeqn{\delta}{delta} method to 
#' \mjeqn{\hat{N}|I_{1},\dotso I_{N}}{N|I_1,...,I_N}:
#' \mjdeqn{\mathbb{E}\left(\text{var}
#' \left(\hat{N}|I_{1},\dots,I_{n}\right)\right)=
#' \left.\left(\frac{\partial(N|I_1,...,I_N)}{\partial\boldsymbol{\beta}}\right)^{T}
#' \text{cov}\left(\boldsymbol{\beta}\right)
#' \left(\frac{\partial(N|I_1,...,I_N)}{\partial\boldsymbol{\beta}}\right)
#' \right|_{\boldsymbol{\beta}=\hat{\boldsymbol{\beta}}}}{E(var(N|I_1,...,I_N))=
#' (d(N|I_1,...,I_N)/d.beta)'cov(beta)(d(N|I_1,...,I_N)/d.beta)}
#' and the 
#' \mjeqn{\text{var}\left(\mathbb{E}(\hat{N}|I_{1},\dots,I_{n})\right)}{var(E(N|I_1,...,I_N))}
#' term may be derived analytically (if we assume independence of
#' observations) since 
#' \mjeqn{\hat{N}|I_{1},\dots,I_{n}}{N|I_1,...,I_N} is just a constant. 
#' In general this gives us:
#' \mjdeqn{
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
#' }{(var(E(N|I_1,...,I_N)) = var(sum_k=1^N I_k / P(Y_k > 0))
#' = sum_k=1^N var(I_k / P(Y_k > 0))
#' = sum_k=1^N var(I_k) / P(Y_k > 0)^2
#' = sum_k=1^N P(Y_k > 0)(1-P(Y_k > 0)) / P(Y_k > 0)^2
#' = sum_k=1^N (1-P(Y_k > 0)) / P(Y_k > 0)
#' ~~ sum_k=1^N (1-P(Y_k > 0)) I_k / P(Y_k > 0)^2
#' = sum_k=1^N_obs (1-P(Y_k > 0)) / P(Y_k > 0)^2)}
#' Where the approximation on 6th line appears because in 5th line we sum over all
#' units, that includes unobserved units, since \mjeqn{I_{k}}{I_k} are independent
#' and \mjeqn{I_{k}\sim b(\mathbb{P}(Y_{k}>0))}{
#' I_k follows bernoulli distribution with p = \mathbb{P}(Y_{k}>0)} the 6th line
#' is an unbiased estimator of the 5th line.
#' 
#' The other method for estimating variance is \code{"bootstrap"}, but since
#' \mjeqn{N_{obs}=\sum_{k=1}^{N}I_{k}}{N_obs = sum_k=1^N I_k} is also a random
#' variable bootstrap will not be as simple as just drawing \mjeqn{N_{obs}}{N_obs}
#' units from data with replacement and just computing \mjeqn{\hat{N}}{N}.
#' 
#' Method described above is refered to in literature as \code{"nonparametric"}
#' bootstrap (see [controlPopVar()]), due to ignoring variability in observed
#' sample size it is likely to underestimate variance.
#' 
#' A more sophisticated bootstrap procedure may be described as follows:
#' 1. Compute the probability distribution as: 
#' \mjdeqn{\frac{\hat{\boldsymbol{f}}_{0}}{\hat{N}}, \frac{\boldsymbol{f}_{1}}{\hat{N}}, \dotso, \frac{\boldsymbol{f}_{\max{y}}}{\hat{N}}}{f_0 / N, f_1 / N, ..., f_max(y) / N}
#' where \mjeqn{\boldsymbol{f}_{n}}{f_n} denotes observed marginal frequency of
#' units being observed exactly n times, round the quantitites above to nearest 
#' integer if necessary.
#' 2. Draw \mjeqn{\hat{N}}{N} units from the distribution above 
#' (if \mjeqn{\hat{N}}{N} is not an integer than draw \mjeqn{\lfloor\hat{N}\rfloor + b(\hat{N}-\lfloor\hat{N}\rfloor)}{floor(N) + b(N-floor(N))}).
#' 3. Truncated units with \mjeqn{y=0}{y=0}.
#' 4. If there are covariates draw them from original data with replacement from 
#' uniform distribution. For example if unit drawn to new data has 
#' \mjeqn{y=2}{y=2} choose one of covariate vectors from original data that 
#' was associated with unit for which was observed 2 times.
#' 5. Regress \mjeqn{\boldsymbol{y}_{new}}{y_new} on \mjeqn{\boldsymbol{X}_{vlm new}}{X_vlmNew}
#' and obtain \mjeqn{\hat{\boldsymbol{\beta}}_{new}}{beta_new}, with starting 
#' point \mjeqn{\hat{\boldsymbol{\beta}}}{beta} to make it slightly faster, 
#' use them to compute \mjeqn{\hat{N}_{new}}{N_new}.
#' 6. Repeat 2-5 unit there are at least \code{B} statistics are obtained.
#' 7. Compute confidence intervals based on \code{alpha} and \code{confType} 
#' specified in [controlPopVar()].
#' 
#' This procedure is known in literature as \code{"semiparametric"} bootstrap
#' it is necessary to assume that the have a correct estimate \mjeqn{\hat{N}}{N}
#' in order to use this type of bootstrap.
#' 
#' Lastly there is \code{"paramteric"} bootstrap where we assume that the 
#' probabilistic model used to obtain \mjeqn{\hat{N}}{N} is correct the 
#' bootstrap procedure may then be described as:
#' 1. Draw \mjeqn{\hat{N}}{N} covariate information vectors with replacement from
#' data according to probability distribution 
#' \mjdeqn{\frac{\lfloor N_{k}\rfloor + M_{k}}{\lfloor\hat{N}\rfloor}}{([N_k] + M_k)/[N]}
#' where \mjeqn{M_{k}\sim b(N_{k}-\lfloor N_{k}\rfloor)}{M_k ~ b(N_k - [N_k])}, 
#' \mjeqn{N_{k}}{N_{k}} is the contribution of kth unit i.e. 
#' \mjeqn{\frac{I_{k}}{\mathbb{P}(Y_{k}>0)}}{I_l/P(Y_k>0)} and
#' \mjeqn{\lfloor \cdot\rfloor}{[]} is the floor function.
#' 2. Determine \mjeqn{\boldsymbol{\eta}}{eta} matrix using estimate 
#' \mjeqn{\hat{\boldsymbol{\beta}}}{beta}.
#' 3. Generate \mjeqn{\boldsymbol{y}}{y} (dependent variable) vector using
#' \mjeqn{\boldsymbol{\eta}}{eta} and probability mass function associated with
#' chosen model.
#' 4. Truncated units with \mjeqn{y=0}{y=0} and construct 
#' \mjeqn{\boldsymbol{y}_{new}}{y_new} and 
#' \mjeqn{\boldsymbol{X}_{vlm new}}{X_vlmNew}.
#' 5. Regress \mjeqn{\boldsymbol{y}_{new}}{y_new} on 
#' \mjeqn{\boldsymbol{X}_{vlm new}}{X_vlmNew}
#' and obtain \mjeqn{\hat{\boldsymbol{\beta}}_{new}}{beta_new} 
#' use them to compute \mjeqn{\hat{N}_{new}}{N_new}.
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
#' Lastly confidence intervals for \mjeqn{\hat{N}}{N} are computed (in analytic case)
#' either by assuming that it follows a normal distribution or that variable
#' \mjeqn{\ln(N-\hat{N})}{log_e(N_true - N_estimate)} follows a normal 
#' distribution. These estimates may be found using either \code{summary.singleR}
#' method or \code{popSizeEst.singleR} function. They're labelled as 
#' \code{normal} and \code{logNormal} respectively.
#'
#' @references General single source capture recapture literature:
#' 
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
#' Vector generalised linear models: 
#' 
#' 
#' Yee, T. W. (2015). Vector Generalized Linear and Additive Models: 
#' With an Implementation in R. New York, USA: Springer. ISBN 978-1-4939-2817-0.
#'
#' @returns Returns an object of class \code{c("singleR", "glm", "lm")} containing:\cr
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
#'  \item{\code{iter} -- Numbers of iterations performed in fitting or if \code{stats::optim} was used number of call to loglikelihhod function.}
#'  \item{\code{dfResiduals} -- Residual degrees of freedom.}
#'  \item{\code{dfNull} -- Null degrees of freedom.}
#'  \item{\code{fittValues} -- Data frame of fitted values for both mu (the expected value) and lambda (Poisson parameter).}
#'  \item{\code{populationSize} -- A list containing information of population size estimate.}
#'  \item{\code{modelFrame} -- Model frame if specified at call.}
#'  \item{\code{linearPredictors} -- Vector of fitted linear predictors.}
#'  \item{\code{trcount} -- Number of truncated observations.}
#'  \item{\code{sizeObserved} -- Number of observations in original model frame.}
#'  \item{\code{terms} -- terms attribute of model frame used.}
#'  \item{\code{contrasts} -- contrasts specified in function call.}
#'  \item{\code{naAction} -- naAction used.}
#'  \item{\code{which} -- list indicating which observations were used in regression/population size estimation.}
#'  \item{\code{fittingLog} -- log of fitting information for \code{"IRLS"} fitting if specified in \code{controlMethod}.}
#' }
#' 
#' @seealso 
#' [stats::glm()] -- For more information on generalised linear models.
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
#' [estimatePopsize.fit()] -- For more information on fitting procedure in
#' \code{esitmate_popsize}.
#' 
#' [popSizeEst()] [redoPopEstimation()] -- For extracting population size 
#' estimation results are applying post-hoc procedures.
#' 
#' [summary.singleR()] -- For summarising important information about the
#' model and population size estimation results.
#' 
#' [marginalFreq()] -- For information on marginal frequencies and comparison
#' between observed and fitted quantities.
#' 
#' [VGAM::vglm()] -- For more information on vector generalised linear models.
#' 
#' [singleRmodels()] -- For description of various models.
#' @examples 
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
#' modelSingleRcapture <- estimatePopsize(formula = TOTAL_SUB ~ ., 
#'                                        data = farmsubmission, 
#'                                        model = ztnegbin, 
#'                                        method = "IRLS")
#' # comparison with VGAM package, VGAM uses slightly different parametrisation
#' # so we use negloglink instead of loglink for size parameter
#' # i.e 1 / dispersion parameter
#' if (require(VGAM)) {
#'   modelVGAM <- vglm(
#'      formula = TOTAL_SUB ~ ., 
#'      family = posnegbinomial(lsize = negloglink()), 
#'      data = farmsubmission
#'   )
#'   summary(modelVGAM)
#' }
#' # Results are comparable
#' summary(modelSingleRcapture)
#' summary(marginalFreq(modelSingleRcapture))
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
#' # additional formula specification.
#' Model <- estimatePopsize(
#'    formula = TOTAL_SUB ~ ., data = farmsubmission, 
#'    model = oiztgeom, 
#'    method = "IRLS", 
#'    controlMethod = controlMethod(stepsize = .05, 
#'                                  momentumFactor = 1.1, 
#'                                  epsilon = 1e-12, 
#'                                  silent = TRUE),
#'    controlModel = controlModel(omegaFormula = ~ .)
#' )
#' summary(marginalFreq(Model), df = 18 - length(Model$coefficients) - 1)
#' summary(Model)
#' @importFrom stats model.frame model.matrix model.response
#' @export
estimatePopsize <- function(formula,
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
                             weights = NULL,
                             subset = NULL,
                             naAction = NULL,
                             method = c("optim", 
                                        "IRLS", 
                                        "maxLik"),
                             popVar = c("analytic",
                                         "bootstrap",
                                         "noEst"),
                             controlMethod = NULL,
                             controlModel = NULL,
                             controlPopVar = NULL,
                             modelFrame = TRUE,
                             x = FALSE,
                             y = TRUE,
                             contrasts = NULL,
                             ...) {
  if (missing(method)) method <- "optim"
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
  # adding control parameters that may possibly be missing
  # since passing simple lists as control arguments is allowed
  m1 <- controlPopVar
  m1 <- m1[sapply(m1, is.null) == FALSE]
  m2 <- controlPopVar(
    fittingMethod = match.arg(method), 
    bootstrapFitcontrol = controlMethod(
      epsilon = 1e-3, 
      maxiter = 20, 
      optimMethod = if (grepl(x = family$family, pattern = "negbin") || 
                        grepl(x = family$family, pattern = "^ztoi")  || 
                        grepl(x = family$family, pattern = "^oizt")) 
        "Nelder-Mead" 
      else 
        "L-BFGS-B", 
      silent = TRUE
    )
  )
  m2 <- m2[names(m2) %in% names(m1) == FALSE]
  controlPopVar <- append(m1, m2)
  
  m1 <- controlMethod
  m2 <- controlMethod(
    optimMethod = if (grepl(x = family$family, pattern = "negbin") || 
                      grepl(x = family$family, pattern = "^ztoi")) 
      "Nelder-Mead" else "L-BFGS-B"
  )
  m2 <- m2[names(m2) %in% names(m1) == FALSE]
  controlMethod <- append(m1, m2)
  
  m1 <- controlModel
  m2 <- controlModel()
  m2 <- m2[names(m2) %in% names(m1) == FALSE]
  controlModel <- append(m1, m2)
  
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
  
  modelFrame <- stats::model.frame(combinedFromula, data,  ...)
  variables <- stats::model.matrix(combinedFromula, modelFrame, contrasts = contrasts, ...)
  terms <- attr(modelFrame, "terms")
  contrasts <- attr(variables, "contrasts")
  observed <- model.response(modelFrame)
  
  if (NCOL(observed) > 1) 
    stop("Single source capture-recapture models support only single dependent variable.")
  
  sizeObserved <- nrow(data) + controlPopVar$trcount

  if (!is.null(weights)) {
    priorWeights <- as.numeric(weights)
  } else {
    priorWeights <- rep(1, nrow(modelFrame))
  }
  weights <- 1
  
  if(!all(observed > 0)) {
    stop("Error in function estimatePopsize, data contains zero-counts.")
  }

  wch <- singleRcaptureinternalDataCleanupSpecialCases(family = family, 
                                                       observed = observed, 
                                                       popVar = popVar)

  controlPopVar$trcount <- controlPopVar$trcount + wch$trr
  
  Xvlm <- singleRinternalGetXvlmMatrix(X = subset(
    modelFrame, 
    select = colnames(modelFrame)[-(attr(terms, "response"))], 
    subset = wch$reg
  ), nPar = family$parNum, formulas = formulas, parNames = family$etaNames)
  
  
  start <- controlMethod$start #TODO:: Re-add use ztpoisson as start
  if (isTRUE(controlMethod$useZtpoissonAsStart)) 
    stop("useZtpoissonAsStart option is temporarily removed.")
  
  if (is.null(start)) {
    eval(family$getStart)
  }
  
  names(start) <- colnames(Xvlm)
  
  FITT <- estimatePopsize.fit(
    y = observed[wch$reg],
    X = Xvlm,
    family = family,
    control = controlMethod,
    method = method,
    priorWeights = priorWeights[wch$reg],
    start = start
  )
  
  coefficients <- FITT$beta
  names(coefficients) <- names(start)
  iter <- FITT$iter
  dfReduced <- nrow(Xvlm) - length(coefficients)
  IRLSlog <- FITT$logg
  
  
  logLike <- family$makeMinusLogLike(y = observed[wch$reg], X = Xvlm,
  weight = priorWeights[wch$reg])
  
  grad <- family$makeMinusLogLike(y = observed[wch$reg], X = Xvlm, 
  weight = priorWeights[wch$reg], deriv = 1)
  
  hessian <- family$makeMinusLogLike(y = observed[wch$reg], X = Xvlm,
  weight = priorWeights[wch$reg], deriv = 2)

  eta <- matrix(as.matrix(Xvlm) %*% coefficients, ncol = family$parNum)
  colnames(eta) <- family$etaNames
  rownames(eta) <- rownames(variables[wch$reg])
  weights <- FITT$weights
  
  if (family$family == "zelterman") {
    eta <- matrix(as.matrix(variables) %*% coefficients, ncol = 1)
    colnames(eta) <- family$etaNames
    rownames(eta) <- rownames(variables)
  }

  fitt <- data.frame(family$mu.eta(eta = eta),
  family$mu.eta(eta = eta, type = "nontrunc"))
  colnames(fitt) <- c("mu", "link")
  
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
  LOG <- -logLike(coefficients)
  resRes <- priorWeights * (observed[wch$reg] - fitt)
  if (family$family %in% c("zelterman", "chao")) {resRes <- resRes - 1}

  deviance <- sum(family$devResids(
    y = observed[wch$reg], 
    wt = priorWeights[wch$reg],
    eta = if (family$family == "zelterman") 
            eta[wch$reg] 
          else 
            eta
  ) ^ 2)
  
  if (popVar == "noEst") {
    Pop <- NULL #TODO:: make sure methods are prepared for this
  } else {
    POP <- singleRcaptureinternalpopulationEstimate(
      y = observed[wch$est],
      formulas = formulas,
      X = variables[wch$est, ],
      grad = grad,
      hessian = hessian,
      popVar = popVar,
      weights = priorWeights[wch$est],
      eta = eta,
      family = family,
      beta = coefficients,
      control = controlPopVar,
      Xvlm = if (family$family %in% c("zelterman", "chao") && popVar == "bootstrap") 
               variables 
             else 
               Xvlm,
      W = if (method == "IRLS")
            weights 
          else 
            family$Wfun(prior = priorWeights, eta = eta),
      sizeObserved = sizeObserved,
      modelFrame = modelFrame,
      cov = NULL
    )
  }
  structure(
    list(
      y = if(isTRUE(returnElements[[1]])) as.numeric(observed) else NULL, # drop names
      X = if(isTRUE(returnElements[[2]])) variables else NULL,
      formula = formulas,
      call = match.call(),
      coefficients = coefficients,
      control = list(controlModel = controlModel,
                     controlMethod = controlMethod),
      nullDeviance = nullDeviance,
      model = family,
      deviance = deviance,
      priorWeights = priorWeights,
      weights = weights,
      residuals = resRes,
      logL = LOG,
      iter = iter,
      dfResidual = dfReduced,
      dfNull = length(observed) - 1,
      fittValues = fitt,
      populationSize = POP,
      modelFrame = if (isTRUE(returnElements[[3]])) modelFrame else NULL,
      linearPredictors = eta,
      trcount = controlPopVar$trcount,
      sizeObserved = sizeObserved,
      terms = terms,
      contrasts = contrasts,
      naAction = naAction,
      which = wch,
      fittingLog = if (is.null(IRLSlog)) "IRLS logs were not saved." else IRLSlog
    ),
    class = c("singleR", "glm", "lm")
  )
}

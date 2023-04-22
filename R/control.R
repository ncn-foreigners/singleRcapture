#' @title Control parameters for regression
#' @author Piotr Chlebicki, Maciej Beręsewicz
#' \loadmathjax
#' @description \code{controlMethod} constructs a list with all necessary control parameters
#' for regression fitting in \code{estimatePopsize.fit} and \code{estimatePopsize}.
#'
#' @param epsilon Tolerance for fitting algorithms by default \code{1e-8}.
#' @param maxiter Maximum number of iterations.
#' @param verbose Value indicating whether to trace steps of fitting algorithm for 
#' \code{IRLS} fitting method different values of verbose give the following information:
#' \itemize{
#'   \item 1 -- Returns information on the number of current iteration and current log-likelihood.
#'   \item 2 -- Returns information on vector of regression parameters at current iteration (and all of the above).
#'   \item 3 -- Returns information on reduction of log-likelihood at current iteration (and all of the above).
#'   \item 4 -- Returns information on value of log-likelihood function gradient at current iteration (and all of the above).
#'   \item 5 -- Returns information on convergence criterion and values that are taken into account when considering convergence (and all of the above).
#' }
#' if \code{optim} method was chosen verbose will be passed to [stats::optim()] as trace.
#' @param printEveryN Integer value indicating how often to print information
#' specified in \code{verbose}, by default set to \code{1}.
#' @param start initial parameters for regression coefficients
#' if \code{NULL} they will be derived from simple poisson regression.
#' @param silent Logical, indicating whether warnings in \code{IRLS} method should be suppressed.
#' @param optimPass Optional list of parameters passed to \code{stats::optim(..., control = optimPass)}
#' if FALSE then list of control parameters will be inferred from other parameters.
#' @param optimMethod method of [stats::optim()] used L-BFGS-B is the default 
#' except for negative binomial and one inflated models where \code{"Nelder-Mead"} is used.
#' @param stepsize Only for \code{IRLS}, scaling of updates to \code{beta} vector 
#' lower value means slower convergence but more accuracy by default 1. 
#' In general if fitting algorithm fails lowering this value tends to 
#' be most effective at correcting it.
#' @param checkDiagWeights Logical value indicating whether to check if diagonal 
#' elements of working weights matrixes in \code{IRLS} are sufficiently positive 
#' so that these matrixes are positive defined. By default \code{TRUE}.
#' @param weightsEpsilon Small number to ensure positivity of weights matrixes. 
#' Only matters if \code{checkDiagWeights} is set to \code{TRUE}. 
#' By default \mjeqn{= 1\cdot 10^{-8}}{= 1 * 10^-10}
#' @param momentumFactor Experimental parameter in \code{IRLS} only allowing for 
#' taking previous step into account at current step, i.e instead of 
#' updating regression parameters as:
#' \mjdeqn{\boldsymbol{\beta}_{(a)} = \boldsymbol{\beta}_{(a-1)} + \text{stepsize} 
#' \cdot \text{step}_{(a)}}{beta_a = beta_a-1 + stepsize * step_a}
#' the update will be made as:
#' \mjdeqn{\boldsymbol{\beta}_{(a)} = \boldsymbol{\beta}_{(a-1)} + \text{stepsize} 
#' \cdot (\text{step}_{(a)} + \text{momentum}\cdot\text{step}_{(a-1)})}{
#' beta_a = beta_a-1 + stepsize * (step_a + momentum * step_a-1)}
#' @param useZtpoissonAsStart boolean value indicating whether to chose starting 
#' parameters from \code{ztpoisson} regression this one is useful mostly for 
#' various one inflated models.
#' @param momentumActivation the value of log-likelihood reduction bellow 
#' which momentum will apply.
#' @param criterion Criterion used to determine convergence in \code{IRLS}, 
#' multiple values may be provided. By default \code{c("coef", "abstol")}.
#' @param saveIRLSlogs Logical value indicating if information specified in
#' \code{verbose} should be saved to output object, by default \code{FALSE}.
#'
#' @return List with selected parameters, it is also possible to call list directly.
#' @seealso [singleRcapture::estimatePopsize()] [singleRcapture::estimatePopsize.fit()] [singleRcapture::controlModel()] [singleRcapture::controlPopVar()]
#' @export
controlMethod <- function(epsilon = 1e-8,
                          maxiter = 1000,
                          verbose = 0,
                          printEveryN = 1,
                          start = NULL,
                          optimMethod = "L-BFGS-B",
                          silent = FALSE,
                          optimPass = FALSE,
                          stepsize = 1,
                          checkDiagWeights = TRUE,
                          weightsEpsilon = 1e-8,
                          momentumFactor = 0,
                          saveIRLSlogs = FALSE,
                          useZtpoissonAsStart = FALSE,
                          momentumActivation = 5,
                          criterion = c("coef", "abstol", "reltol")) {
  if (!missing(criterion) && all(c("abstol", "reltol") %in% criterion)) {
    stop("Choosing both absolute tolerance and relative tolerance for convergence criterion is not allowed.")
  }
  if (!missing(criterion) && all(!(c("coef", "abstol", "reltol") %in% criterion))) {
    stop("At least one convergence criterion has to be chosen")
  }
  list(
    epsilon = epsilon,
    maxiter = maxiter,
    verbose = verbose,
    printEveryN = printEveryN,
    start = start,
    optimMethod = optimMethod,
    silent = silent,
    optimPass = optimPass,
    stepsize = stepsize,
    checkDiagWeights = checkDiagWeights,
    weightsEpsilon = weightsEpsilon,
    momentumFactor = momentumFactor,
    momentumActivation = momentumActivation,
    saveIRLSlogs = saveIRLSlogs,
    useZtpoissonAsStart = useZtpoissonAsStart,
    criterion = if (missing(criterion)) c("coef", "abstol") else criterion
  )
}
#' @title Control parameters specific to some models
#' @author Piotr Chlebicki, Maciej Beręsewicz
#' 
#' @description \code{controlModel} constructs a list with all necessary control parameters
#' in \code{estimatePopsize} that are either specific to selected model or don't fit
#' anywhere else.
#' 
#' Specifying additional formulas should be done by using only right hand side of
#' the formula also for now all variables from additional formulas should also be
#' included in the "main" formula.
#' 
#' @param weightsAsCounts For now does nothing. The plan is to have this indicate whether
#' \code{prior.weights} are to be treated as counts for sub populations and adjust all
#' necessary methods and functionalities, like adjustments in bootstrap or
#' decreasing weights in \code{dfbeta} instead or deleting rows from data, 
#' to accommodate this form of data.
#' @param omegaFormula Formula for inflation parameter in one inflated zero 
#' truncated and zero truncated one inflated models.
#' @param alphaFormula Formula for dispersion parameter in negative binomial
#' based models.
#' @param piFormula Formula for probability parameter in pseudo hurdle zero 
#' truncated and zero truncated pseudo hurdle models.
#' @return A list with selected parameters, it is also possible to call list directly.
#' @seealso [singleRcapture::estimatePopsize()] [singleRcapture::controlMethod()] [singleRcapture::controlPopVar()]
#' @export
controlModel <- function(weightsAsCounts = FALSE,
                         omegaFormula = ~ 1,
                         alphaFormula = ~ 1,
                         piFormula = ~ 1) {
  # This is not fully completed yet.
  list(
    weightsAsCounts = weightsAsCounts,
    omegaFormula = omegaFormula,
    alphaFormula = alphaFormula,
    piFormula = piFormula
  )
}
#' @title  Control parameters for population size estimation
#' @author Piotr Chlebicki, Maciej Beręsewicz
#'
#' @description Creating control parameters for population size estimation and 
#' respective standard error and variance estimation.
#'
#' @param alpha Significance level, 0.05 used by default.
#' @param trcount Truncated count - a number to be added to point estimator 
#' and both sides of confidence intervals.
#' @param bootType Bootstrap type. Default is \code{"parametric"}, 
#' other possible values are: \code{"semiparametric"} and \code{"nonparametric"}.
#' @param B Number of bootstrap samples to be performed (default 500).
#' @param confType Type of confidence interval for bootstrap confidence interval, 
#' \code{"percentile"} by default. 
#' Other possibility: \code{"studentized"} and \code{"basic"}.
#' @param keepbootStat Boolean value indicating whether to keep a vector of statistics produced by bootstrap.
#' @param traceBootstrapSize Boolean value indicating whether to print size of 
#' bootstrapped sample after truncation for semi- and fully parametric bootstraps.
#' @param bootstrapVisualTrace Boolean value indicating whether to plot bootstrap 
#' statistics in real time.
#' @param fittingMethod Method used for fitting models from bootstrap samples.
#' @param bootstrapFitcontrol Control parameters for each regression works exactly 
#' like \code{controlMethod} but for fitting models from bootstrap samples.
#' @param sd Indicates how to compute standard deviation of population 
#' size estimator either as:
#' \loadmathjax
#' \mjdeqn{\hat{\sigma}=\sqrt{\hat{\text{var}}(\hat{N})}}{sd=sqrt(var(N))}
#' for sqrt or for normalMVUE as the unbiased minimal variance estimator for normal distribution:
#' \mjdeqn{\hat{\sigma}=\sqrt{\hat{\text{var}}(\hat{N})}\frac{\Gamma\left((N_{obs}-1)/2\right)}{\Gamma\left(N_{obs}/2\right)}\sqrt{\frac{N_{obs}}{2}}}{sd=sqrt(var(N))sqrt(N_obs/2)Gamma(N_obs-1/2)/Gamma(N_obs/2)}
#' where the ration involving gamma functions is computed by loggamma function.
#' @param covType type of covariance matrix for regression parameters by default 
#' observed information matrix, more options will be here in the future.
#'
#' @return A list with selected parameters, it is also possible to call list directly.
#' @seealso [singleRcapture::estimatePopsize()] [singleRcapture::controlModel()] [singleRcapture::controlMethod()]
#' @export
controlPopVar <- function(alpha = .05,
                          trcount = 0,
                          bootType = c("parametric",
                                       "semiparametric",
                                       "nonparametric"),
                          B = 500,
                          confType = c("percentilic",
                                       "normal",
                                       "basic"), # TODO: add all
                          keepbootStat = TRUE,
                          traceBootstrapSize = FALSE,
                          bootstrapVisualTrace = FALSE,
                          fittingMethod = NULL,
                          bootstrapFitcontrol = NULL,
                          sd = c("sqrtVar", "normalMVUE"),
                          covType = c("observedInform",
                                      "Fisher")) {
  list(
    alpha = alpha,
    trcount = trcount,
    bootType = if(missing(bootType)) "parametric" else bootType,
    B = B,
    confType = if(missing(confType)) "percentilic" else confType,
    keepbootStat = keepbootStat,
    traceBootstrapSize = traceBootstrapSize,
    bootstrapVisualTrace = bootstrapVisualTrace,
    fittingMethod = fittingMethod,
    bootstrapFitcontrol = bootstrapFitcontrol,
    sd = if(missing(sd)) "sqrtVar" else sd,
    covType = if (missing(covType)) "observedInform" else covType
  )
}
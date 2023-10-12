#' @title Control parameters for regression
#' @author Piotr Chlebicki, Maciej Beręsewicz
#' \loadmathjax
#' 
#' @description \code{controlMethod} constructs a list with all necessary 
#' control parameters for regression fitting in 
#' \code{estimatePopsizeFit} and \code{estimatePopsize}.
#'
#' @param epsilon tolerance for fitting algorithms by default \code{1e-8}.
#' @param maxiter maximum number of iterations.
#' @param verbose value indicating whether to trace steps of fitting algorithm for 
#' \code{IRLS} fitting method different values of verbose give the following information:
#' \itemize{
#'   \item 1 -- Returns information on the number of current 
#'   iteration and current log-likelihood.
#'   \item 2 -- Returns information on vector of regression parameters 
#'   at current iteration (and all of the above).
#'   \item 3 -- Returns information on reduction of log-likelihood 
#'   at current iteration (and all of the above).
#'   \item 4 -- Returns information on value of log-likelihood function gradient 
#'   at current iteration (and all of the above).
#'   \item 5 -- Returns information on convergence criterion and values that are 
#'   taken into account when considering convergence (and all of the above).
#' }
#' if \code{optim} method was chosen verbose will be passed to [stats::optim()] as trace.
#' @param printEveryN integer value indicating how often to print information
#' specified in \code{verbose}, by default set to \code{1}.
#' @param coefStart,etaStart initial parameters for regression coefficients
#' or linear predictors if \code{NULL}. For \code{IRLS} fitting only \code{etaStart}
#' is needed so if \code{coefStart} is provided it will be converted to \code{etaStart},
#' for \code{optim} fitting \code{coefStart} is necessary and argument \code{etaStart}
#' will be ignored.
#' @param silent logical, indicating whether warnings in \code{IRLS} method should be suppressed.
#' @param optimPass optional list of parameters passed to \code{stats::optim(..., control = optimPass)}
#' if FALSE then list of control parameters will be inferred from other parameters.
#' @param optimMethod method of [stats::optim()] used  \code{"Nelder-Mead"} is the default .
#' @param stepsize only for \code{IRLS}, scaling of updates to \code{beta} vector 
#' lower value means slower convergence but more accuracy by default 1. 
#' In general if fitting algorithm fails lowering this value tends to 
#' be most effective at correcting it.
#' @param checkDiagWeights logical value indicating whether to check if diagonal 
#' elements of working weights matrixes in \code{IRLS} are sufficiently positive 
#' so that these matrixes are positive defined. By default \code{TRUE}.
#' @param weightsEpsilon small number to ensure positive definedness of weights matrixes. 
#' Only matters if \code{checkDiagWeights} is set to \code{TRUE}. 
#' By default \code{1e-8}.
#' @param momentumFactor experimental parameter in \code{IRLS} only allowing for 
#' taking previous step into account at current step, i.e instead of 
#' updating regression parameters as:
#' \mjsdeqn{\boldsymbol{\beta}_{(a)} = 
#' \boldsymbol{\beta}_{(a-1)} + \text{stepsize} \cdot \text{step}_{(a)}}
#' the update will be made as:
#' \mjsdeqn{
#' \boldsymbol{\beta}_{(a)} = \boldsymbol{\beta}_{(a-1)} + \text{stepsize} 
#' \cdot (\text{step}_{(a)} + \text{momentum}\cdot\text{step}_{(a-1)})}
#' @param momentumActivation the value of log-likelihood reduction bellow 
#' which momentum will apply.
#' @param criterion criterion used to determine convergence in \code{IRLS}, 
#' multiple values may be provided. By default \code{c("coef", "abstol")}.
#' @param saveIRLSlogs logical value indicating if information specified in
#' \code{verbose} should be saved to output object, by default \code{FALSE}.
#'
#' @return List with selected parameters, it is also possible to call list directly.
#' @seealso [singleRcapture::estimatePopsize()] [singleRcapture::estimatePopsizeFit()] 
#' [singleRcapture::controlModel()] [singleRcapture::controlPopVar()]
#' @export
controlMethod <- function(epsilon             = 1e-8,
                          maxiter             = 1000,
                          verbose             = 0,
                          printEveryN         = 1L,
                          coefStart           = NULL,
                          etaStart            = NULL,
                          optimMethod         = "Nelder-Mead",
                          silent              = FALSE,
                          optimPass           = FALSE,
                          stepsize            = 1,
                          checkDiagWeights    = TRUE,
                          weightsEpsilon      = 1e-8,
                          momentumFactor      = 0,
                          saveIRLSlogs        = FALSE,
                          momentumActivation  = 5,
                          criterion           = c("coef", 
                                                  "abstol", 
                                                  "reltol")) {
  
  if (!missing(criterion) && all(c("abstol", "reltol") %in% criterion)) 
    stop("Choosing both absolute tolerance and relative tolerance for convergence criterion is not allowed.")
  
  if (!missing(criterion) && all(!(c("coef", "abstol", "reltol") %in% criterion))) 
    stop("At least one convergence criterion has to be chosen")
  
  #if (isFALSE(is.integer(printEveryN) & printEveryN > 0))
  #  stop("printEveryN argument is either negative or not an integer")
  
  if (!isTRUE(is.finite(epsilon)) || !isTRUE(0 < epsilon) || isTRUE(length(epsilon) > 1))
    stop("Argument epsilon has to be a positive numeric value (of length 1).")
  
  if (!isTRUE(is.finite(maxiter)) || !isTRUE(0 < maxiter) || isTRUE(length(maxiter) > 1))
    stop("Argument maxiter has to be a positive numeric value (of length 1).")
  
  if (!isTRUE(is.finite(verbose)) || isTRUE(length(verbose) > 1))
    stop("Argument verbose has to be a numeric value (of length 1).")
  
  if (!isTRUE(is.finite(printEveryN)) || !isFALSE(0 > printEveryN) || isTRUE(length(printEveryN) > 1)) {
    stop("Argument printEveryN has to be a nonnegative numeric value (of length 1).")
    # If it's not an int
    printEveryN <- as.integer(printEveryN)
  }
  
  if (!is.null(etaStart) && !isTRUE(is.numeric(etaStart)))
    stop("Argument etaStart has to be either a numeric vector or NULL.")
  
  if (!is.null(coefStart) && !isTRUE(is.numeric(coefStart)))
    stop("Argument coefStart has to be either a numeric vector or NULL.")
  
  if (!isTRUE(is.logical(silent)) || isTRUE(length(silent) > 1))
    stop("Argument silent should be logical value (of length 1).")
  
  if (!isTRUE(is.logical(checkDiagWeights)) || isTRUE(length(checkDiagWeights) > 1))
    stop("Argument checkDiagWeights should be logical value (of length 1).")
  
  if (!isTRUE(is.logical(saveIRLSlogs)) || isTRUE(length(saveIRLSlogs) > 1))
    stop("Argument saveIRLSlogs should be logical value (of length 1).")
  
  if (!isTRUE(is.finite(stepsize)) || !isTRUE(0 < stepsize) || isTRUE(length(stepsize) > 1))
    stop("Argument stepsize has to be a positive numeric value (of length 1).")
  
  if (!isTRUE(is.finite(momentumActivation)) || !isTRUE(0 < momentumActivation) || 
       isTRUE(length(momentumActivation) > 1))
    stop("Argument momentumActivation has to be a positive numeric value (of length 1).")
  
  if (!isTRUE(is.finite(momentumFactor)) || !isFALSE(0 > momentumFactor) || 
       isTRUE(length(momentumFactor) > 1))
    stop("Argument momentumFactor has to be a non negative numeric value (of length 1).")
  
  list(
    epsilon             = epsilon,
    maxiter             = maxiter,
    verbose             = verbose,
    printEveryN         = printEveryN,
    etaStart            = etaStart,
    coefStart           = coefStart,
    optimMethod         = optimMethod,
    silent              = silent,
    optimPass           = optimPass,
    stepsize            = stepsize,
    checkDiagWeights    = checkDiagWeights,
    weightsEpsilon      = weightsEpsilon,
    momentumFactor      = momentumFactor,
    momentumActivation  = momentumActivation,
    saveIRLSlogs        = saveIRLSlogs,
    criterion           = if (missing(criterion)) c("coef", "abstol") else criterion
  )
}
#' @title Control parameters specific to some models
#' @author Piotr Chlebicki, Maciej Beręsewicz
#' 
#' @description \code{controlModel} constructs a list with all necessary 
#' control parameters in \code{estimatePopsize} that are either specific to 
#' selected model or do not fit anywhere else.
#' 
#' Specifying additional formulas should be done by using only right hand side of
#' the formula also for now all variables from additional formulas should also be
#' included in the "main" formula.
#' 
#' @param weightsAsCounts boolean value indicating whether to treat \code{weights}
#' argument as number of occurrences for each row in the \code{data} and adjust
#' necessary methods and functionalities, like adjustments in bootstrap or
#' decreasing weights in \code{dfbeta} instead or deleting rows from data, 
#' to accommodate this form of model specification.
#' @param omegaFormula formula for inflation parameter in one inflated zero 
#' truncated and zero truncated one inflated models.
#' @param alphaFormula formula for dispersion parameter in negative binomial
#' based models.
#' @param piFormula formula for probability parameter in pseudo hurdle zero 
#' truncated and zero truncated pseudo hurdle models.
#' 
#' @return A list with selected parameters, it is also possible to call list directly.
#' @seealso [singleRcapture::estimatePopsize()] [singleRcapture::controlMethod()] [singleRcapture::controlPopVar()] [singleRcapture::singleRmodels()]
#' @export
controlModel <- function(weightsAsCounts = FALSE,
                         omegaFormula = ~ 1,
                         alphaFormula = ~ 1,
                         piFormula = ~ 1) {
  # This is not fully completed yet.
  list(
    weightsAsCounts = weightsAsCounts,
    omegaFormula    = omegaFormula,
    alphaFormula    = alphaFormula,
    piFormula       = piFormula
  )
}
#' @title  Control parameters for population size estimation
#' @author Piotr Chlebicki, Maciej Beręsewicz
#' \loadmathjax
#'
#' @description Creating control parameters for population size estimation and 
#' respective standard error and variance estimation.
#'
#' @param alpha significance level, 0.05 used by default.
#' @param cores For bootstrap only, number of processor cores to be used,
#' any number greater than 1 activates code designed with \code{doParallel}, 
#' \code{foreach} and \code{parallel} packages. Note that for now using parallel
#' computing makes tracing impossible so \code{traceBootstrapSize} and 
#' \code{bootstrapVisualTrace} parameters are ignored in this case.
#' @param bootType bootstrap type. Default is \code{"parametric"}, 
#' other possible values are: \code{"semiparametric"} and \code{"nonparametric"}.
#' @param B number of bootstrap samples to be performed (default 500).
#' @param confType type of confidence interval for bootstrap confidence interval, 
#' \code{"percentile"} by default. 
#' Other possibilities: \code{"studentized"} and \code{"basic"}.
#' @param keepbootStat boolean value indicating whether to keep a vector of 
#' statistics produced by bootstrap.
#' @param traceBootstrapSize boolean value indicating whether to print size of 
#' bootstrapped sample after truncation for semi- and fully parametric bootstraps.
#' @param bootstrapVisualTrace boolean value indicating whether to plot bootstrap 
#' statistics in real time.
#' @param fittingMethod method used for fitting models from bootstrap samples.
#' @param bootstrapFitcontrol control parameters for each regression works exactly 
#' like \code{controlMethod} but for fitting models from bootstrap samples.
#' @param sd indicates how to compute standard deviation of population 
#' size estimator either as:
#' \mjsdeqn{\hat{\sigma}=\sqrt{\hat{\text{var}}(\hat{N})}}
#' for \code{sqrt} (which is slightly biased if \mjseqn{\hat{N}}
#' has a normal distribution) or for \code{normalMVUE} as the unbiased 
#' minimal variance estimator for normal distribution:
#' \mjsdeqn{\hat{\sigma}=\sqrt{\hat{\text{var}}(\hat{N})}
#' \frac{\Gamma\left(\frac{N_{obs}-1}{2}\right)}{\Gamma\left(\frac{N_{obs}}{2}\right)}
#' \sqrt{\frac{N_{obs}}{2}}}
#' where the ration involving gamma functions is computed by log gamma function.
#' @param covType type of covariance matrix for regression parameters by default 
#' observed information matrix.
#'
#' @return A list with selected parameters, it is also possible to call list directly.
#' @seealso [singleRcapture::estimatePopsize()] [singleRcapture::controlModel()] [singleRcapture::controlMethod()]
#' @export
controlPopVar <- function(alpha = .05,
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
                          fittingMethod = c("optim", "IRLS"),
                          bootstrapFitcontrol = NULL,
                          sd = c("sqrtVar", "normalMVUE"),
                          covType = c("observedInform", "Fisher"),
                          cores = 1L) {
  
  if (missing(fittingMethod)) fittingMethod <- "IRLS"
  if (missing(bootType)) bootType <- "parametric"
  if (missing(confType)) confType <- "percentilic"
  if (missing(covType))  covType <- "observedInform"
  if (missing(sd))       sd <- "sqrtVar"
  
  # I'm using !isTRUE instead of isFALSE because I do not wan't nulls to pass the tests
  if (!isTRUE(is.finite(alpha)) || isTRUE(0 > alpha || 1 < alpha) || isTRUE(length(alpha) > 1))
    stop("Argument alpha should be a numeric value between 0 and 1 (of length 1).")
  
  if (!isTRUE(sd %in% c("sqrtVar", "normalMVUE")) && !missing(sd))
    stop("Argument sd should be a character with value either sqrtVar or normalMVUE.")
  
  if ((!isTRUE(covType %in% c("observedInform", "Fisher")) && !missing(covType)) || isTRUE(length(covType) > 1))
    stop("Argument covType should be a character with value either observedInform or Fisher.")
  
  
  if (!missing(bootType) && !isTRUE(bootType %in% c("parametric",
                                                    "semiparametric",
                                                    "nonparametric")))
    stop("Argument bootType should be a character with value either in c(parametric, semiparametric, nonparametric).")
  
  if (!isTRUE(is.finite(B)) || isTRUE(B < 0) || isTRUE(length(B) > 1))
    stop("Argument B should be a numeric value greater than 0 (of length 1).")
  
  if (!isTRUE(all(c(is.logical(traceBootstrapSize), is.logical(bootstrapVisualTrace)))) || 
       isTRUE(length(bootstrapVisualTrace) > 1) || isTRUE(length(traceBootstrapSize) > 1))
    stop("Arguments traceBootstrapSize and bootstrapVisualTrace should be logical values (of length 1).")
    
    
  list(
    bootstrapVisualTrace = bootstrapVisualTrace,
    bootstrapFitcontrol  = bootstrapFitcontrol,
    traceBootstrapSize   = traceBootstrapSize,
    fittingMethod        = fittingMethod,
    keepbootStat         = keepbootStat,
    bootType             = bootType,
    confType             = confType,
    covType              = covType,
    cores                = cores,
    alpha                = alpha,
    B                    = B,
    sd                   = sd
  )
}

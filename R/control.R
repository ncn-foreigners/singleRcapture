#' Control parameters for regression
#'
#' @param epsilon relative tolerance (absolute for robust method) for fitting algorithms by default 1e-8
#' @param dispEpsilon tolerance for estimation of dispersion parameter only matters for negbins in robust
#' @param maxiter Maximal number of iterations
#' @param verbose Value indicating whether to trace steps of fitting algorithm for robust its either 0 (for no tracing), 1 (for tracing logarithm likelihhod) or 2 (for tracing logarithm likelihood and vector of regression parameters) for mle it is passed to  stats::optim as trace
#' @param start initial parameters for regression if NULL they will be derived from simple poisson regression
#' @param dispersionstart initial parameters for dispersion parameter if applies
#' @param omegastart initial parameters for inflation parameter if applies
#' @param mleMethod method of stats::optim used L-BFGS-B is the default except for negative binomial models where Nelder-Mead is used
#' @param silent logical indicating whether warnings in robust method should be suppressed
#' @param optimPass optional parameter allowing for passing list into as stats::optim(..., control = optimPass) if FALSE then list of controlparameters will be infered from other parameters
#' @param dispGiven logical indicates whether to estimate dispersion parameter or assume the one provided is already correct
#' @param thetaGiven TODO
#' @param stepsize Only for robust, scaling of stepsize lower value means slower convergence but more accuracy by default 1
#'
#' @return A list with selected parameters, it is also possible to call list directly.
#' @export
control.method <- function(epsilon = 1e-8,
                           dispEpsilon = 1e-5,
                           maxiter = 1000,
                           verbose = 0,
                           start = NULL,
                           dispersionstart = NULL,
                           omegastart = NULL,
                           mleMethod = "L-BFGS-B",
                           silent = FALSE,
                           optimPass = FALSE,
                           dispGiven = FALSE,
                           thetaGiven = FALSE,
                           stepsize = 1) {
  list(
    epsilon = epsilon,
    dispEpsilon = dispEpsilon,
    maxiter = maxiter,
    verbose = verbose,
    start = start,
    dispersionstart = dispersionstart,
    omegastart = omegastart,
    mleMethod = mleMethod,
    silent = silent,
    optimPass = optimPass,
    dispGiven = dispGiven,
    stepsize = stepsize,
    thetaGiven = thetaGiven
  )
}
#' control.model
#' TODO
#' @param weightsAsCounts TODO
#' @param omegaIntercept TODO
#' @param alphaIntercept TODO
#' for now does nothing
#' @return control.model
#' @export
control.model <- function(weightsAsCounts = FALSE,
                          omegaIntercept = NULL,# ????
                          alphaIntercept = TRUE # ????
                          # I suspect no other parameters will be used in the whole package so maybe its better to just specify them in control
                          # instead of making formula argument in main function a list or something like that.
                          # In VGAM there is a zero argument when calling vgam family class functions it specifies which argument
                          # is to be modeled as intercept only eg. in ztnegbin zero is by default on alpha so their fitting is 
                          # the same as ours. The only other aptions are for lambda or nothing i.e. there is no manual formula specification
                          ) {
  # TODO
  list(
    weightsAsCounts = weightsAsCounts,
    omegaIntercept = omegaIntercept,
    alphaIntercept = alphaIntercept
  )
}
#' Control parameters for population size estimation
#'
#' @param alpha Significance level, a number from range (0, 1) 5% by default.
#' @param trcount Truncated count - a number to be added to point estimator and both sides of confidence intervals.
#' @param bootType Type of bootstrap performed, by default Parametric, other posible values are: Semiparametric and Nonparametric.
#' @param B Number of bootstrap samples to be performed.
#' @param confType Type of confidence interval for bootstrap confidence interval, percentilic by default.
#' @param keepbootStat Boolean value indicating whether to keep a vector of statistic produced by bootstrap, for large values of B it may significant amount of size.
#' @param traceBootstrapSize Boolean value indicating whether to print size of bootstrapped sample after truncation for semi and fully parametric boostraps.
#' @param fittingMethod method used for fitting models from boostrap samples either "robust" or "mle" is left as NULL will be chosen automatically.
#' @param bootstrapFitcontrol control parameters for each regression works exactly like control.method.
#' @param covType type of covariance matrix for regression parameters by default observed information matrix, more options will be here in the future.
#'
#' @return A list with selected parameters, it is also possible to call list directly.
#' @export
control.pop.var <- function(alpha = .05,
                            trcount = 0,
                            bootType = c("parametric",
                                         "semiparametric",
                                         "nonparametric"),
                            B = 500,
                            confType = c("percentilic",
                                         "studentized",
                                         "basic"), # TODO: add all
                            keepbootStat = TRUE,
                            traceBootstrapSize = FALSE,
                            fittingMethod = NULL,
                            bootstrapFitcontrol = NULL,
                            covType = c("observedInform",
                                        "Fisher")
                            ) {
  list(
    alpha = alpha,
    trcount = trcount,
    bootType = if(missing(bootType)) "parametric" else bootType,
    B = B,
    confType = if(missing(confType)) "percentilic" else confType,
    keepbootStat = keepbootStat,
    traceBootstrapSize = traceBootstrapSize,
    fittingMethod = fittingMethod,
    bootstrapFitcontrol = bootstrapFitcontrol,
    covType = if (missing(covType)) "observedInform" else covType
  )
}
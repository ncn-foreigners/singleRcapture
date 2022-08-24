#' Control parameters for regression
#'
#' @param epsilon relative tolerance (absolute for robust method) for fitting algorithms by default 1e-8.
#' @param maxiter Maximal number of iterations.
#' @param verbose Value indicating whether to trace steps of fitting algorithm for robust fitting method different values of verbose give the following information:
#' \itemize{
#'   \item 1 -- Returns information on the number of current iteration and current log-likelihood.
#'   \item 2 -- Returns information on vector of regression parameters at current iteration (and all of the above).
#'   \item 3 -- Returns information on reduction of log-likelihood at current iteration (and all of the above).
#'   \item 4 -- Returns information on value of log-likelihood function gradient at current iteration (and all of the above).
#'   \item 5 -- Returns information on convergence criterion and values that are taken into account when considering convergence (and all of the above).
#' }
#' if mle method was chosen verbose will be passed to [stats::optim()] as trace.
#' @param start initial parameters for regression associated with main formula specified in function call if NULL they will be derived from simple poisson regression.
#' @param alphaStart initial parameters for dispersion parameter if applies.
#' @param omegaStart initial parameters for inflation parameter if applies.
#' @param piStart initial parameters for probability parameter if applies.
#' @param mleMethod method of [stats::optim()] used L-BFGS-B is the default except for negative binomial and one inflated models where Nelder-Mead is used.
#' @param silent logical indicating whether warnings in robust method should be suppressed.
#' @param optimPass optional parameter allowing for passing list into as stats::optim(..., control = optimPass) if FALSE then list of control parameters will be inferred from other parameters.
#' @param stepsize only for robust, scaling of stepsize lower value means slower convergence but more accuracy by default 1. In general if fitting algorithm fails lowering this value tends to be most effective at correcting it.
#' @param momentumFactor experimental parameter in robust only allowing for taking previous step into account at current step, i.e instead of updating regression parameters as:
#' \loadmathjax
#' \mjsdeqn{\boldsymbol{\beta}_{(a)} = \boldsymbol{\beta}_{(a-1)} + \text{stepsize} \cdot \text{step}_{(a)}}
#' the update will be made as:
#' \mjsdeqn{\boldsymbol{\beta}_{(a)} = \boldsymbol{\beta}_{(a-1)} + \text{stepsize} \cdot (\text{step}_{(a)} + \text{step}_{(a-1)})}
#' @param momentumActivation the value of log-likelihood reduction bellow which momentum will apply.
#'
#' @return A list with selected parameters, it is also possible to call list directly.
#' @export
control.method <- function(epsilon = 1e-8,
                           maxiter = 1000,
                           verbose = 0,
                           start = NULL,
                           alphaStart = NULL,
                           omegaStart = NULL,
                           piStart = NULL,
                           mleMethod = "L-BFGS-B",
                           silent = FALSE,
                           optimPass = FALSE,
                           stepsize = 1,
                           momentumFactor = 0,
                           momentumActivation = 5) {
  list(
    epsilon = epsilon,
    maxiter = maxiter,
    verbose = verbose,
    start = start,
    alphaStart = alphaStart,
    omegaStart = omegaStart,
    piStart = piStart,
    mleMethod = mleMethod,
    silent = silent,
    optimPass = optimPass,
    stepsize = stepsize,
    momentumFactor = momentumFactor,
    momentumActivation = momentumActivation
  )
}
#' control.model
#' TODO
#' @param weightsAsCounts TODO
#' @param omegaFormula TODO
#' @param alphaFormula TODO
#' for now does nothing
#' @return control.model
#' @export
control.model <- function(weightsAsCounts = FALSE,
                          omegaFormula = ~ 1,
                          alphaFormula = ~ 1,
                          piFormula = ~ .
                          ) {
  # TODO
  list(
    weightsAsCounts = weightsAsCounts,
    omegaFormula = omegaFormula,
    alphaFormula = alphaFormula,
    piFormula = piFormula
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
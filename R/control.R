#' Control parameters for regression
#'
#' @param epsilon Relative tolerance for fitting algorithms by default 1e-8.
#' @param maxiter Maximum number of iterations.
#' @param verbose Value indicating whether to trace steps of fitting algorithm for robust its either 0 (for no tracing), 1 (for tracing logarithm likelihood) or 2 (for tracing logarithm likelihood and vector of regression parameters) for mle it is passed to  stats::optim as trace.
#' @param start Initial parameters for regression coefficients if NULL they will be derived from simple poisson regression.
#' @param dispersionstart Initial parameters for dispersion parameter if applicable.
#' @param mleMethod Passed to stats::optim. L-BFGS-B is the default except for negative binomial models where Nelder-Mead is used.
#' @param silent Logical, indicating whether warnings in robust method should be suppressed.
#' @param optimPass Optional list of parameters passed to stats::optim(..., control = optimPass) if FALSE then list of control parameters will be inferred from other parameters.
#' @param disp.given Logical, indicates whether to estimate dispersion parameter or assume the one provided is already correct.
#' @param stepsize Only for robust, scaling of stepsize lower value means slower convergence but more accuracy by default 1.
#'
#' @return List with selected parameters, it is also possible to call list directly.
#' @export
control.method <- function(epsilon = 1e-8,
                           dispEpsilon = 1e-5,
                           maxiter = 1000,
                           verbose = 0,
                           start = NULL,
                           dispersionstart = NULL,
                           mleMethod = "L-BFGS-B",
                           silent = FALSE,
                           optimPass = FALSE,
                           disp.given = FALSE,
                           stepsize = 1) {
  list(
    epsilon = epsilon,
    dispEpsilon = dispEpsilon,
    maxiter = maxiter,
    verbose = verbose,
    start = start,
    dispersionstart = NULL,
    mleMethod = mleMethod,
    silent = silent,
    optimPass = optimPass,
    stepsize = stepsize
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
                          omegaFormula = NULL,# ????
                          alphaFormula = NULL # ????
                          # I suspect no other parameters will be used in the whole package so maybe its better to just specify them in control
                          # instead of making formula argument in main function a list or something like that.
                          # In VGAM there is a zero argument when calling vgam family class functions it specifies which argument
                          # is to be modeled as intercept only eg. in ztnegbin zero is by default on alpha so their fitting is 
                          # the same as ours. The only other aptions are for lambda or nothing i.e. there is no manual formula specification
                          ) {
  # TODO
  list(
    weightsAsCounts = weightsAsCounts,
    omegaFormula = omegaFormula,
    alphaFormula = alphaFormula
  )
}
#' Control parameters for population size estimation
#'
#' @param alpha Significance level, 0.05 used by default.
#' @param trcount Truncated count -- a number to be added to point estimator and both sides of confidence intervals.
#' @param bootType Bootstrap type. Default is \code{"parametric"}, other possible values are: \code{"semiparametric"} and \code{"nonparametric"}.
#' @param B Number of bootstrap samples to be performed (default 500).
#' @param confType Type of confidence interval for bootstrap confidence interval, \code{"percentile"} by default. Other possibility: \code{"studentized"} and \code{"basic"}.
#' @param keepbootStat Boolean value indicating whether to keep a vector of statistics produced by bootstrap.
#' @param traceBootstrapSize Boolean value indicating whether to print size of bootstrapped sample after truncation for semi- and fully parametric bootstraps.
#' @param fittingMethod Method used for fitting models from bootstrap samples.
#' @param bootstrapFitcontrol Control parameters for each regression works exactly like \code{control.method}
#' @param covType Type of covariance matrix for regression parameters by default observed information matrix, more options will be here in the future.
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
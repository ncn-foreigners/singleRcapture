#' Control parameters for regression
#'
#' @param epsilon relative tolerance for fitting algorithms by default 1e-8
#' @param maxiter Maximal number of iterations
#' @param verbose Value indicating whether to trace steps of fitting algorithm for robust its either 0 (for no tracing), 1 (for tracing logarithm likelihhod) or 2 (for tracing logarithm likelihood and vector of regression parameters) for mle it is passed to  stats::optim as trace
#' @param start initial parameters for regression if NULL they will be derived from simple poisson regression
#' @param dispersionstart initial parameters for dispersion parameter if applies
#' @param mleMethod method of stats::optim used L-BFGS-B is the default except for negative binomial models where Nelder-Mead is used
#' @param silent logical indicating whether warnings in robust method should be suppressed
#' @param optimPass optional parameter allowing for passing list into as stats::optim(..., control = optimPass) if FALSE then list of controlparameters will be infered from other parameters
#' @param disp.given logical indicates whether to estimate dispersion parameter or assume the one provided is already correct
#'
#' @return A list with selected parameters, it is also possible to call list directly.
#' @export
control.method <- function(epsilon = 1e-8,
                           maxiter = 1000,
                           verbose = 0,
                           start = NULL,
                           dispersionstart = NULL,
                           mleMethod = "L-BFGS-B",
                           silent = FALSE,
                           optimPass = FALSE,
                           disp.given = FALSE) {
  list(
    epsilon = epsilon,
    maxiter = maxiter,
    verbose = verbose,
    start = start,
    dispersionstart = NULL,
    mleMethod = mleMethod,
    silent = silent,
    optimPass = optimPass
  )
}
control.model <- function(typefitted = "link") {
  # TODO
  list(typefitted = typefitted)
}
#' Control parameters for variance estimation
#'
#' @param alpha Significance level, a number from frange (0, 1) 5% by default.
#' @param trcount Truncated count - a number to be added to point estimator and both sides of confidence intervals.
#' @param bootType Type of bootstrap performed, by default Parametric, other posible values are: Semiparametric and Nonparametric.
#' @param B Number of bootstrap samples to be performed.
#' @param confType Type of confidence interval for bootstrap confidence interval, Percentile by default, may be change to studentized.
#' @param keepbootStat Boolean value indicating whether to keep a vector of statistic produced by bootstrap, for large balues of strapNumber it may significant amount of size.
#' @param traceBootstrapSize Boolean value indicating whether to print size of bootstrapped sample after truncation for semi and fully parametric boostraps.
#' @param fittingMethod method used for fitting models from boostrap samples
#' @param bootstrapFitcontrol control parameters for each regression works exactly like control.method
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
                                         "basic"),
                            keepbootStat = TRUE,
                            traceBootstrapSize = FALSE,
                            fittingMethod = NULL,
                            bootstrapFitcontrol = NULL
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
    bootstrapFitcontrol = bootstrapFitcontrol
  )
}
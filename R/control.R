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
#'
#' @return List with selected parameters, it is also possible to call list directly.
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
#' @param alpha Significance level, 0.05 used by default.
#' @param trcount Truncated count -- a number to be added to point estimator and both sides of confidence intervals.
#' @param bootType Bootstrap type. Default is \code{"parametric"}, other possible values are: \code{"semiparametric"} and \code{"nonparametric"}.
#' @param B Number of bootstrap samples to be performed (default 500).
#' @param confType Type of confidence interval for bootstrap confidence interval, \code{"percentile"} by default. Other possibility: \code{"studentized"} and \code{"basic"}.
#' @param keepbootStat Boolean value indicating whether to keep a vector of statistics produced by bootstrap.
#' @param traceBootstrapSize Boolean value indicating whether to print size of bootstrapped sample after truncation for semi- and fully parametric bootstraps.
#' @param fittingMethod Method used for fitting models from bootstrap samples.
#' @param bootstrapFitcontrol Control parameters for each regression works exactly like \code{control.method}
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
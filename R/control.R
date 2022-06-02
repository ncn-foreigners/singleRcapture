control.method <- function() {
  # TODO
  list()
}
control.model <- function(typefitted = "link") {
  # TODO
  list(typefitted = typefitted)
}
#' Control parameters for variance estimation
#'
#' @param signiflevel Significance level, a number from frange (0, 1) 5% by default.
#' @param trcount Truncated count - a number to be added to point estimator and both sides of confidence intervals.
#' @param bootType Type of bootstrap performed, by default Parametric, other posible values are: Semiparametric and Nonparametric.
#' @param strapNumber Number of bootstrap samples to be performed.
#' @param confType Type of confidence interval for bootstrap confidence interval, Percentile by default, may be change to studentized.
#' @param keepbootStat Boolean value indicating whether to keep a vector of statistic produced by bootstrap, for large balues of strapNumber it may significant amount of size.
#' @param traceBootstrapSize Boolean value indicating whether to print size of bootstrapped sample after truncation for semi and fully parametric boostraps.
#'
#' @return A list with selected parameters, it is also possible to call list directly
#' @export
control.pop.var <- function(signiflevel = .05,
                            trcount = 0,
                            bootType = c("parametric",
                                         "semiparametric",
                                         "nonparametric"),
                            strapNumber = 500,
                            confType = c("percentilic",
                                         "studentized"),
                            keepbootStat = TRUE,
                            traceBootstrapSize = FALSE
                            #bootstrapFitcontrol = control.method, Add after control.method
                            ) {
  list(signiflevel = signiflevel,
       trcount = trcount,
       bootType = if(missing(bootType)) "parametric" else bootType,
       strapNumber = strapNumber,
       confType = if(missing(bootType)) "percentilic" else bootType,
       keepbootStat = keepbootStat,
       traceBootstrapSize = traceBootstrapSize)
}
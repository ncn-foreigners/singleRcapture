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
#' @param signiflevel Significance level, a number from frange (0, 1) 95% by default.
#' @param trcount Truncated count - a number to be added to point estimator and both sides of confidence intervals.
#' @param bootType Type of bootstrap performed, by default Parametric, other posible values are: Semiparametric and Nonparametric.
#' @param strapNumber Number of bootstrap resamplings to be performed.
#' @param confType a type of confidence interval for bootstrap confidence intervals, Percentilic by default, may be change to studentized.
#'
#' @return A list with selected parameters, it is also possible to call list directly
#' @export
control.pop.var <- function(signiflevel = .95,
                            trcount = 0,
                            bootType = "Parametric",
                            strapNumber = 500,
                            confType = "Percentilic") {
  list(signiflevel = signiflevel,
       trcount = trcount,
       bootType = bootType,
       strapNumber = strapNumber)
}
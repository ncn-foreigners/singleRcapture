#' Total Population Point and Interval estimate
#'
#' Creates a Horvitz-Thompson/Chao/Zelterman
#' point and interval estimate for total
#' and missing population from zero truncated poisson model.
#'
#' @param y Observed values
#' @param X A matrix of covariates
#' @param grad A gradient of a model with respect to regression parameters
#' @param parameter An estimated Parameter for the model
#' @param beta Fitted regression values
#' @param family Model Family
#' @param weights If model is weighted weights of particular observations
#' @param hessian Hessian of a model
#' @param dispersion Estimated dispersion parameter
#' for truncated Negative binomial distributions
#' @param method A method of constructing confidence interval either analytic
#' to use formula for analytic CI or bootstrap where bootstraped confidence
#' interval may either be based on 2.5%-97.5% percientiles ("bootstrapPerc")
#' or by estimating SD ("bootstrapSD")
#' @param trcount Optional parameter for Zero-one truncated models, if population estimate
#' is for one inflated model then it specifies one counts and includes them in
#' final population estimate both point and interval, and for zeltermann/chao
#' estimator where it specifies counts of not used in estimate
#'
#' @return Returns a list of size 3 with:
#' Point estimate, Interval estimate, and Variance
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom stats quantile
#' @export
#'
populationEstimate <- function(y,
                               X,
                               grad,
                               parameter,
                               beta,
                               weights = 1,
                               hessian,
                               family,
                               trcount,
                               dispersion,
                               method = "analytic") {
  if (method == "analytic") {
    strappedStatistic <- "No bootstrap performed"
    
    N <- family$pointEst(disp = dispersion,
                          pw = weights,
                          lambda = parameter) + trcount

    variation <- family$popVar(beta = beta, pw = weights,
                                lambda = parameter,
                                disp = dispersion,
                                hess = hessian, X = X)

    confidenceInterval <- c(lowerBound = max(N - 1.96 * sqrt(variation),
                                            (length(y) + trcount)),
                             upperBound = N + 1.96 * sqrt(variation))
  } else if (grepl("bootstrap", method, fixed = TRUE)) {
    
    N <- family$pointEst(disp = dispersion,
                         pw = weights,
                         lambda = parameter) + trcount
    
    if (!is.null(dispersion)) {
      beta <- beta[-1]
    }
    
    strappedStatistic <- noparBoot(family = family,
                                   y = y, X = X,
                                   dispersion = dispersion,
                                   beta = beta,
                                   weights = weights,
                                   trcount = trcount,
                                   numboot = 10000)

    if (method == "bootstrapSD") {
      variation <- stats::var(strappedStatistic)
      confidenceInterval <- c(lowerBound = max(N - 1.96 * sqrt(variation),
                                               (length(y) + trcount)),
                              upperBound = N + 1.96 * sqrt(variation))

    } else {
      variation <- stats::var(strappedStatistic)
      confidenceInterval <- stats::quantile(strappedStatistic,
                                            c(0.025, 0.975))
      names(confidenceInterval) <- c("lowerBound", "upperBound")
    }
  }

  list(pointEstimate = N,
       variance = variation,
       confidenceInterval = confidenceInterval,
       boot = strappedStatistic)
}

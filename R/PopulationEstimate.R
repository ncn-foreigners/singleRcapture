#' Function to estimate population size based on single-source capture-recapture models
#'
#' Provides point and interval estimates for the population size based on variety of models
#'
#' @param y Observed values
#' @param X A matrix of covariates
#' @param grad A gradient of a model with respect to regression parameters
#' @param lambda An estimated lambda parameter for the model 
#' @param beta Estimated parameters of the regression model (optional)
#' @param family Model family
#' @param weights Case weights if applied
#' @param weights0 Weights for all observations
#' @param hessian Hessian of a model
#' @param dispersion Estimated dispersion parameter
#' for truncated Negative binomial distributions
#' @param method A method of constructing confidence interval either analytic
#' to use formula for analytic CI or bootstrap where bootstraped confidence
#' interval may either be based on 2.5%-97.5% percentiles or by bootstrap SE
#' estimation
#' @param control List argument, same as control.pop.var in estimate_popsize
#'
#' @return Returns a list of size 3 with:
#' Point estimate, Interval estimate, and Variance
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom stats quantile
#' @importFrom stats qnorm
#' @export
populationEstimate <- function(y,
                               X,
                               grad,
                               lambda,
                               beta,
                               weights = 1,
                               weights0 = NULL,
                               hessian,
                               family,
                               dispersion,
                               method = "analytic",
                               control) {
  siglevel <- control$alpha
  trcount <- control$trcount
  numboot <- control$B
  sc <- qnorm(p = 1 - siglevel / 2)
  funBoot <- switch(control$bootType,
                    "parametric" = parBoot,
                    "semiparametric" = semparBoot,
                    "nonparametric" = noparBoot)
  if (method == "analytic") {
    strappedStatistic <- "No bootstrap performed"

    N <- family$pointEst(disp = dispersion,
                         pw = if (family$family == "zelterman") weights0 else weights,
                         lambda = lambda) + trcount

    variation <- as.numeric(family$popVar(beta = beta, 
                                          pw = if (family$family == "zelterman") weights0 else weights,
                                          lambda = lambda,
                                          disp = dispersion,
                                          hess = hessian, X = X))

    G <- exp(sc * sqrt(log(1 + variation / ((N - length(y)) ** 2))))
    confidenceInterval <- data.frame(t(data.frame(
      "Studentized" = c(lowerBound = max(N - sc * sqrt(variation), 
                                         length(y) + trcount), 
                        upperBound = N + sc * sqrt(variation)),
      "Logtransform" = c(lowerBound = length(y) + (N - length(y)) / G, 
                         upperBound = length(y) + (N - length(y)) * G)
    )))
  } else if (grepl("bootstrap", method, fixed = TRUE)) {
    N <- family$pointEst(disp = dispersion,
                         pw = if (family$family != "zelterman") {weights} else {weights0},
                         lambda = lambda) + trcount

    if (!is.null(dispersion)) {
      beta <- beta[-1]
    }

    strappedStatistic <- funBoot(family = family,
                                 y = y, 
                                 X = X,
                                 dispersion = dispersion,
                                 beta = beta,
                                 weights = list(weights, weights0),
                                 trcount = trcount,
                                 numboot = numboot,
                                 lambda = lambda,
                                 trace = control$traceBootstrapSize,
                                 method = control$fittingMethod,
                                 control.bootstrap.method = control$bootstrapFitcontrol)

    if (control$confType == "percentilic") {
      variation <- stats::var(strappedStatistic)
      confidenceInterval <- stats::quantile(strappedStatistic,
                                            c(siglevel / 2,
                                              1 - siglevel / 2))
      names(confidenceInterval) <- c("lowerBound", "upperBound")
    } else {
      variation <- stats::var(strappedStatistic)
      G <- exp(sc * sqrt(log(1 + variation / ((N - length(y)) ** 2))))
      
      confidenceInterval <- data.frame(t(data.frame(
        "Studentized" = c(lowerBound = max(N - sc * sqrt(variation), 
                                           (length(y) + trcount)), 
                          upperBound = N + sc * sqrt(variation)),
        "Logtransform" = c(lowerBound = length(y) + (N - length(y)) / G, 
                           upperBound = length(y) + (N - length(y)) * G)
      )))
    }
  }

  list(pointEstimate = N,
       variance = variation,
       confidenceInterval = confidenceInterval,
       boot = if (isTRUE(control$keepbootStat)) strappedStatistic else NULL,
       control = control)
}

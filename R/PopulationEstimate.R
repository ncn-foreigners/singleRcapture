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
#' @param control list, same as control.pop.var in estimate_popsize
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
                               parameter,
                               beta,
                               weights = 1,
                               hessian,
                               family,
                               dispersion,
                               method = "analytic",
                               control) {
  siglevel <- control$signiflevel
  boot.type <- control$bootType
  trcount <- control$trcount
  numboot <- control$strapNumber
  sc <- qnorm(p = 1 - (1 - siglevel) / 2)
  funBoot <- ifelse(boot.type == "Parametric",
                    parBoot,
                    ifelse(boot.type == "Semiparametric",
                           semparBoot,
                           noparBoot))
  if (method == "analytic") {
    strappedStatistic <- "No bootstrap performed"
    
    N <- family$pointEst(disp = dispersion,
                          pw = weights,
                          lambda = parameter) + trcount

    variation <- as.numeric(family$popVar(beta = beta, pw = weights,
                                          lambda = parameter,
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
                         pw = weights,
                         lambda = parameter) + trcount
    
    if (!is.null(dispersion)) {
      beta <- beta[-1]
    }

    if (is.null(numboot)){
      numboot <- ifelse(grepl("negbin", family$family, fixed = TRUE),
                        2000,
                        10000)
    }
    strappedStatistic <- funBoot(family = family,
                                 y = y, X = X,
                                 dispersion = dispersion,
                                 beta = beta,
                                 weights = weights,
                                 trcount = trcount,
                                 numboot = numboot,
                                 lambda = parameter)

    if (control$confType == "Percentilic") {
      variation <- stats::var(strappedStatistic)
      confidenceInterval <- stats::quantile(strappedStatistic,
                                            c((1 - siglevel) / 2, 
                                              1 - (1 - siglevel) / 2))
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
       boot = strappedStatistic,
       control = control)
}

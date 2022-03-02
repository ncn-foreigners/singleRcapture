#' summary.Model
#'
#' Describes an object of model class representing generalized linear regression
#' @param object Object to be summarised
#' @param ... other arguments to be passed to other methods
#' @return Easy to read summary of regression and most important data from a 'SingleR' class
#' @export
summary.singleR <- function(object, ...) {
  ob <- cbind(round(object$coefficients, digits = 3),
              round(object$standard_errors, digits = 3),
              round(object$wValues, digits = 2),
              round(object$pValues, digits = 3))
  colnames(ob) <- c('Coefficients','StandardErrors','W-values','P(>|t|)')

  print(object$call)
  cat("\nResponse Residuals:\n")
  print(summary(c(object$residuals)))
  cat("\nCoefficients:\n")
  print(ob)
  cat('\nAkane information criterion', object$aic,
      '\n\nLog-likelihood:', object$logL, 'on', object$df.residual, 'degrees of freedom\n',
      '\nPop size estimate:',
      '\nPoint estimate', object$populationSize$pointEstimate,
      '\nVariance', object$populationSize$variance,
      '\nConficence interval', object$populationSize$confidenceInterval)
}

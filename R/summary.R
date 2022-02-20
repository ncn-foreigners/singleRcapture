#' summary.Model
#'
#' Describes an object of model class representing generalized linear regression
#' @param object Object to be summarised
#' @param ... other arguments to be passed to other methods
#' @return Easy to read summary of regression and most important data from a 'SingleR' class
#' @export
summary.SingleR <- function(object, ...) {
  ob <- cbind(object$coefficients,
              object$Standard_Errors,
              signif(object$Tvalues, digits = 2),
              signif(object$Pvalues, digits = 2))
  colnames(ob) <- c('Coefficients','StandardErrors','Tvalues','P(>|t|)')

  print(object$call)
  cat("\nResponse Residuals:\n")
  #print("Response Residuals:\n")
  print(summary(c(object$residuals)))
  cat("\nCoefficients:\n")
  print(ob)
  cat('\nAkane information criterion', object$aic,
      '\n\nLog-likelihood:', object$LogL, 'on', object$df.residual, 'degrees of freedom\n',
      '\nPop size estimate:',
      '\nPoint estimate', object$Population_Size$Point_estimate,
      '\nVariance', object$Population_Size$Variance,
      '\nConficence interval', object$Population_Size$Confidence_Interval)
}

#' summary.Model
#'
#' @details Describes an object of model class representing generalized linear regression
#' @param object Object to be summarised
#' @param ... other arguments to be passed to other methods
#' @return Easy to read summary of regression and most important data from a 'SingleR' class
#' @export
summary.singleR <- function(object, ...) {
  signif <- c()
  for (k in object$pValues) {
    if (k <= 2e-16) {
      signif <- c(signif, "****")
    } else if (k <= .001) {
      signif <- c(signif, "***")
    } else if (k <= .01) {
      signif <- c(signif, "**")
    } else if (k <= .05){
      signif <- c(signif, "*")
    } else if (k <= .1){
      signif <- c(signif, ".")
    } else {
      signif <- c(signif, "")
    }
  }
  ob <- data.frame(round(object$coefficients, digits = 3),
                   round(object$standard_errors, digits = 3),
                   round(object$wValues, digits = 2),
                   signif(object$pValues, digits = 2),
                   signif)
  colnames(ob) <- c("Estimate", "Std. Error", "W-values", "P(>|t|)", "")

  print(object$call)
  cat("\nResponse Residuals:\n")
  print(summary(c(object$residuals)))
  cat("\nCoefficients:\n")
  print(ob)
  cat("\nSignif. codes:  0 \'****\' 0.001 \'***\' 0.01 \'**\' 0.05 \'*\' 0.1 \'.\' 1 \' \'")
  cat("\nAIC:", object$aic,
      "\n\nLog-likelihood:", object$logL, "on", object$df.residual,
      "Degrees of freedom",
      "\nNumber of iteration:", object$iter[1],
      "\n\nPopulation size estimation",
      "\nPoint estimate", object$populationSize$pointEstimate,
      "\nVariance", object$populationSize$variance,
      "\n95\% CI", object$populationSize$confidenceInterval)
}

#' Residuals of regresion
#'
#' @details S3 method for singleR class
#'
#' @param object TODO
#' @param type TODO
#' @param ... TODO
#'
#' @return returns a vector of residuals of
#' selected type
#' @export
residuals.singleR <- function(object,
                              type = c("pearson",
                                       "response",
                                       "working",
                                       "partial"),
                              ...) {
  # This is currently a placeholder to be written later
  return(object)
}

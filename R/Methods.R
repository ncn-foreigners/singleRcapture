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
  cat("-----------------------",
      "Signif. codes:  0 \'****\' 0.001 \'***\' 0.01 \'**\' 0.05 \'*\' 0.1 \'.\' 1 \' \'",
      sep = "\n")
  cat("\nAIC:", object$aic,
      "\nBIC:", object$bic,
      "\n\nLog-likelihood:", object$logL, "on", object$df.residual,
      "Degrees of freedom",
      "\nNumber of iteration:", object$iter[1],
      "\n\nPopulation size estimation results:",
      "\nPoint estimate", object$populationSize$pointEstimate,
      "\nVariance", object$populationSize$variance,
      "\n95% CI", object$populationSize$confidenceInterval)
}

#' residuals.singleR
#'
#' @details S3 method for singleR class
#'
#' @param object TODO
#' @param type TODO
#' @param ... TODO
#'
#' @return returns a vector of residuals of
#' selected type
#' @importFrom stats residuals
#' @export
residuals.singleR <- function(object,
                              type = c("pearson",
                                       "response",
                                       "working",
                                       "partial"),
                              ...) {
  res <- object$residuals
  disp <- object$dispersion
  wts <- object$prior.weights
  mu <- object$fitt.values
  y <- object$y
  if (type == "pearson") {
    rs <- res * sqrt(wts) / object$model$variance(mu, disp)
  } else if (type == "working") {
    rs <-  res / mu
  } else if (type == "response") {
    rs <- res
  } else {
    # partial residuals??
  }
  rs
}

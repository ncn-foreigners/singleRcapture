#' summary.Model
#'
#' Describes an object of model class representing generalized linear regression
#' @param object Object to be summarised
#' @return Easy to read summary of regression and most important data from a 'Model' class
#' @export
summary.Model <- function(object){
  ob <- cbind(object$Coefficients,
              object$StandardErrors,
              signif(object$Tvalues,digits = 3),
              signif(object$Pvalues,digits = 3))
  colnames(ob) <- c('Coefficients','StandardErrors','Tvalues','P(>|t|)')

  cat('Call:\nDependent variable =',object$DependentVariable,
      '  Expalantory variables =',object$ExplanatoryVariables,
      '\n\nCoefficients:\n')
  print(ob)
  cat('\nAkane information criterion',object$AIC,
      '\n\nResponse Residuals\n',summary(object$ResponseResiduals),
      '\n\nLog-likelihood:',object$LogL,'on ',object$df.reduced,'degrees of freedom\n',
      'Population size estimate:',
      '\nPoint estimate',object$PopulationSize$Point_estimate,
      '\nVariation',object$PopulationSize$Variation,
      '\nConficence interval',object$PopulationSize$ConfidenceInterval)
}

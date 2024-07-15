#' @method print singleRStaticCountData
#' @exportS3Method 
print.singleRStaticCountData <- function(x, ...) {
  cat("Call: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(coef(x))
  cat("\nDegrees of Freedom:", x$df.null, 
      "Total (i.e NULL);", x$dfResidual, 
      "Residual")
  
  cat("\nAIC: ", signif(AIC(x)), 
      "\nBIC: ", signif(BIC(x)), 
      "\nResidual deviance: ", signif(x$deviance))
  cat("\n-----------------------------------\n",
      "Population size estimation results:\n",
      sep = "")
  print(x$populationSize)

  invisible(x)
}

#' @method print singleRfamily
#' @exportS3Method 
print.singleRfamily <- function(x, hideFormulas = c(FALSE, TRUE), ...) {
  cat("Family of distributions:", x$family,
      "\nNames of parameters:", x$etaNames, 
      "\nLinks:", attr(x$links, "linkNames"), 
      sep = " ")
  if (!is.null(x$extraInfo)) {
    
    nn <- c("mean", "variance", "popSizeEst", "meanTr", "varianceTr")
    if (any(xx <- is.na(x$extraInfo[nn])))
      warning(cat("The: ", nn[is.na(xx)], sep = " ",
                  "slots are empty for family: ", x$family))
    
    if (!hideFormulas[1])
      cat("\n--\nFormula for mean in this distribution: ", x$extraInfo["mean"],
          "\nFormula for variance in this distribution: ", x$extraInfo["variance"],
          "\nFormula for population size estimation in this distribution: sum(", 
          x$extraInfo["popSizeEst"], ")",
          sep = "")
    
    if (!hideFormulas[2])
      cat("\n\nFormula for mean in truncated distribution: ", x$extraInfo["meanTr"],
          "\nFormula for variance in truncated distribution: ", x$extraInfo["varianceTr"],
          sep = "")
  }
  invisible(x)
}
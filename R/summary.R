#' \loadmathjax
#' @title Summary statistics for model of singleRStaticCountData class.
#' 
#' @description
#' A \code{summary} method for \code{singleRStaticCountData} class
#' 
#' 
#' @details Works 
#' analogically to \code{summary.glm} but includes population size estimation 
#' results. If any additional statistics, such as confidence intervals for
#' coefficients or coefficient correlation, are specified they will be printed.
#' 
#' @param object object of singleRStaticCountData class.
#' @param test type of test for significance of parameters \code{"t"} for t-test 
#' and \code{"z"} for normal approximation of students t distribution, by 
#' default \code{"z"} is used if there are more than 30 degrees of freedom
#' and \code{"t"} is used in other cases.
#' @param resType type of residuals to summarize any value that is allowed in
#' \code{residuals.singleRStaticCountData} except for \code{"all"} is allowed. By default 
#' pearson residuals are used.
#' @param correlation logical value indicating whether correlation matrix should
#' be computed from covariance matrix by default \code{FALSE}.
#' @param confint logical value indicating whether confidence intervals for
#' regression parameters should be constructed. By default \code{FALSE}.
#' @param cov covariance matrix corresponding to regression parameters. 
#' It is possible to give \code{cov} argument as a function of \code{object}.
#' If not specified it will be constructed using \code{vcov.singleRStaticCountData} method.
#' (i.e using Cramer-Rao lower bound)
#' @param popSizeEst a \code{popSizeEstResults} class object.
#' If not specified population size estimation results will be drawn from
#' \code{object}. If any post-hoc procedures, such as sandwich covariance matrix 
#' estimation or bias reduction, were taken it is possible to include them in 
#' population size estimation results by calling \code{redoPopEstimation}.
#' @param ... additional optional arguments passed to the following functions:
#' \itemize{
#' \item \code{vcov.singleRStaticCountData} -- if no \code{cov} argument was provided.
#' \item \code{cov} -- if \code{cov} parameter specified at call was a function.
#' \item \code{confint.singleRStaticCountData} -- if \code{confint} parameter was set to \code{TRUE} at function call.
#' In particular it is possible to set confidence level in \code{...}.
#' }
#' 
#' @return An object of \code{summarysingleRStaticCountData} class containing:
#' \itemize{
#' \item \code{call} -- A call which created \code{object}.
#' \item \code{coefficients} -- A dataframe with estimated regression coefficients
#' and their summary statistics such as standard error Wald test statistic and
#' p value for Wald test.
#' \item \code{residuals} -- A vector of residuals of type specified at call.
#' \item \code{aic} -- Akaike's information criterion.
#' \item \code{bic} -- Bayesian (Schwarz's) information criterion.
#' \item \code{iter} -- Number of iterations taken in fitting regression.
#' \item \code{logL} -- Logarithm of likelihood function evaluated at coefficients.
#' \item \code{deviance} -- Residual deviance.
#' \item \code{populationSize} -- Object with population size estimation results.
#' \item \code{dfResidual} -- Residual degrees of freedom.
#' \item \code{sizeObserved} -- Size of observed population.
#' \item \code{correlation} -- Correlation matrix if \code{correlation} parameter was set to \code{TRUE}
#' \item \code{test} -- Type of statistical test performed.
#' \item \code{model} -- Family class object specified in call for \code{object}.
#' \item \code{skew} -- If bootstrap sample was saved contains estimate of skewness.
#' }
#' @method summary singleRStaticCountData
#' @importFrom stats pt
#' @importFrom stats coef
#' @importFrom stats sd
#' @seealso [redoPopEstimation()] [stats::summary.glm()]
#' @exportS3Method 
summary.singleRStaticCountData <- function(object, 
                                           test = c("t", "z"), 
                                           resType = "pearson", 
                                           correlation = FALSE, 
                                           confint = FALSE, 
                                           cov, 
                                           popSizeEst, 
                                           ...) {
  if (resType == "all") {stop("Can't use 'resType = all' in summary.singleRStaticCountData method, if you wish to obtain all aviable types of residuals call residuals.singleRStaticCountData method directly.")}
  dfResidual <- object$dfResidual
  if (missing(test)) {if (dfResidual > 30) test <- "z" else test <- "t"}
  if (missing(cov)) {
    cov <- vcov(object, ...)
  } else if (is.function(cov)) {
    cov <- cov(object, ...)
  }
  
  pers <- residuals(object, type = resType)
  
  cf <- stats::coef(object)
  se <- sqrt(diag(cov))
  
  wValues <- cf / se
  
  pValues <- switch (test,
    "t" = 2 *    stats::pt(q = -abs(wValues), df = dfResidual),
    "z" = 2 * stats::pnorm(q = abs(wValues), lower.tail = FALSE)
  )
  
  crr <- if (isFALSE(correlation)) {NULL} else {cov / outer(se, se)}
  
  if(isTRUE(correlation)) {rownames(crr) <- colnames(crr) <- names(cf)}
  
  cnfint <- if(isTRUE(confint)) {confint(object, ...)} else {NULL}
  
  if (is.numeric(object$populationSize$boot)) {
    n <- length(object$populationSize$boot)
    m <- sum((object$populationSize$boot - mean(object$populationSize$boot)) ^ 3) / n
    s <- sd(object$populationSize$boot)
    skew <- m / (s ^ 3)
  } else {
    skew <- NULL
  }
  
  ob <- data.frame(cf, se, wValues, pValues)
  
  colnames(ob) <- switch(test,
    "t" = c("Estimate", "Std. Error", "t value", "P(>|t|)"),
    "z" = c("Estimate", "Std. Error", "z value", "P(>|z|)")
  )
  if (isTRUE(confint)) {
    ob[, 4] <- cnfint[, 1]
    ob[, 5] <- cnfint[, 2]
    ob[, 6] <- pValues
    colnames(ob)[4:6] <- c(colnames(cnfint), colnames(ob)[4])
  }
  
  cf <- ob
  structure(
    list(
      call = object$call,
      coefficients = cf,
      residuals = pers,
      aic = AIC(object, ...),
      bic = BIC(object, ...),
      iter = object$iter,
      logL = object$logL,
      deviance = object$deviance,
      populationSize = if (missing(popSizeEst)) object$populationSize else popSizeEst,
      dfResidual = dfResidual,
      sizeObserved = object$sizeObserved,
      correlation = crr,
      test = test,
      model = object$model,
      skew = skew
    ),
    class = "summarysingleRStaticCountData"
  )
}
#' @method print summarysingleRStaticCountData
#' @importFrom stats printCoefmat
#' @exportS3Method 
print.summarysingleRStaticCountData <- function(x, 
                                                signif.stars = getOption("show.signif.stars"), 
                                                digits = max(3L, getOption("digits") - 3L), 
                                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Pearson Residuals:\n")
  print(summary(c(x$residuals[, 1])))
  
  cat("\nCoefficients:\n")
  
  cond <- matrix(sapply(
    x$model$etaNames, 
    FUN = function(k) {
      sapply(strsplit(rownames(x$coefficients), 
                      split = ":"), 
             FUN = function(x) x[[length(x)]] == k)
    }
  ), 
  ncol = length(x$model$etaNames), 
  dimnames = list(
    NULL,
    x$model$etaNames
    )
  )
  
  for (k in x$model$etaNames) {
    if (!all(cond[, k] == FALSE)) {
      cat("-----------------------\nFor linear predictors associated with:", k, "\n")
      # Base library has no str_length and no proper sub_string this is a 
      # slightly convoluted way of making do without them
      toPrint <- subset(x$coefficients, cond[,k])
      lengths <- sapply(rownames(toPrint), 
                        function(x) {length(strsplit(x, split = "")[[1]])})
      lK <- length(unlist(strsplit(k, split = "")))
      rownames(toPrint) <- sapply(
        1:nrow(toPrint), 
        function(x) {substr(
          x = rownames(toPrint)[x], 
          start = 1, 
          stop = lengths[[x]] - (1 + lK)
        )}
      )
      printCoefmat(
        toPrint, digits = digits, 
        signif.stars = signif.stars, 
        signif.legend = if (k == x$model$etaNames[length(x$model$etaNames)]) signif.stars else FALSE, 
        P.values = TRUE, has.Pvalue = TRUE, 
        na.print = "NA", 
        ...
      )
    } else {
      cat("-----------------------\nFor linear predictors associated with:", k, "\n")
      printCoefmat(
        subset(x$coefficients, rowSums(cond) == 0), 
        digits = digits, 
        signif.stars = signif.stars, 
        signif.legend = if (k == x$model$etaNames[length(x$model$etaNames)]) signif.stars else FALSE, 
        P.values = TRUE, has.Pvalue = TRUE, 
        na.print = "NA", 
        ...
      )
    }
  }
  
  cat("\n")
  
  if (!is.null(x$correlation)) {
    corr <- round(x$correlation, digits = 2)
    corr[!lower.tri(corr)] <- ""
    print(corr[-1, -dim(corr)[2]], quote = FALSE)
    cat("\n")
  }
  
  sd <- sqrt(x$populationSize$variance)
  
  if (x$populationSize$control$sd == "normalMVUE") {
    sd <- sd / (sqrt(2 / (x$sizeObserved - 1)) * 
    exp(lgamma(x$sizeObserved / 2) - lgamma((x$sizeObserved - 1) / 2)))
  }
  
  cat(
    "AIC: ", x$aic,
    "\nBIC: ", x$bic,
    "\nResidual deviance: ", x$deviance,
    "\n\nLog-likelihood: ", x$logL, " on ", x$dfResidual, " Degrees of freedom ",
    # optim does not allow for accessing information
    # on number of iterations performed only a number 
    # of calls for gradient and objective function
    if (isTRUE(x$call$method == "optim")) "\nNumber of calls to log-likelihood function: "
    else "\nNumber of iterations: " , x$iter[1], 
    "\n-----------------------",
    "\nPopulation size estimation results: ",
    "\nPoint estimate ", x$populationSize$pointEstimate, 
    "\nObserved proportion: ", round(100 * x$sizeObserved / x$populationSize$pointEstimate, 
                                     digits = 1), 
    "% (N obs = ", x$sizeObserved, ")",
    if (!is.null(x$skew)) "\nBoostrap sample skewness: ", 
    if (!is.null(x$skew)) x$skew, 
    if (!is.null(x$skew)) "\n0 skewness is expected for normally distributed variable\n---",
    if (isTRUE(x$call$popVar == "bootstrap")) "\nBootstrap Std. Error " else "\nStd. Error ", 
    sd, "\n", 
    (1 - x$populationSize$control$alpha) * 100, 
    "% CI for the population size:\n", 
    sep = ""
  )
  
  print(x$populationSize$confidenceInterval)
  
  cat((1 - x$populationSize$control$alpha) * 100, 
      "% CI for the share of observed population:\n", 
      sep = "")
  
  dd <- as.data.frame(x$populationSize$confidenceInterval)
  
  if (ncol(dd) == 1) {
    vctpop <- sort(x$populationSize$confidenceInterval, decreasing = TRUE)
    names(vctpop) <- rev(names(vctpop))
    print(100 * x$sizeObserved / vctpop)
  } else {
    print(data.frame(
      lowerBound = 100 * x$sizeObserved / dd[, 2], 
      upperBound = 100 * x$sizeObserved / dd[, 1],
      row.names = rownames(dd)
    ))
  }
  
  
  invisible(x)
}
#' summary.Model
#'
#' @details Describes an object of model class representing generalized linear regression
#' @param object Object for which method is applied summarized
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
  colnames(ob) <- c("Estimate", "Std. Error", "z value", "P(>|z|)", "")

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
      "\nNumber of iterations:", object$iter[1],
      "\n\nPopulation size estimation results:",
      "\nPoint estimate", object$populationSize$pointEstimate,
      "\nVariance", object$populationSize$variance,
      "\nStd. Error", sqrt(object$populationSize$variance),
      "\n", object$populationSize$control$signiflevel * 100, "% CI:\n")
  print(object$populationSize$confidenceInterval)
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
                                       "partial",
                                       "deviance"),
                              ...) {
  type <- match.arg(type)
  res <- object$residuals
  disp <- object$dispersion
  wts <- object$weights
  mu <- object$fitt.values
  y <- object$y
  rs <- switch(type,
               working = res / mu,
               response = res,
               raw = res,
               pearson = res * sqrt(wts / object$model$variance(mu = mu, 
                                                                disp = disp)),
               deviance = object$model$dev.resids(y = y, mu = mu,
                                                  disp = disp,
                                                  wt = wts))
  
  #if (type == "pearson") {
  #  rs <- res * sqrt(wts / object$model$variance(mu, disp))
  #} else if (type == "working") {
  #  rs <-  res / mu
  #} else if (type == "response") {
  #  rs <- res
  #} else if (type == "raw") {
  #  rs <- res
  #} else if (type == "deviance") {
     # partial residuals??
  #  rs <- sqrt(d) * sign(y - mu)
  #  stop()
  #}
  rs
}
#' Summary for marginal frequencies
#'
#' @param object object of singleRmargin class
#' @param df degrees of freedom are sometimes not possible to automatically obtain if so this overwrites df of an object
#' @param dropl5 boolean value indicating whether to group bins with frequencies < 5, drop them or do nothing
#' @param ... Currently does nothing
#'
#' @return A chi squared test for comparison between fitted and observed marginal frequencies
#' @export
summary.singleRmargin <- function(object, df = NULL,
                                  dropl5 = c("drop", 
                                             "group", 
                                             "no"), 
                                  ...) {
  if (length(dropl5) > 1) {dropl5 <- "no"}
  y <- object$y
  A <- object$table[names(y)]
  if ((is.null(df)) && (object$df < 1)) {
    warning("Degrees of freedom may be inacurate")
    df <- 1
  } else if (is.null(df)) {
    df <- object$df
  }
  if(dropl5 == "group") {
    l <- ((y < 5) | (A < 5))
    if (object$df == df) {df <- df - length(y) + length(y[!l]) + 1}
    y <- c(y[!l], sum(y[l]))
    A <- c(A[!l], sum(A[l]))
  } else if(dropl5 == "drop") {
    l <- ((y < 5) | (A < 5))
    if (object$df == df) {df <- df - length(y) + length(y[!l])}
    y <- y[!l]
    A <- A[!l]
  }
  X2 <- sum(((A - y) ** 2) / A)
  G <- 2 * sum(y * log(y / A))
  pval <- stats::pchisq(q = c(X2, G), df = df, lower.tail = FALSE)
  vect <- data.frame(round(c(X2, G), digits = 2),
                     rep(df, 2), signif(pval, digits = 2))
  rownames(vect) <- c("Chi-squared test", "G-test")
  colnames(vect) <- c("Test statistics", "df", "P(>X^2)")
  vect
}

#' vcov method for singleR class
#' @title vcov method for singleR class
#' @param object object of clas singleRclass
#' @param ... variables to pass to other methods
#'
#' @method vcov singleR
#' @return A covariance matrix for fitted coefficients obtained by inverting 
#' analitical hesian at estimated coefficients, i.e. using Cramér–Rao bound
#' with observed information matrix.
#' @export
vcov.singleR <- function(object, ...) {
  solve(
    object$model$make_hessian(y = object$y, X = object$X, 
                              weight = object$prior.weights)(object$coefficients)
  )
}
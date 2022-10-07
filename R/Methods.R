#' dfpopsize
#'
#' @param model Model for which leave one out diagnostic of popsize will be done.
#' @param dfbeta If dfbeta was already obtained it is possible to pass them into 
#' function so that they need not be computed for the second time.
#' @param observedPop TODO
#' @param ... Arguments to be passed down to other methods such as dfbeta.
#'
#' @return TODO: MAKE BETTER DOCUMENTATION
#' @export
dfpopsize <- function(model, dfbeta = NULL, observedPop = FALSE, ...) {
  UseMethod("dfpopsize")
}

# gridSearch.singleR <- function(object, useOptim, updateCalculations = TRUE, ...) {
#   # TODO::
#   ################
#   ## Important ###
#   ################
#   # check this for new function
#   # Make this a method or a function??
#   # TODO search on all start parameters
#   if (isFALSE(grepl("^ztoi", object$model$family))) {
#     stop("For now only one iflated models have this method implemented")
#   }
#   logL <- object$model$makeMinusLogLike(y = object$y, X = object$X, weight = object$prior.weights)
#   FUN <- function(CHNG) {
#     A <- estimate_popsize.fit(
#       y = object$y,
#       X = object$X,
#       family = object$model,
#       control = object$control$control.method,
#       method = object$call$method,
#       prior.weights = object$prior.weights,
#       start = c(CHNG, object$coefficients[-1]),
#       dispersion = object$dispersion,
#       omegaTheta = CHNG,
#       ...
#     )$beta
#     return(logL(A))
#   }
#   opt <- stats::optim(par = object$omegaTheta,fn = FUN, lower = -.95, upper = 10, method = "Brent", ...)$par
#   if (isTRUE(updateCalculations)) {
#     FITT <- estimate_popsize.fit(
#       y = object$y,
#       X = object$X,
#       family = object$model,
#       control = object$control$control.method,
#       method = object$call$method,
#       prior.weights = object$prior.weights,
#       start = c(opt, object$coefficients[-1]),
#       dispersion = object$dispersion,
#       omegaTheta = opt,
#       ...
#     )
#     object$coefficients <- FITT$beta
#     eta <- as.matrix(object$X) %*% object$coefficients[-1]
#     rownames(eta) <- rownames(object$X)
#     object$iter <- FITT$iter
#     object$weights <- FITT$weights
#     object$linear.predictors <- eta
#     object$fitt.values <- data.frame("mu" = object$model$mu.eta(eta, disp = object$dispersion, theta = opt),
#                                      "link" = object$model$linkinv(eta))
#     object$omegaTheta <- opt
#     LOG <- -logL(FITT$beta)
#     object$logL <- LOG
#     object$resRes <- object$prior.weights * (object$y - object$fitt.values)
#     object$aic <- 2 * (length(FITT$beta) - LOG)
#     object$bic <- length(FITT$beta) * log(length(object$y)) - 2 * LOG
#     object$deviance <- sum(object$model$dev.resids(y = object$y, 
#                                                    mu = object$model$linkinv(eta),
#                                                    disp = object$dispersion,
#                                                    wt = object$prior.weights, 
#                                                    theta = opt) ** 2)
#     object$populationSize <- singleRcaptureinternalpopulationEstimate(
#       y = object$y,
#       X = object$X,
#       grad = object$model$makeGradient(y = object$y, X = object$X, weight = object$prior.weights),
#       hessian = object$model$makeHessian(y = object$y, X = object$X, weight = object$prior.weights),
#       method = object$call$pop.var,
#       weights = object$prior.weights,
#       weights0 = object$prior.weights,
#       lambda = object$fitt.values$link,
#       family = object$model,
#       dispersion = object$dispersion,
#       omegaTheta = object$omegaTheta,
#       beta = object$coefficients,
#       control = object$populationSize$control
#     )
#     return(object)
#   } else {
#     return(opt)
#   }
# }
#' Summary for marginal frequencies
#'
#' @param object object of singleRmargin class.
#' @param df degrees of freedom are sometimes not possible to automatically obtain if so this overwrites df of an object.
#' @param dropl5 a character indicating treatment of cells with frequencies < 5 either grouping them, droping or leaving them as is. Defaults to drop.
#' @param ... Currently does nothing
#'
#' @method summary singleRmargin
#' @return A chi squared test for comparison between fitted and observed marginal frequencies
#' @exportS3Method
summary.singleRmargin <- function(object, df = NULL,
                                  dropl5 = c("drop", 
                                             "group", 
                                             "no"), 
                                  ...) {
  if (missing(dropl5)) {dropl5 <- "drop"}
  y <- object$y
  if (grepl("zot", object$name) & (1 %in% names(y))) {y <- y[-1]}
  A <- object$table[names(y)]
  if ((is.null(df)) && (object$df < 1)) {
    warning("Degrees of freedom may be inacurate.")
    df <- 1
  } else if (is.null(df)) {
    df <- object$df
  }
  if(dropl5 == "group") {
    l <- (A < 5)
    y <- c(y[!l], sum(y[l]))
    A <- c(A[!l], sum(A[l]))
  } else if(dropl5 == "drop") {
    l <- (A < 5)
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
  structure(
    list(Test = vect,
         l5 = switch(dropl5,
                     drop = "dropped",
                     group = "grouped",
                     no = "preserved"),
         y = y),
    class = "summarysingleRmargin"
  )
}
#' vcov method for singleR class
#' @title Obtain Covariance Matrix from Fitted singleR class Object
#' @param object object of class singleRclass.
#' @param type type of estimate for variance covariance matrix for now either
#' expected (fisher) information matrix or observed information matrix.
#' @param ... variables to pass to solve.
#' @description Returns a estimated covariance matrix for model coefficients
#' calculated from analytic hessian or fisher information matrix. 
#' Covariance type is taken from control parameter that have been provided
#' on call that created object.
#' @details For now the covariance matrix is calculated 
#' @method vcov singleR
#' @return A covariance matrix for fitted coefficients, rows and columns of which 
#' correspond to parameters returned by \code{coef} method.
#' @exportS3Method
vcov.singleR <- function(object, type = c("Fisher", "observedInform"), ...) {
  if (missing(type) ){type <- object$populationSize$control$covType}
  X <- singleRinternalGetXvlmMatrix(X = subset(object$modelFrame, select = attr(object$modelFrame, "names")[-1], subset = object$which$reg), nPar = object$model$parNum, 
                                    formulas = object$formula, parNames = object$model$etaNames)
  hwm <- X[[2]]
  X <- X[[1]]
  res <- switch(
    type,
    "observedInform" = solve(
      -object$model$makeHessian(y = object$y[object$which$reg], X = X, lambdaPredNumber = hwm[1],
                                weight = object$prior.weights[object$which$reg])(object$coefficients),
      ...
    ),
    "Fisher" = {
      if (object$call$method == "robust") {W <- object$weights} else {W <- object$model$Wfun(prior = object$prior.weights[object$which$reg], eta = object$linear.predictors)};
      solve(
      singleRinternalMultiplyWeight(X = X, W = W, hwm = hwm) %*% X,
      ...
    )}
  )
  dimnames(res) <- list(names(object$coefficients), names(object$coefficients))
  res
}
#' Hat values for singleRclass
#' @title Hat values for singleRclass
#' @param model object of clas singleRclass
#' @param ... additional parameters to pass to other methods
#' @description TODO
#' 
#' @method hatvalues singleR
#' @importFrom stats hatvalues
#' @return TODO
#' @exportS3Method 
hatvalues.singleR <- function(model, ...) {
  X <- subset(model$modelFrame, select = attr(model$modelFrame, "names")[-1], subset = model$which$reg)
  #X <- model$modelFrame[model$which$reg, attr(model$modelFrame, "names")[-1]]
  X <- singleRinternalGetXvlmMatrix(X = X, nPar = model$model$parNum, formulas = model$formula, parNames = model$model$etaNames)
  hwm <- X[[2]]
  X <- X[[1]]
  if (model$call$method == "robust") {
    W <- model$weights
  } else {
    W <- model$model$Wfun(prior = model$prior.weights[model$which$reg], eta = if (model$model$family == "zelterman") model$linear.predictors[model$which$reg, ] else model$linear.predictors)
  }
  mlt <- singleRinternalMultiplyWeight(X = X, W = W, hwm = hwm)
  hatvalues <- diag(X %*% solve(mlt %*% X) %*% mlt)
  hatvalues <- matrix(hatvalues, ncol = model$model$parNum, dimnames = list(1:(length(hatvalues) / model$model$parNum), model$model$etaNames))
  hatvalues
}
#' dfbeta for singleRclass
#' @title TODO
#' @param model Fitted object of singleR class
#' @param maxit.new maximal number of iterations for regression
#' @param ... Arguments to pass to other methods
#' @description TODO
#' 
#' @return TODO
dfbetasingleR <- function(model,
                          maxit.new = 1,
                          ...) {
  # formula method removed since it doesn't give good results will reimplement if we find better formula
  X <- subset(model$modelFrame, subset = model$which$reg, select = attr(model$modelFrame, "names")[-1])
  y <- model$y[model$which$reg]
  cf <- model$coefficients
  pw <- model$prior.weights[model$which$reg]
  res <- matrix(nrow = nrow(X), ncol = length(cf))
  hwm <- singleRinternalGetXvlmMatrix(X = X, nPar = model$model$parNum, formulas = model$formula, parNames = model$model$etaNames)[[2]]
  for (k in 1:nrow(X)) {
    res[k, ] <- cf - estimate_popsize.fit(
      control = control.method(
        silent = TRUE, 
        start = cf,
        maxiter = maxit.new + 1,
        ...
      ),
      y = y[-k],
      X = singleRinternalGetXvlmMatrix(X = subset(X, rownames(X) != rownames(X)[k]), nPar = model$model$parNum, formulas = model$formula, parNames = model$model$etaNames)[[1]],
      start = cf,
      family = model$model,
      prior.weights = pw[-k],
      method = model$call$method,
      hwm = hwm
    )$beta
  }
  colnames(res) <- names(model$coefficients)
  res
}
#' @title Confidence Intervals for Model Parameters
#' 
#' @description A function that computes studentized confidence intervals
#' for model coefficients
#' 
#' @param object a fitted model object.
#' @param parm names of parameters for which confidence intervals are to be 
#' computed, if missing all parameters will be considered
#' @param level confidence level for intervals.
#' @param ... additional argument(s) for methods.
#' 
#' @method confint singleR
#' @return An object with named columns that include upper and 
#' lower limit of confidence intervals
#' @exportS3Method
confint.singleR <- function(object,
                            parm, 
                            level = 0.95, 
                            ...) {
  std <- sqrt(diag(vcov.singleR(object)))
  if (missing(parm)) {
    coef <- object$coefficients
  } else {
    coef <- object$coefficients[parm]
    std <- std[parm]
  }
  sc <- qnorm(p = 1 - (1 - level) / 2)
  res <- data.frame(coef - sc * std, coef + sc * std)
  colnames(res) <- c(paste0(100 * (1 - level) / 2, "%"),
                     paste0(100 * (1 - (1 - level) / 2), "%"))
  res
}

# There is no need for doccumenting the following methods:
#' @method residuals singleR
#' @importFrom stats residuals
#' @exportS3Method
residuals.singleR <- function(object,
                              type = c("pearson",
                                       "pearsonSTD",
                                       "response",
                                       "working",
                                       "deviance",
                                       "all"),
                              ...) {
  type <- match.arg(type)
  res <- object$residuals
  disp <- object$dispersion
  omegaTheta <- object$omegaTheta
  wts <- object$prior.weights
  mu <- object$fitt.values
  y <- object$y
  if (!(all(object$which$reg == object$which$est)) && type == "all") stop("type = all is not aviable for some models")
  if (type == "pearsonSTD" && object$model$parNum > 1) {stop("Standardized pearson residuals not yet implemented for models with multiple linear predictors")}
  rs <- switch(
    type,
    working = as.data.frame(object$model$funcZ(eta = if (object$model$family == "zelterman") object$linear.predictors[object$which$reg, ] else object$linear.predictors, weight = object$weights, y = y[object$which$reg]), col.names = paste0("working:", object$model$etaNames)),
    response = res,
    pearson = data.frame("pearson" = (if (object$model$family == "zelterman") res$mu[object$which$reg] else res$mu) / sqrt(object$model$variance(eta = if (object$model$family == "zelterman") object$linear.predictors[object$which$reg, ] else object$linear.predictors, type = "trunc"))),
    pearsonSTD = data.frame("pearsonSTD" = (if (object$model$family == "zelterman") res$mu[object$which$reg] else res$mu) / sqrt((1 - hatvalues(object)) * object$model$variance(eta = if (object$model$family == "zelterman") object$linear.predictors[object$which$reg, ] else object$linear.predictors, type = "trunc"))),
    deviance = data.frame("deviance" = object$model$dev.resids(y = y[object$which$reg], eta = object$linear.predictors, wt = wts[object$which$reg])),
    all = {colnames(res) <- c("muResponse", "linkResponse");
      data.frame(
      as.data.frame(object$model$funcZ(eta = if (object$model$family == "zelterman") object$linear.predictors[object$which$reg, ] else object$linear.predictors, weight = object$weights, y = y[object$which$reg]), col.names = paste0("working:", object$model$etaNames)),
      res,
      "pearson" = as.numeric((if (object$model$family == "zelterman") res$mu[object$which$reg] else res$mu) / sqrt(object$model$variance(eta = if (object$model$family == "zelterman") object$linear.predictors[object$which$reg, ] else object$linear.predictors, type = "trunc"))),
      "pearsonSTD" = if (object$model$parNum == 1) as.numeric((if (object$model$family == "zelterman") res$mu[object$which$reg] else res$mu) / sqrt((1 - hatvalues(object)) * object$model$variance(eta = if (object$model$family == "zelterman") object$linear.predictors[object$which$reg, ] else object$linear.predictors, type = "trunc"))) else 0,
      "deviance" = as.numeric(object$model$dev.resids(y = y[object$which$reg], eta = object$linear.predictors, wt = wts[object$which$reg])),
      row.names = rownames(object$linear.predictors)
    )}
  )
  rs
}
#' @method family singleR
#' @importFrom stats family
#' @exportS3Method
family.singleR <- function(object, ...) {
  object$model
}
#' @method print summarysingleRmargin
#' @exportS3Method 
print.summarysingleRmargin <- function(x, ...) {
  cat("Test for Goodness of fit of a regression model:\n",
      "\n", sep = "")
  print(x$Test)
  cat("\n--------------------------------------------------------------",
      "\nCells with fitted frequencies of < 5 have been", x$l5, "\nNames of cells used in calculating test(s) statistic:", names(x$y), "\n", sep = " ")
}
#' @method AIC singleR
#' @importFrom stats AIC
#' @exportS3Method 
AIC.singleR <- function(object, ...) {
  object$aic
}
#' @method BIC singleR
#' @importFrom stats BIC
#' @exportS3Method 
BIC.singleR <- function(object, ...) {
  object$bic
}
#' @method extractAIC singleR
#' @importFrom stats extractAIC
#' @exportS3Method 
extractAIC.singleR <- function(fit, scale, k = 2, ...) {
  -2 * fit$logL + k * length(fit$coefficients)
}
#' @method dfbeta singleR
#' @importFrom stats dfbeta
#' @exportS3Method 
dfbeta.singleR <- function(model, ...) {
  dfbetasingleR(model, ...)
}
#' @method logLik singleR
#' @importFrom stats logLik
#' @exportS3Method 
logLik.singleR <- function(object, ...) {
  val <- object$logL
  attr(val, "nobs") <- dim(residuals(object))[1]
  attr(val, "df") <- length(object$coefficients)
  class(val) <- "logLik"
  val
}
#' @method model.matrix singleR
#' @importFrom stats model.matrix
#' @exportS3Method 
model.matrix.singleR <- function(object, type = c("lm", "vlm"), ...) {
  if (missing(type)) type <- "lm"
  switch (type,
    lm = {
      X <- object$X;
      X[object$which$reg, ]
      },
    vlm = {
      X <- subset(object$modelFrame, subset = object$which$reg, select = attr(object$modelFrame, "names")[-1]);
      singleRinternalGetXvlmMatrix(X = X, nPar = object$model$parNum, formulas = object$formula, parNames = object$model$etaNames);
      }
  )
}

#' @method dfpopsize singleR
#' @exportS3Method 
dfpopsize.singleR <- function(model, dfbeta = NULL, observedPop = FALSE, ...) {
  if (isTRUE(model$call$pop.var == "bootstrap")) warning("dfpopsize may (in some cases) not work correctly when bootstrap was chosen as population variance estimate.")
  dfb <- if (is.null(dfbeta)) {dfbeta(model, ...)} else {dfbeta}
  if (model$model$family == "zelterman") {
    dfbnew <- matrix(0, ncol = ncol(dfb), nrow = model$sizeObserved)
    rownames(dfbnew) <- as.character(1:model$sizeObserved)
    dfbnew[rownames(dfb), ] <- dfb
    dfb <- dfbnew
  }
  X <- model$modelFrame[model$which$est, attr(model$modelFrame, "names")[-1]]
  X <- subset(model$modelFrame, subset = model$which$est, select = attr(model$modelFrame, "names")[-1])
  N <- model$populationSize$pointEstimate
  res <- NULL
  range <- 1:sum(model$which$est)
  pw <- model$prior.weights[model$which$est]
  for (k in range) {
    cf <- model$coefficients - dfb[k, ]
    res <- c(res, model$model$pointEst(
      eta = matrix(singleRinternalGetXvlmMatrix(X = subset(X, rownames(X) != rownames(X)[k]), nPar = model$model$parNum, formulas = model$formula, parNames = model$model$etaNames)[[1]] %*% cf, ncol = model$model$parNum),
      pw = pw[-k]) + model$trcount)
  }
  
  if(isTRUE(observedPop) & (grepl("zot", model$model$family) | model$model$family == "chao")) {
    res1 <- vector(mode = "numeric", length = model$sizeObserved)
    names(res1) <- 1:model$sizeObserved
    res1[model$which$est] <- N - res
    res1[!model$which$est] <- 1
    
  } else {
    res1 <- N - res
  }
  
  res1
}
#' @method summary singleR
#' @importFrom stats pt
#' @importFrom stats coef
#' @exportS3Method 
summary.singleR <- function(object, test = c("t", "z"), resType = "pearson", correlation = FALSE, confint = FALSE, ...) {
  if (resType == "all") {stop("Can't use 'resType = all' in summary.singleR method, if you wish to obtain all aviable types of residuals call residuals.singleR method directly.")}
  if (missing(test)) {test <- "z"}
  df.residual <- object$df.residual
  cov <- vcov.singleR(object, ...)
  pers <- residuals.singleR(object, type = resType, ...)
  cf <- object$coefficients
  se <- sqrt(diag(cov))
  wValues <- cf / se
  pValues <- switch (test,
                     "t" = 2 * stats::pt(q = -abs(wValues), df = df.residual),
                     "z" = 2 * stats::pnorm(q =  abs(wValues), lower.tail = FALSE)
  )
  crr <- if (isFALSE(correlation)) {NULL} else {cov / outer(se, se)}
  if(isTRUE(correlation)) {rownames(crr) <- colnames(crr) <- names(cf)}
  cnfint <- if(isTRUE(confint)) {confint.singleR(object, ...)} else {NULL}
  structure(
    list(
      call = object$call,
      coefficients = cf,
      standardErrors = se,
      wValues = wValues,
      pValues = pValues,
      residuals = pers,
      aic = object$aic,
      bic = object$bic,
      iter = object$iter,
      logL = object$logL,
      deviance = object$deviance,
      populationSize = object$populationSize,
      df.residual = df.residual,
      sizeObserved = object$sizeObserved,
      correlation = crr,
      test = test,
      cnfint = cnfint
    ),
    class = "summarysingleR"
  )
}
#' @importFrom stats cooks.distance
#' @method cooks.distance singleR
#' @exportS3Method 
cooks.distance.singleR <- function(model, ...) {
  if (model$model$parNum > 1) stop("Cooks distance is only implemented for single parameter families.")
  res <- residuals(model, type = "pearsonSTD") ** 2
  res <- res[, 1]
  res <- (res * (hatvalues(model) / (length(model$coefficients))))
  names(res) <- rownames(model$linear.predictors)
  res
}
#' @method print summarysingleR
#' @exportS3Method 
print.summarysingleR <- function(x, ...) {
  # ifelse is faster than if(_) {} else {}, sapply is faster than for hence the change
  signifCodes <- sapply(x$pValues, function(k) {
    ifelse(k <= 0, "****",
           ifelse(k <= .001, "***",
                  ifelse(k <= .01, "**", 
                         ifelse(k <= .05, "*",
                                ifelse(k <= .1, ".", "")))))
  })
  ob <- data.frame(round(x$coefficients, digits = 3),
                   round(x$standardErrors, digits = 3),
                   round(x$wValues, digits = 2),
                   signif(x$pValues, digits = 2),
                   signifCodes)
  colnames(ob) <- switch(x$test,
                         "t" = c("Estimate", "Std. Error", "t value", "P(>|t|)", ""),
                         "z" = c("Estimate", "Std. Error", "z value", "P(>|z|)", ""))
  if (!is.null(x$cnfint)) {
    ob[, 5] <- x$cnfint[, 1]
    ob[, 6] <- x$cnfint[, 2]
    ob[, 7] <- signifCodes
    colnames(ob)[5:7] <- c(colnames(x$cnfint), "")
  }
  
  print(x$call)
  cat("\nPearson Residuals:\n")
  print(summary(c(x$residuals[, 1])))
  cat("\nCoefficients:\n")
  print(ob)
  cat("-----------------------",
      "Signif. codes:  0 \'****\' 0.001 \'***\' 0.01 \'**\' 0.05 \'*\' 0.1 \'.\' 1 \' \'\n",
      sep = "\n")
  if (!is.null(x$correlation)) {
    corr <- round(x$correlation, digits = 2)
    corr[!lower.tri(corr)] <- ""
    print(corr[-1, -dim(corr)[2]], quote = FALSE)
    cat("\n")
  }
  sd <- sqrt(x$populationSize$variance)
  if (x$populationSize$control$sd == "normalMVUE") {
    sd <- sd / (sqrt(2 / (x$sizeObserved - 1)) * exp(lgamma(x$sizeObserved / 2) - lgamma((x$sizeObserved - 1) / 2)))
  }
  cat("AIC: ", x$aic,
      "\nBIC: ", x$bic,
      "\nDeviance: ", x$deviance,
      "\n\nLog-likelihood: ", x$logL, " on ", x$df.residual, " Degrees of freedom ",
      if (isTRUE(x$call$method == "robust")) {
        "\nNumber of iterations: "
      } else {
        "\nNumber of calls to log-likelihood function: " # optim does not allow for accesing information
        # on number of iterations performed only a number of calls for gradient and objective function
      }, x$iter[1], 
      "\n-----------------------",
      "\nPopulation size estimation results: ",
      "\nPoint estimate ", x$populationSize$pointEstimate, 
      "\nObserved proportion: ", round(100 * x$sizeObserved / x$populationSize$pointEstimate, digits = 1), "% (N obs = ", x$sizeObserved, ")",
      if (isTRUE(x$call$pop.var == "bootstrap")) {"\nBootstrap Std. Error "} else {"\nStd. Error "}, sd,
      "\n", (1 - x$populationSize$control$alpha) * 100, "% CI for the population size:\n", sep = "")
  print(x$populationSize$confidenceInterval)
  cat((1 - x$populationSize$control$alpha) * 100, "% CI for the share of observed population:\n", sep = "")
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
}

#' @importFrom stats fitted
#' @method fitted singleR
#' @exportS3Method 
fitted.singleR <- function(object,
                           type = c("mu", 
                                    "link"), # maybe add marginal frequencies/other options
                           ...) {
  if (missing(type)) type <- "mu"
  object$fitt.values[[type]] # fitted should return either E(Y) or E(Y|Y>0) otherwise we're breaking R conventions
}
# TODO:: update
#' #' simulate
#' #' 
#' #' An S3class for \code{stats::simulate} to handle \code{singleR} objects.
#' #' 
#' #' @param object an object representing a fitted model.
#' #' @param nsim number of response vectors to simulate. Defaults to \code{1}.
#' #' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' #' @param ... additional optional arguments.
#' #' @return a \code{data.frame} with \code{n} rows and \code{nsim} columns.
#' #' @seealso [stats::simulate()]
#' #' @examples 
#' #' set.seed(1)
#' #' N <- 10000
#' #' gender <- rbinom(N, 1, 0.2)
#' #' eta <- -1 + 0.5*gender
#' #' counts <- rpois(N, lambda = exp(eta))
#' #' df <- data.frame(gender, eta, counts)
#' #' df2 <- subset(df, counts > 0)
#' #' mod1 <-  estimate_popsize(formula = counts ~ 1 + gender, data = df2, 
#' #' model = "ztpoisson", method = "mle", pop.var = "analytic")
#' #' mod1_sims <- simulate(mod1, nsim=10)
#' #' colMeans(mod1_sims) 
#' #' @importFrom stats simulate
#' #' @method simulate singleR
#' #' @exportS3Method
#' simulate.singleR <- function(object, nsim=1, seed = NULL, ...) {
#'   ###################
#'   ### TODO update ###
#'   ###################
#'   stop("Not yet updated")
#'   if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
#'     runif(1)
#'   if (is.null(seed))
#'     RNGstate <- get(".Random.seed", envir = .GlobalEnv)
#'   else {
#'     R.seed <- get(".Random.seed", envir = .GlobalEnv)
#'     set.seed(seed)
#'     RNGstate <- structure(seed, kind = as.list(RNGkind()))
#'     on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
#'   }
#'   
#'   linpred <- object$linear.predictors
#'   n <- NROW(linpred)
#'   if (grepl("negbin", object$model$family)) {
#'     val <- object$model$simulate(n*nsim, exp(linpred), exp(-object$dispersion)) # dispersion argument will be removed in next update remember to change that to accomodate that
#'   } else {
#'     val <- object$model$simulate(n*nsim, exp(linpred))
#'   }
#'   
#'   dim(val) <- c(n, nsim)
#'   val <- as.data.frame(val)
#'   names(val) <- paste0("sim_", seq_len(nsim))
#'   return(val)
#' }

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
    warning("Degrees of freedom may be inacurate")
    df <- 1
  } else if (is.null(df)) {
    df <- object$df
  }
  if(dropl5 == "group") {
    l <- (A < 5)
    if (object$df == df) {df <- df - length(y) + length(y[!l]) + 1}
    y <- c(y[!l], sum(y[l]))
    A <- c(A[!l], sum(A[l]))
  } else if(dropl5 == "drop") {
    l <- (A < 5)
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
  structure(
    list(Test = vect,
         l5 = switch(dropl5,
                     drop = "dropped",
                     group = "grouped",
                     no = "preserved")),
    class = "summarysingleRmargin"
  )
}
#' vcov method for singleR class
#' @title vcov method for singleR class.
#' @param object object of clas singleRclass.
#' @param ... variables to pass to solve.
#' @description Returns a estimated covariance matrix for model coefficients
#' calculated from analytic hessian or fisher information matrix. 
#' Covariance type is taken from control parameter that have been provided
#' on call that created object.
#' 
#' @method vcov singleR
#' @return A covariance matrix for fitted coefficients obtained by inverting 
#' analitical hesian at estimated coefficients, i.e. using Cramér–Rao bound
#' with observed information/fisher information matrix.
#' @exportS3Method
vcov.singleR <- function(object, ...) {
  if(grepl(x = object$model$family, pattern = "^zot.*")) {X <- object$X[rownames(object$linear.predictors), ]} else {X <- object$X}
  switch(
    object$populationSize$control$covType,
    "observedInform" = solve(
      -object$model$makeHessian(y = object$y, X = X, 
                                weight = object$prior.weights)(object$coefficients),
      ...
    ),
    "Fisher" = solve(
      crossprod(x = as.matrix(X) * as.numeric(object$model$Wfun(prior = object$prior.weights, disp = object$dispersion, eta = object$linear.predictors)), 
                y = as.matrix(X)),
      ...
    )
  )
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
  if (grepl("zot", model$model$family) | model$model$family %in% c("chao", "zelterman")) {
    X <- model$X[rownames(model$linear.predictors), ]
  } else {
    X <- model$X
  }
  if (model$call$method == "robust") {
    W <- model$weights
  } else {
    W <- model$model$Wfun(prior = model$prior.weights, mu = model$fitt.values$mu, eta = model$linear.predictors, disp = model$dispersion)[, 1]
  }
  
  hatvector <- diag(
    tcrossprod(
      x = as.matrix(X) %*% 
        solve(crossprod(x = as.matrix(X), 
                        y = as.matrix(X * W))),
      y = as.matrix(W * X))
  )
  hatvector
}
#' dfbeta for singleRclass
#' @title TODO
#' @param model Fitted object of singleR class
#' @param method Method for finding dfbetas, either "formula" where the formula for 
#' approximate difference will be used or "simulation" where for each observation new 
#' object will be fitted, with already obtained coefficients as starting values and with 
#' just 1 iteration (or other number specified in maxit.new parameter), with i'th row 
#' removed for all observations.
#' @param maxit.new maximal number of iterations for regression
#' @param ... Arguments to pass to other methods
#' @description TODO
#' 
#' @return TODO
dfbetasingleR <- function(model, 
                          method = c("formula", "simulation"),
                          maxit.new = 1,
                          ...) {
  if (missing(method)) {method = "simulation"}
  switch (method,
          "formula" = {
            if (model$call$method == "robust") {
              W <- model$weights
            } else {
              W <- model$model$Wfun(prior = model$prior.weights[rownames(model$linear.predictors), ], mu = model$fitt.values$mu, eta = model$linear.predictors, disp = model$dispersion)
            }
            X <- model$X[rownames(model$X) %in% rownames(model$linear.predictors),]
            hatvector <- hatvalues.singleR(model, ...)
            rp <- residuals.singleR(object = model, type = "pearson")$pearson
            res <- t(
              solve(crossprod(x = X, (X * W))) %*% (t(X) * sqrt(W) * rp / sqrt(1 - hatvector))
            )},
          "simulation" = {
            X <- model$X[rownames(model$X) %in% rownames(model$linear.predictors),]
            if (model$model$family %in% c("chao", "zelterman")) {
              y <- model$y[model$y %in% c(1, 2)]
            } else {
              y <- model$y
            }
            res <- matrix(nrow = nrow(X), ncol = ncol(X) + !(is.null(model$dispersion)))
            cf <- model$coefficients
            for (k in 1:nrow(X)) {
              #cat("Iter nr.", k, "\n")
              res[k, ] <- cf - estimate_popsize.fit(
                control = control.method(
                  silent = TRUE, 
                  start = cf, 
                  maxiter = maxit.new + 1
                ),
                y = y[-k],
                X = X[-k, ],
                start = cf,
                dispersion = model$dispersion,
                family = model$model,
                prior.weights = if (length(model$prior.weights) == 1) model$prior.weights else model$prior.weights[-k],
                method = model$call$method
              )$beta
            }
          }
  )
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
  if (missing(parm)) {
    coef <- object$coefficients
    std <- object$standard_errors
  } else {
    coef <- object$coefficients[parm]
    std <- object$standard_errors[parm]
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
                                       "response",
                                       "working",
                                       "deviance",
                                       "all"),
                              ...) {
  type <- match.arg(type)
  res <- object$residuals
  disp <- object$dispersion
  wts <- object$prior.weights
  mu <- object$fitt.values
  y <- object$y
  if (object$model$family %in% c("chao", "zelterman")) {
    indx <- (y %in% 1:2)
    #mu <- mu[indx, ]
    #res <- res[indx, ]
    #if (length(wts) != 1) {wts <- wts[indx]}
    y <- y[indx]
  }
  rs <- switch(
    type,
    working = data.frame("working" = object$model$funcZ(eta = object$linear.predictors,
                                                        weight = object$weights,
                                                        y = y, mu = mu$link,
                                                        disp = object$dispersion)),
    response = res,
    pearson = data.frame("pearson" = res$mu / sqrt((1 - hatvalues(object)) * object$model$variance(mu = if (object$model$family %in% c("chao", "zelterman")) {mu$mu} else {mu$link}, disp = object$dispersion, type = "trunc"))),
    deviance = data.frame("deviance" = object$model$dev.resids(y = y, mu = mu$link, disp = disp, wt = wts)),
    all = {colnames(res) <- c("muResponse", "linkResponse");
    data.frame(
      "working" = object$model$funcZ(eta = object$linear.predictors,
                                     weight = object$weights,
                                     y = y, mu = mu$link,
                                     disp = object$dispersion),
      res,
      "pearson" = res$mu / sqrt((1 - hatvalues(object)) * object$model$variance(mu = if (object$model$family %in% c("chao", "zelterman")) {mu$mu} else {mu$link}, disp = object$dispersion, type = "trunc")),
      "deviance" = object$model$dev.resids(y = y, mu = mu$mu, disp = disp, wt = wts),
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
  cat("\n--------------------------------------------------------\n",
      "Cells with fitted frequencies of < 5 have been ", x$l5, "\n", sep = "")
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
model.matrix.singleR <- function(object, ...) {
  # TODO :
  ## - add Xvlm
  object$X
}

#' @method dfpopsize singleR
#' @exportS3Method
dfpopsize.singleR <- function(model, dfbeta = NULL, observedPop = FALSE, ...) {
  dfb <- if (is.null(dfbeta)) {dfbeta(model, ...)} else {dfbeta}
  if (model$model$family == "zelterman") {
    dfbnew <- matrix(0, ncol = ncol(dfb), nrow = model$sizeObserved)
    rownames(dfbnew) <- as.character(1:model$sizeObserved)
    dfbnew[rownames(dfb), ] <- dfb
    dfb <- dfbnew
    X <- model$X
  } else {
    X <- model$X[rownames(model$linear.predictors),]
  }
  N <- model$populationSize$pointEstimate
  res <- NULL
  range <- if(model$model$family == "zelterman") {1:model$sizeObserved} else {1:nrow(X)}
  for (k in range) {
    cf <- model$coefficients - dfb[k, ]
    disp <- model$dispersion
    if (grepl("negbin", model$model$family)) {
      disp <- cf[1]
      cf <- cf[-1]
    }
    res <- c(res, model$trcount + model$model$pointEst(disp = disp,
                                                       pw = if (length(model$prior.weights) == 1) {model$prior.weights} else {model$prior.weights[-k]},
                                                       lambda = model$model$linkinv(as.matrix(X[-k, ]) %*% cf)))
  }
  
  if(isTRUE(observedPop) & (grepl("zot", model$model$family) | model$model$family == "chao")) {
    res1 <- rep(Inf, model$sizeObserved)
    names(res1) <- 1:model$sizeObserved
    res1[rownames(model$linear.predictors)] <- N - res
    res1[is.infinite(res1)] <- -1
  } else {
    res1 <- N - res
    if (model$model$family != "zelterman") {
      names(res1) <- rownames(X)
    }
  }
  
  res1
}
#' @method summary singleR
#' @importFrom stats pt
#' @importFrom stats coef
#' @exportS3Method
summary.singleR <- function(object, test = c("t", "z"), correlation = FALSE, ...) {
  if (missing(test)) {test <- "z"}
  df.residual <- object$df.residual
  cov <- vcov.singleR(object, ...)
  pers <- residuals.singleR(object, type = "pearson", ...)
  cf <- object$coefficients
  se <- sqrt(diag(cov))
  wValues <- cf / se
  pValues <- switch (test,
                     "t" = 2 * stats::pt(q = -abs(wValues), df = df.residual),
                     "z" = 2 * stats::pnorm(q =  abs(wValues), lower.tail = FALSE)
  )
  crr <- if (isFALSE(correlation)) {NULL} else {cov / outer(se, se)}
  if(isTRUE(correlation)) {rownames(crr) <- colnames(crr) <- names(cf)}
  structure(
    list(
      call = object$call,
      coefficients = cf,
      standard_errors = se,
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
      test = test
    ),
    class = "summarysingleR"
  )
}
#' @importFrom stats cooks.distance
#' @method cooks.distance singleR
#' @exportS3Method
cooks.distance.singleR <- function(model, ...) {
  res <- ((residuals(model, type = "pearson") ** 2) * (hatvalues(model) / (length(model$coefficients))))$pearson
  names(res) <- rownames(model$linear.predictors)
  res
}
#' @method print summarysingleR
#' @exportS3Method
print.summarysingleR <- function(x, ...) {
  # ifelse is faster than if(_) {} else {}, sapply is faster than for hence the change
  signif <- sapply(x$pValues, function(k) {
    ifelse(k <= 0, "****",
           ifelse(k <= .001, "***",
                  ifelse(k <= .01, "**", 
                         ifelse(k <= .05, "*",
                                ifelse(k <= .1, ".", "")))))
  })
  ob <- data.frame(round(x$coefficients, digits = 3),
                   round(x$standard_errors, digits = 3),
                   round(x$wValues, digits = 2),
                   signif(x$pValues, digits = 2),
                   signif)
  colnames(ob) <- switch(x$test,
                         "t" = c("Estimate", "Std. Error", "t value", "P(>|t|)", ""),
                         "z" = c("Estimate", "Std. Error", "z value", "P(>|z|)", ""))
  
  print(x$call)
  cat("\nStandardised Pearson Residuals:\n")
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
  cat("AIC: ", x$aic,
      "\nBIC: ", x$bic,
      "\nDeviance: ", x$deviance,
      "\n\nLog-likelihood: ", x$logL, " on ", x$df.residual, " Degrees of freedom ",
      if (x$call$method == "robust") {
        "\nNumber of iterations: "
      } else {
        "\nNumber of calls to log-likelihood function: " # optim does not allow for accesing information
        # on number of iterations performed only a number of calls for gradient and objective function
      }, x$iter[1], 
      "\n-----------------------",
      "\nPopulation size estimation results: ",
      "\nPoint estimate ", x$populationSize$pointEstimate, 
      "\nObserved proportion: ", round(100 * x$sizeObserved / x$populationSize$pointEstimate, digits = 1), "% (N obs = ", x$sizeObserved, ")",
      if (x$call$pop.var == "bootstrap") {"\nBootstrap Std. Error "} else {"\nStd. Error "}, sqrt(x$populationSize$variance),
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
fitted.singleR <- function(object,...) {
  object$linear.predictors[, 1]
}

#' simulate
#' 
#' @param object an object representing a fitted model.
#' @param nsim number of response vectors to simulate. Defaults to \code{1}.
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param ... additional optional arguments.
#' @return a \code{data.frame} with \code{n} rows and \code{nsim} columns.
#' @seealso [stats::simulate()]
#' @examples 
#' set.seed(1)
#' N <- 10000
#' gender <- rbinom(N, 1, 0.2)
#' eta <- -1 + 0.5*gender
#' counts <- rpois(N, lambda = exp(eta))
#' df <- data.frame(gender, eta, counts)
#' df2 <- subset(df, counts > 0)
#' mod1 <-  estimate_popsize(formula = counts ~ 1 + gender, data = df2, 
#' model = "ztpoisson", method = "mle", pop.var = "analytic")
#' mod1_sims <- simulate(mod1, nsim=10)
#' colMeans(mod1_sims) 
#' @importFrom stats simulate
#' @method simulate singleR
#' @exportS3Method
simulate.singleR <- function(object, nsim=1, seed = NULL, ...) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  linpred <- fitted(object)
  n <- NROW(linpred)
  if (grepl("negbin", object$model$family)) {
    val <- object$model$simulate(n*nsim, exp(linpred), exp(-object$dispersion))
  } else {
    val <- object$model$simulate(n*nsim, exp(linpred))
  }
  
  dim(val) <- c(n, nsim)
  val <- as.data.frame(val)
  names(val) <- paste0("sim_", seq_len(nsim))
  return(val)
}

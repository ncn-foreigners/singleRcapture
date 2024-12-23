#' @title Extract population size estimation results.
#' 
#' @description An extractor function with \code{singleRStaticCountData} method for extracting
#' important information regarding pop size estimate.
#'
#' @param object object with population size estimates.
#' @param ... additional optional arguments, currently not used in \code{singleRStaticCountData} class method. 
#'
#' @return An object of class \code{popSizeEstResults} containing population size estimation results.
#' @export
popSizeEst <- function(object, ...) {
  UseMethod("popSizeEst")
}

#' @method family singleRStaticCountData
#' @importFrom stats family
#' @exportS3Method
family.singleRStaticCountData <- function(object, ...) {
  object$model
}
#' @method extractAIC singleRStaticCountData
#' @importFrom stats extractAIC
#' @exportS3Method 
extractAIC.singleRStaticCountData <- function(fit, scale, k = 2, ...) {
  -2 * fit$logL + k * length(fit$coefficients)
}
# CODE MODIFIED FROM stats:::logLik.glm
#' @method logLik singleRStaticCountData
#' @importFrom stats logLik
#' @exportS3Method 
logLik.singleRStaticCountData <- function(object, 
                                          type = c("value",
                                                   "function"), 
                                          deriv = 0:2,
                                          ...) {
  if (missing(type))
    type <- "value"
  
  if (missing(deriv))
    deriv <- 0
  
  if (type == "value") {
    val <- object$logL
    attr(val, "nobs") <- nobs(object)
    attr(val, "df") <- length(stats::coef(object))
    class(val) <- "logLik"
    val
  } else {
    object$model$makeMinusLogLike(
      y = as.numeric(model.response(model.frame(object))),
      X = model.matrix(object, type = "vlm"),
      weight = object$priorWeights,
      deriv = deriv
    )
  }
}
# CODE MODIFIED FROM stats:::model.frame.glm
#' @method model.frame singleRStaticCountData
#' @importFrom stats glm
#' @importFrom stats model.frame
#' @importFrom stats update
#' @exportS3Method 
model.frame.singleRStaticCountData <- function(formula, ...) {
  dots <- list(...)
  dotargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  
  if (length(dotargs) || is.null(formula$modelFrame)) {
    # fcall <- formula$call
    # fcall$method <- "model.frame"
    # fcall[[1L]] <- quote(stats::glm)
    # fcall[names(nargs)] <- nargs
    # env <- environment(formula$terms)
    # eval(fcall, env)
    # TODO:: low priority add na action and subset here
    combinedFromula <- singleRinternalMergeFormulas(formula$formula)
    if (!is.null(dotargs$data)) {
      jj <- all.vars(combinedFromula)[attr(terms(combinedFromula), "response")]
      if (!(jj %in% colnames(dotargs$data))) {
        combinedFromula <- update(combinedFromula, NULL ~ .)
      }
    }
    stats::model.frame(
      combinedFromula,
      data = if (is.null(dotargs$data)) eval(formula$call$data) else dotargs$data
    )
  }
  else formula$modelFrame
}

#' @method model.matrix singleRStaticCountData
#' @importFrom stats model.matrix
#' @exportS3Method 
model.matrix.singleRStaticCountData <- function(object, 
                                                type = c("lm", "vlm"), 
                                                ...) {
  if (missing(type)) type <- "lm"
  
  switch (type,
    lm = {
      X <- model.frame(object);
      X <- model.matrix(object$terms, X)
    },
    vlm = {
      X <- model.frame(object, ...);
      singleRinternalGetXvlmMatrix(
        X = X, 
        formulas = object$formula, 
        parNames = object$model$etaNames
      );
    }
  )
}

#' @method popSizeEst singleRStaticCountData
#' @rdname popSizeEst
#' @exportS3Method
popSizeEst.singleRStaticCountData <- function(object, ...) {
  object$populationSize
}

#' @method print popSizeEstResults
#' @exportS3Method 
print.popSizeEstResults <- function(x, ...) {
  cat("Point estimate: ", x$pointEstimate, 
      "\nVariance: ", x$variance, "\n", 
      (1 - x$control$alpha) * 100, "% confidence intervals:\n", 
      sep = "")
  print(x$confidenceInterval)
  
  invisible(x)
}

#' @importFrom stats fitted
#' @method fitted singleRStaticCountData
#' @exportS3Method 
fitted.singleRStaticCountData <- function(object,
                                          type = c("truncated", 
                                                   "nontruncated",
                                                   "all"),
                                          ...) {
  if (missing(type)) type <- "truncated"
  switch (type,
          truncated    = object$fittValues$truncated,
          nontruncated = object$fittValues$nontruncated,
          all          = object$fittValues
  )
}

#' @importFrom stats nobs
#' @method nobs singleRStaticCountData
#' @exportS3Method 
nobs.singleRStaticCountData <- function(object, ...) {
  object$sizeObserved
}

#' @importFrom stats df.residual
#' @method df.residual singleRStaticCountData
#' @exportS3Method 
df.residual.singleRStaticCountData <- function(object, ...) {
  object$dfResidual
}

#' @importFrom stats sigma
#' @method sigma singleRStaticCountData
#' @rdname regDiagSingleR
#' @exportS3Method 
sigma.singleRStaticCountData <- function(object, ...) {
  predict(object, type = "mean", se = TRUE)[c(3, 4)]
}

#' @importFrom stats influence
#' @method influence singleRStaticCountData
#' @rdname regDiagSingleR
#' @exportS3Method 
influence.singleRStaticCountData <- function(model, do.coef = FALSE, ...) {
  res <- list()
  hat <- hatvalues(model)
  if (NCOL(hat) > 1) {
    for (k in 1L:NCOL(hat)) {
      res[[paste0("hat:", colnames(hat)[k])]] <- hat[, k]
    }
  } else {
    res[["hat"]] <- hat[, 1]
  }
  
  if (isTRUE(do.coef)) {
    dfb <- dfbeta(model, ...)
    res[["coefficients"]] <- dfb
  }
  
  sigma <- sigma(model)
  res[["sigma:truncated"]]    <- sigma[, 1]
  res[["sigma:nontruncated"]] <- sigma[, 2]
  
  res[["dev.res"]] <- residuals(model, type = "deviance")[, 1]
  
  res[["pear.res"]] <- residuals(model, type = "pearson")[, 1]
  
  res
}

#' @importFrom stats rstudent
#' @method rstudent singleRStaticCountData
#' @rdname regDiagSingleR
#' @exportS3Method 
rstudent.singleRStaticCountData <- function(model, ...) {
  res <- residuals(model, type = "pearson")[, 1]
  hat <- hatvalues(model)[, 1]
  
  res <- res / sqrt((1 - hat))
  
  res[is.infinite(res)] <- NaN
  res
}

#' @importFrom stats rstandard
#' @method rstandard singleRStaticCountData
#' @rdname regDiagSingleR
#' @exportS3Method 
rstandard.singleRStaticCountData <- function(model,
                                             type = c("deviance", "pearson"), 
                                             ...) {
  type <- match.arg(type)
  res <- switch (type,
    pearson  = residuals(model, type = "pearsonSTD")[, 1],
    deviance = residuals(model, type = "deviance")[, 1] / 
      (sqrt(1 - hatvalues(model)[, 1]) * sigma(model)[, 1]),
  )
  
  res[is.infinite(res)] <- NaN
  res
}

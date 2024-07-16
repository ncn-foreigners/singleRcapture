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

#' @method AIC singleRStaticCountData
#' @importFrom stats AIC
#' @exportS3Method 
AIC.singleRStaticCountData <- function(object, ...) {
  2 * (length(object$coefficients) - object$logL)
}
#' @method BIC singleRStaticCountData
#' @importFrom stats BIC
#' @exportS3Method 
BIC.singleRStaticCountData <- function(object, ...) {
  length(object$coefficients) * log(nobs(object, ...)) - 2 * object$logL
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

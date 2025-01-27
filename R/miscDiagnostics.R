#' \loadmathjax

#' @title Obtain Covariance Matrix estimation
#' 
#' @description A \code{vcov} method for \code{singleRStaticCountData} class.
#' 
#' @param object object of singleRStaticCountData class.
#' @param type type of estimate for covariance matrix for now either
#' expected (Fisher) information matrix or observed information matrix.
#' @param ... additional arguments for method functions
#' 
#' @details  Returns a estimated covariance matrix for model coefficients
#' calculated from analytic hessian or Fisher information matrix usually 
#' utilizing asymptotic effectiveness of maximum likelihood estimates.
#' Covariance type is taken from control parameter that have been provided
#' on call that created \code{object} if arguments \code{type} was not specified.
#' 
#' @method vcov singleRStaticCountData
#' @return A covariance matrix for fitted coefficients, rows and columns of which 
#' correspond to parameters returned by \code{coef} method.
#' @seealso [vcovHC.singleRStaticCountData()] [sandwich::sandwich()]
#' @exportS3Method
vcov.singleRStaticCountData <- function(object, 
                         type = c("Fisher", 
                                  "observedInform"), 
                         ...) {
  if (missing(type) ){type <- object$populationSize$control$covType}
  X <- model.matrix(object, "vlm")
  res <- switch(
    type,
    "observedInform" = solve(
      -object$model$makeMinusLogLike(y = object$y, X = X,
      weight = object$priorWeights, deriv = 2, 
      offset = object$offset)(object$coefficients),
      ...
    ),
    "Fisher" = {
      if (isTRUE(object$call$method == "IRLS")) {W <- object$weights} 
        else {
          W <- object$model$Wfun(
            prior = object$priorWeights, 
            eta   = object$linearPredictors, 
            y     = as.numeric(model.response(model.frame(object)))
          )
        };
      if (isTRUE(object$control$controlMethod$checkDiagWeights)) {
        W[, (1:length(family(object)$etaNames)) ^ 2] <- ifelse(
          W[, (1:length(family(object)$etaNames)) ^ 2] < 
          object$control$controlMethod$weightsEpsilon, 
          object$control$controlMethod$weightsEpsilon, 
          W[, (1:length(family(object)$etaNames)) ^ 2]
        )
      };
      solve(
      singleRinternalMultiplyWeight(X = X, W = W) %*% X,
      ...
    )}
  )
  dimnames(res) <- list(names(object$coefficients), names(object$coefficients))
  res
}

#' @title Confidence Intervals for Model Parameters
#' 
#' @description A function that computes studentized confidence intervals
#' for model coefficients.
#' 
#' @param object object of singleRStaticCountData class.
#' @param parm names of parameters for which confidence intervals are to be 
#' computed, if missing all parameters will be considered.
#' @param level confidence level for intervals.
#' @param ... currently does nothing.
#' 
#' @method confint singleRStaticCountData
#' @return An object with named columns that include upper and 
#' lower limit of confidence intervals.
#' @exportS3Method
confint.singleRStaticCountData <- function(object,
                            parm, 
                            level = 0.95, 
                            ...) {
  std <- sqrt(diag(vcov.singleRStaticCountData(object)))
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

# functions with documented elsewhere

#' @rdname regDiagSingleR
#' @method hatvalues singleRStaticCountData
#' @importFrom stats hatvalues
#' @exportS3Method 
hatvalues.singleRStaticCountData <- function(model, ...) {
  X <- model.frame.singleRStaticCountData(model, ...)
  X <- singleRinternalGetXvlmMatrix(
    X = X, 
    formulas = model$formula, 
    parNames = model$model$etaNames
  )
  if (isTRUE(model$call$method == "IRLS")) {
    W <- model$weights
  } else {
    W <- model$model$Wfun(
      prior = model$priorWeights, 
      eta   = model$linearPredictors,
      y     = if (is.null(model$y)) model.response(model.frame(model)) 
      else model$y
    )
  }
  
  mlt <- singleRinternalMultiplyWeight(X = X, W = W)
  hatvalues <- diag(X %*% solve(mlt %*% X) %*% mlt)
  hatvalues <- matrix(
    hatvalues, 
    ncol = length(model$model$etaNames), 
    dimnames = list(
      1:(length(hatvalues) / length(model$model$etaNames)), 
      model$model$etaNames
    )
  )
  hatvalues
}

#' @method residuals singleRStaticCountData
#' @importFrom stats residuals
#' @rdname regDiagSingleR
#' @exportS3Method
residuals.singleRStaticCountData <- function(object,
                                             type = c(
                                               "pearson",
                                               "pearsonSTD",
                                               "response",
                                               "working",
                                               "deviance",
                                               "all"
                                             ),
                                             ...) {
  type <- match.arg(type)
  res <- object$residuals
  wts <- object$priorWeights
  y <- if (is.null(object$y)) stats::model.response(model.frame(object)) 
    else object$y
  
  if (type == "pearsonSTD" && length(object$model$etaNames) > 1)
    stop(paste0("Standardized pearson residuals not yet",
                "implemented for models with multiple linear predictors"))
  
  rs <- switch(
    type,
    working = as.data.frame(
      object$model$funcZ(eta = object$linearPredictors, 
                         weight = object$weights, 
                         y = y, prior = wts), 
      col.names = paste0("working:", object$model$etaNames)
    ),
    response = res,
    pearson = data.frame(
      "pearson" = res$truncated / sqrt(
        object$model$variance(eta = object$linearPredictors, type = "trunc")
      )
    ),
    pearsonSTD = data.frame(
      "pearsonSTD" = res$truncated / sqrt(
        (1 - hatvalues(object)) * 
        object$model$variance(eta = object$linearPredictors, type = "trunc")
      )
    ),
    deviance = data.frame(
      "deviance" = object$model$devResids(y = y, 
                                          eta = object$linearPredictors, 
                                          wt = wts)
    ),
    all = {colnames(res) <- c("truncatedResponse", 
                              "nontruncatedResponse");
      data.frame(
        as.data.frame(object$model$funcZ(eta = object$linearPredictors, 
                                         weight = object$weights, 
                                         y = y, prior = wts),
                      col.names = paste0("working:", object$model$etaNames)),
        res,
        "pearson" = as.numeric(res$truncated / sqrt(object$model$variance(
          eta = object$linearPredictors, type = "trunc"
        ))),
        "pearsonSTD" = if (length(object$model$etaNames) == 1) as.numeric(
          res$truncated / sqrt(
            (1 - hatvalues(object)) * 
            object$model$variance(object$linearPredictors, type = "trunc")
          )
        ) else NA,
        "deviance" = as.numeric(object$model$devResids(
          y = y, eta = object$linearPredictors, wt = wts
        )),
        row.names = rownames(object$linearPredictors)
    )}
  )
  rs
}

#' @importFrom stats cooks.distance
#' @method cooks.distance singleRStaticCountData
#' @rdname regDiagSingleR
#' @exportS3Method 
cooks.distance.singleRStaticCountData <- function(model, ...) {
  if (length(model$model$etaNames) > 1) 
    stop("Cooks distance is only implemented for single parameter families.")
  
  res <- residuals(model, type = "pearsonSTD") ^ 2
  res <- res[, 1]
  
  ht <- hatvalues(model)
  res <- (res * (ht / (length(coef(model)))))
  rownames(res) <- rownames(ht)
  res
}
#' \loadmathjax
#' @title Updating population size estimation results.
#'
#' @description A function that applies all post-hoc procedures that were taken
#' (such as heteroscedastic consistent covariance matrix estimation or bias
#' reduction) to population size estimation and standard error estimation.
#' 
#' @param object object for which update of population size estimation results will be done.
#' @param newdata optional \code{data.frame} with new data for pop size estimation.
#' @param cov an updated covariance matrix estimate.
#' @param coef optional vector of coefficients of regression on which to base 
#' population size estimation. If missing it is set to \code{coef(object)}.
#' @param weights optional vector of weights to use in population size estimation. 
#' @param control similar to \code{controlPopVar} in [estimatePopsize()].
#' If missing set to controls provided on call to \code{object}.
#' @param popVar similar to \code{popVar} in [estimatePopsize()].
#' If missing set to \code{"analytic"}.
#' @param offset offset argument for new data
#' @param weightsAsCounts for \code{singleRStaticCountData} method used to specify
#' whether weights should be treated as number of occurrences for rows in data
#' @param ... additional optional arguments, currently not used in \code{singleRStaticCountData} class method.
#'
#' @details
#' Any non specified arguments will be inferred from the \code{object}
#' 
#'
#' @return An object of class \code{popSizeEstResults} containing updated 
#' population size estimation results.
#' 
#' @examples
#' # Create simple model
#' Model <- estimatePopsize(
#'   formula = capture ~ nation + gender, 
#'   data = netherlandsimmigrant, 
#'   model = ztpoisson, 
#'   method = "IRLS"
#' )
#' # Apply heteroscedasticity consistent covariance matrix estimation
#' require(sandwich)
#' cov <- vcovHC(Model, type = "HC3")
#' summary(Model, cov = cov,
#' popSizeEst = redoPopEstimation(Model, cov = cov))
#' # Compare to results with usual covariance matrix estimation
#' summary(Model)
#' 
#' ## get confidence interval with larger significance level
#' redoPopEstimation(Model, control = controlPopVar(alpha = .000001))
#' @export
redoPopEstimation <- function(object, newdata, ...) {
  UseMethod("redoPopEstimation")
}

#' @method redoPopEstimation singleRStaticCountData
#' @rdname redoPopEstimation
#' @exportS3Method
redoPopEstimation.singleRStaticCountData <- function(object, 
                                                     newdata, 
                                                     cov, 
                                                     weights,
                                                     coef,
                                                     control,
                                                     popVar,
                                                     offset,
                                                     weightsAsCounts,
                                                     ...) {
  ### weightsAsPopCount works
  if (missing(cov)) {
    cov <- vcov
  }
  
  if (is.character(cov)) {
    cov <- get(cov, mode = "function", envir = parent.frame())
  }
  if (is.function(cov)) {
    cov <- cov(object, ...)
  }
  
  if (missing(newdata)) {
    Xvlm <- model.matrix(object, "vlm")
    
    pw <- if (missing(weights))
      object$priorWeights
    else weights
    
    offset <- if (missing(offset))
      object$offset
    else offset
    
    etaNew <- if (missing(coef))
      object$linearPredictors
    else matrix(as.matrix(Xvlm) %*% coef, ncol = length(family(object)$etaNames))
    
    if (missing(control))
      control <- object$populationSize$control
    
    if (is.null(control$bootstrapFitcontrol)) {
      control$bootstrapFitcontrol <- object$control$controlMethod
    }
    
    Y <- object$y
    X <- model.matrix(object)
    MM <- model.frame(object, ...)
    nn <- nobs(object)
    
  } else {
    if (missing(control))
      control <- object$populationSize$control
    
    if (is.null(control$bootstrapFitcontrol)) {
      control$bootstrapFitcontrol <- object$control$controlMethod
    }
    
    MM <- model.frame(object, data = newdata, ...)
    X <- model.matrix(object, data = newdata)
    Xvlm <- model.matrix(object, type = "vlm", data = newdata)
    
    Y <- model.response(MM)
    
    nn <- length(Y)
    
    pw <- if (missing(weights))
      rep(1, nn)
    else weights
    
    offset <- if (missing(offset))
      matrix(0, nrow = length(pw), ncol = length(family(object)$etaNames))
    else offset
    
    coef <- if (missing(coef)) stats::coef(object)
    
    etaNew <- matrix(
      as.matrix(Xvlm) %*% coef, 
      ncol = length(family(object)$etaNames),
      dimnames = list(
        rownames(Xvlm),
        family(object)$etaNames
      )
    ) + offset
  }
  
  singleRcaptureinternalpopulationEstimate(
    y = Y,
    formulas = object$formula,
    X = X,
    grad = object$model$makeMinusLogLike(
      y = Y,
      X = Xvlm,
      weight = pw, 
      deriv = 1
    ),
    hessian = object$model$makeMinusLogLike(
      y = Y,
      X = Xvlm, 
      weight = pw, 
      deriv = 2
    ),
    popVar = if (missing(popVar)) 
      "analytic" 
    else popVar,
    weights = pw,
    eta = etaNew,
    family = family(object),
    beta = if (missing(coef))
      stats::coef(object)
    else coef,
    control = if (missing(control))
      object$populationSize$control
    else control,
    Xvlm = Xvlm,
    W = if (isTRUE(object$call$method == "IRLS") & missing(newdata)) 
      object$weights 
    else 
      family(object)$Wfun(prior = pw, eta = etaNew, y = Y),
    sizeObserved = nn,
    modelFrame = MM,
    cov = cov,
    offset = offset,
    weightsFlag = if (missing(weightsAsCounts)) object$control$controlModel$weightsAsCounts else weightsAsCounts
  )
}
#' \loadmathjax
#' @title Predict method for singleRStaticCountData class
#' 
#' @description
#' A method for \code{predict} function, works analogous to \code{predict.glm}
#' but gives the possibility to get standard errors of 
#' mean/distribution parameters and directly get pop size estimates for new data.
#' 
#' 
#' @param object an object of \code{singleRStaticCountData} class.
#' @param newdata an optional \code{data.frame} containing new data.
#' @param type the type of prediction required, possible values are:
#' \itemize{
#'   \item \code{"response"}-- For matrix containing estimated distributions
#'   parameters.
#'   \item \code{"link"}    -- For matrix of linear predictors.
#'   \item \code{"mean"}    -- For fitted values of both \mjseqn{Y} and
#'   \mjseqn{Y|Y>0}.
#'   \item \code{"contr"}   -- For inverse probability weights (here named for 
#'   observation contribution to population size estimate).
#'   \item \code{"popSize"} -- For population size estimation. Note
#'   this results in a call to \code{redoPopEstimation} and it is
#'   usually better to call this function directly.
#' } by default set to \code{"response"}.
#' @param se.fit a logical value indicating whether standard errors should be 
#' computed. Only matters for \code{type} in \code{"response", "mean", "link"}.
#' @param na.action does nothing yet.
#' @param weights optional vector of weights for \code{type} in \code{"contr", "popSize"}.
#' @param cov optional matrix or function or character specifying either
#' a covariance matrix or a function to compute that covariance matrix.
#' By default \code{vcov.singleRStaticCountData} can be set to e.g. \code{vcovHC}.
#' @param ... arguments passed to other functions, for now this only affects
#' \code{vcov.singleRStaticCountData} method and \code{cov} function.
#' 
#' @details Standard errors are computed with assumption of regression
#' coefficients being asymptotically normally distributed, if this assumption
#' holds then each of linear predictors i.e. each row of
#' \mjseqn{\boldsymbol{\eta}=\boldsymbol{X}_{vlm}\boldsymbol{\beta}}
#' is asymptotically normally distributed and their variances are expressed by
#' well known formula. The mean \mjseqn{\mu} and distribution parameters
#' are then differentiable functions of asymptotically normally distributed
#' variables and therefore their variances can be computed using (multivariate)
#' delta method.
#'
#' @return Depending on \code{type} argument if one of \code{"response", "link", "mean"}
#' a matrix with fitted values and possibly standard errors if \code{se.fit}
#' argument was set to \code{TRUE}, if \code{type} was set to \code{"contr"}
#' a vector with inverses of probabilities, finally for \code{"popSize"}
#' an object of class \code{popSizeEstResults} with its own methods containing
#' population size estimation results.
#' 
#' @method predict singleRStaticCountData
#' @seealso [redoPopEstimation()] [stats::summary.glm()] [estimatePopsize()]
#' @exportS3Method
predict.singleRStaticCountData <- function(object,
                                           newdata,
                                           type = c("response", "link", 
                                                    "mean", "popSize", 
                                                    "contr"),
                                           se.fit = FALSE,
                                           na.action = NULL,
                                           weights,
                                           cov,
                                           ...) {
  type <- match.arg(type)
  if (missing(weights)) {
    if (missing(newdata)) {
      weights <- object$priorWeights
    } else {
      weights <- rep(1, NROW(newdata))
    }
  }
  
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
    Xvlm <- model.matrix(
      object, type = "vlm"
    )
    
    eta <- object$linearPredictors
    
    res <- switch (type,
      response = as.data.frame(
        lapply(1:length(family(object)$etaNames), FUN = function(x) {
          family(object)$links[[x]](
            eta[, x],
            inverse = TRUE
          )
      }), col.names = family(object)$etaNames),
      link     = eta,
      mean     = fitted(object, "all"),
      popSize  = popSizeEst(object, ...),
      contr    = family(object)$pointEst(
        pw     = weights,
        eta    = eta,
        contr  = TRUE,
        y      = model.response(model.frame(object))
      )
    )
  } else {
    mf <- model.frame(
      object, data = newdata
    )
    
    Xvlm <- model.matrix(
      object, type = "vlm", data = newdata
    )
    
    eta <- matrix(
      as.matrix(Xvlm) %*% stats::coef(object), 
      ncol = length(family(object)$etaNames),
      dimnames = list(
        rownames(Xvlm),
        family(object)$etaNames
      )
    )

    res <- switch (type,
      response = as.data.frame(
        lapply(1:length(family(object)$etaNames), FUN = function(x) {
        family(object)$links[[x]](
          eta[, x],
          inverse = TRUE
        )
      }), col.names = family(object)$etaNames),
      link = eta,
      mean = data.frame(
        "truncated"    = family(object)$mu.eta(eta = eta),
        "nontruncated" = family(object)$mu.eta(eta = eta, type = "nontrunc")
      ),
      popSize = redoPopEstimation(
        object = object, newdata = newdata,
        weights = if (missing(weights)) rep(1, NROW(mf)) else weights,
        cov = cov,
        ...
      ),
      contr = family(object)$pointEst(
        pw    = weights,
        eta   = eta,
        contr = TRUE,
        y     = model.response(mf) 
      )
    )
  }
  
  if (isTRUE(se.fit)) {
    cov <- vcov(object, ...)
    
    #### beta is asymptotically normal and each eta is just a linear
    #### combination of beta so it is also asymptotically normal
    #### and its variance is computed by quadratic form
    if (type != "mean") {
      se <- matrix(
        sapply(1:NROW(Xvlm), function(x) {
          (t(Xvlm[x,]) %*% cov %*% Xvlm[x,]) ^ .5
        }),
        ncol = length(family(object)$etaNames),
        dimnames = list(
          rownames(eta),
          paste0("se:", family(object)$etaNames)
        )
      )
      
      if (type == "link") 
        res <- cbind(res, se)
      
      ## since beta is asymptotically normal we use delta metod for
      ## geting standard errors of invlink(eta)
      
      if (type == "response") {
        se <- matrix(
          sapply(
            1:length(family(object)$etaNames),
            function (x) {
              se[, x, drop = FALSE] *
                abs(family(object)$links[[x]](
                  eta[, x],
                  inverse = TRUE,
                  deriv = 1
                ))
            }
          ),
          ncol = length(family(object)$etaNames),
          dimnames = dimnames(se)
        )
        res <- cbind(res,se)
      }
    } else if (type == "mean") {
      ## covariance matrix foe each row of linear predictors
      auxVec <- c(0, cumsum(rep(
        NROW(eta), 
        length.out = length(family(object)$etaNames) - 1
      )))
      
      se <- lapply(
        1:NROW(eta),
        function (x) {
          Xvlm[x + auxVec, , drop = FALSE] %*% cov %*%
            t(Xvlm[x + auxVec, , drop = FALSE])
        }
      )
      
      derivMu <- list(
        family(object)$mu.eta(eta, type = "nontrunc", deriv = 1),
        family(object)$mu.eta(eta, type = "trunc", deriv = 1)
      )
      
      res <- data.frame(
        res,
        "se:truncated" = sapply(
          1:NROW(eta),
          function (x) {
            (derivMu[[2]][x, , drop = FALSE] %*% se[[x]] %*%
               t(derivMu[[2]][x, , drop = FALSE])) ^ .5
          }
        ),
        "se:nontruncated" = sapply(
          1:NROW(eta), 
          function (x) {
            (derivMu[[1]][x, , drop = FALSE] %*% se[[x]] %*%
               t(derivMu[[1]][x, , drop = FALSE])) ^ .5
          }
        )
      )
    }
  }

  res
}
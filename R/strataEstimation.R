#' \loadmathjax

#' @title Estimate size of sub populations.
#' 
#' @description A function that estimates sizes of specific sub populations 
#' based on a capture-recapture model for the whole population.
#'
#' @param object an object on which the population size estimates should be based
#' in \code{singleRcapture} package this is a fitter \code{singleRStaticCountData} class object.
#' @param strata a specification of sub populations given by one of:
#' \itemize{
#' \item formula -- a formula to be applied to \code{model.frame} extracted from
#' the object.
#' \item Logical vector with number of entries equal to number of rows in the dataset.
#' \item A (named) list where each element is a logical vector, names of the list
#' will be used to specify names variable in returned object.
#' \item Vector of names of explanatory variables. For \code{singleRStaticCountData} method
#' for this function this specification of \code{strata} parameter will
#' result in every level of explanatory variable having its own sub population
#' for each variable specified.
#' \item If no value was provided the \code{singleRStaticCountData} method for this function 
#' will itself create sub populations based on levels of factor variables
#' in \code{model.frame}.
#' }
#' @param cov for \code{singleRStaticCountData} method an estimate of variance-covariance matrix
#' for estimate of regression parameters. It is possible to pass a function
#' such as for example \code{sandwich::vcovHC} which will be called as:
#' \code{foo(object, ...)} and a user may specify additional arguments of a 
#' function in \code{...} argument. If not provided an estimate for covariance
#' matrix will be set by calling appropriate \code{vcov} method.
#' @param alpha significance level for confidence intervals --
#' Either a single numeric value or a vector of length equal to number of 
#' sub populations specified in \code{strata}. 
#' If missing it is set to \code{.05} in \code{singleRStaticCountData} method.
#' @param ... a vector of arguments to be passed to other functions.
#' For \code{singleRStaticCountData} method for this functions arguments in \code{...} are 
#' passed to either \code{cov} if argument provided was a function or 
#' \code{vcov} if \code{cov} argument was missing at call.
#' 
#' @details In single source capture-recapture models the most frequently used
#' estimate for population size is Horvitz-Thompson type estimate:
#' 
#' \mjsdeqn{\hat{N} = \sum_{k=1}^{N}\frac{I_{k}}{\mathbb{P}(Y_{k}>0)} = 
#' \sum_{k=1}^{N_{obs}}\frac{1}{1-\mathbb{P}(Y_{k}=0)}}
#'
#' where \mjseqn{I_{k}=I_{Y_{k} > 0}} are 
#' indicator variables, with value 1 if kth unit was observed at least once 
#' and 0 otherwise and the inverse probabilistic weights weights for 
#' units observed in the data \mjseqn{\tfrac{1}{\mathbb{P}(Y_{k}>0)}}
#' are estimated using fitted linear predictors.
#' 
#' The estimates for different sub populations are made by changing the
#' \mjseqn{I_{k}=I_{Y_{k} > 0}} indicator variables to 
#' refer not to the population as a whole but to the sub populations that are 
#' being considered i.e. by changing values from 1 to 0 if kth unit is not a 
#' member of sub population that is being considered at the moment.
#' 
#' The estimation of variance for these estimates and estimation of variance for
#' estimate of population size for the whole population follow the same relation
#' as the one described above.
#' 
#' @seealso [vcov.singleRStaticCountData()] [estimatePopsize()]
#'
#' @return A \code{data.frame} object with row names being the names of specified 
#' sub populations either provided or inferred.
#' @export
stratifyPopsize <- function(object, strata, alpha, ...) {
  UseMethod("stratifyPopsize")
}



#' @method stratifyPopsize singleRStaticCountData
#' @rdname stratifyPopsize
#' @importFrom stats vcov
#' @importFrom stats contrasts
#' @exportS3Method
stratifyPopsize.singleRStaticCountData <- function(object, 
                                                   strata,
                                                   alpha, 
                                                   cov = NULL,
                                                   ...) {
  if (missing(strata)) {
    strata <- names(which(attr(object$terms, "dataClasses") == "factor"))
    strata <- strata[strata %in% attr(object$terms, "term.labels")]
    if (!length(strata)) {
      strata <- names(which(attr(object$terms, "dataClasses") == "character"))
      strata <- strata[strata %in% attr(object$terms, "term.labels")]
    }
    if (!length(strata)) {
      stop("No strata argument was provided and no factors or character columns are present in model.frame.")
    }
  }
  if (missing(alpha)) alpha <- .05
  
  # convert strata to list for all viable types of specifying the argument
  if (inherits(strata, "formula")) {
    mf <- model.frame(strata, model.frame(object))
    mmf <- model.matrix(
      strata, data = mf, 
      contrasts.arg = lapply(
        # this makes it so that all levels of factors are encoded
        subset(mf, select = sapply(mf, is.factor)),
        contrasts, contrasts = FALSE
      )
    )
    trm <- attr(mf, "terms")
    strata <- list()
    for (k in attr(trm, "term.labels")) {
      if (k %in% colnames(mf)) {
        if (is.integer(mf[, k]) | is.character(mf[, k])) {
          for (t in unique(mf[,k])) {
            strata[[paste0(k, "==", t)]] <- mf[,k] == t
          }
        } else if (is.factor(mf[, k])) {
          for (t in levels(mf[, k])) {
            strata[[paste0(k, "==", t)]] <- mf[,k] == t
          }
        }
      } else {
        tLevs <- colnames(mmf)[attr(mmf, "assign") == which(attr(trm, "term.labels") == k)]
        for (t in tLevs) {
          strata[[as.character(t)]] <- mmf[, t] == 1
        }
      }
    }
  } else if (is.list(strata)) {
    if (!all(sapply(strata, is.logical)))
      stop("Invalid way of specifying subpopulations in strata. If strata argument is a list then all elements of that list must be logical vectors.")
    
    if (length(strata[[1]]) != object$sizeObserved) 
      stop("Elements of strata object should have length equal to number of observed units.")
    
  } else if (is.logical(strata)) {
    if (length(strata) != object$sizeObserved) 
      stop("strata object should have length equal to number of observed units.")
    
    strata <- list(strata = strata)
  } else if (is.character(strata)) {
    modelFrame <- model.frame(object)
    out <- list()
    for (k in strata) {
      if (!(k %in% colnames(modelFrame))) 
        stop("Variable specified in strata is not present in model frame.")

      if (is.factor(modelFrame[, k])) {
        # this makes a difference on factors that have leves which are not represented in the data
        for (t in levels(modelFrame[, k])) {
          out[[paste0(as.character(k), "==", t)]] <- (modelFrame[, k] == t)
        }
      } else {
        for (t in unique(modelFrame[, k])) {
          out[[paste0(as.character(k), "==", t)]] <- (modelFrame[, k] == t)
        }
      }
    }
    strata <- out
  } else {
    errorMessage <- paste0(
      "Invalid way of specifying subpopulations in strata.\n", 
      "Please provide either:\n",
      "(1) - a list with logical vectors specifying different sub populations\n",
      "(2) - a single logical vector\n",
      "(3) - a formula\n",
      "(4) - a vector with names of variables by which strata will be created\n"
    )
    stop(errorMessage)
  }
  
  family <- family(object = object)
  priorWeights <- object$priorWeights
  eta <- object$linearPredictors
  Xvlm <- model.matrix(object, "vlm")
  y <- if (is.null(object$y)) model.response(model.frame(object)) else object$y
  flagWeighting <- object$control$controlModel$weightsAsCounts
  
  if (is.function(cov)) cov <- cov(object, ...)
  if (is.null(cov)) cov <- vcov(object, ...)
  
  obs <- vector(mode = "numeric", length = length(strata))
  est <- vector(mode = "numeric", length = length(strata))
  stdErr <- vector(mode = "numeric", length = length(strata))
  cnfStudent <- matrix(nrow = length(strata), ncol = 2)
  cnfChao <- matrix(nrow = length(strata), ncol = 2)
  sc <- qnorm(p = 1 - alpha / 2)
  if (length(sc) != length(strata)) sc <- rep(sc, length.out = length(strata))
  
  for (k in 1:length(strata)) {
    cond <- strata[[k]]
    
    if (isTRUE(object$control$controlModel$weightsAsCounts)) {
      obs[k] <- sum(priorWeights[cond])
    } else {
      obs[k] <- sum(cond)
    }
    
    if (obs[k] > 0) {
      est[k] <- family$pointEst(pw  = priorWeights[cond], 
                                eta = eta[cond, , drop = FALSE],
                                y   = y[cond])
      
      stdErr[k] <- family$popVar(
        pw = priorWeights[cond], 
        eta = eta[cond, , drop = FALSE], 
        cov = cov, 
        Xvlm = subset(Xvlm, subset = rep(cond, length(family$etaNames))),
        y   = y[cond]
      ) ^ .5
      
      cnfStudent[k, ] <- est[k] + c(-sc[k] * stdErr[k], sc[k] * stdErr[k])
      
      G <- exp(sc[k] * sqrt(log(1 + (stdErr[k]^2) / ((est[k] - obs[k]) ^ 2))))
      cnfChao[k, ] <- obs[k] + c((est[k] - obs[k]) / G, (est[k] - obs[k]) * G)
    } else {
      est[k] <- 0
      stdErr[k] <- 0
      cnfStudent[k, ] <- c(0, 0)
      cnfChao[k, ] <- c(0, 0)
    }
  }
  
  result <- data.frame(
    names(strata), obs, est, 
    100 * obs / est, stdErr, 
    cnfStudent[, 1], cnfStudent[, 2], 
    cnfChao[, 1], cnfChao[, 2], alpha
  )
  
  bounds <- c("LowerBound", "UpperBound")
  
  colnames(result) <- c(
    "name", "Observed", "Estimated", 
    "ObservedPercentage", "StdError", 
    paste0("normal", bounds), 
    paste0("logNormal", bounds), "confLevel"
  )
  
  result
}
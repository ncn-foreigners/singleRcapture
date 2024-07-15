#' @rdname regDiagSingleR
#' @export
dfpopsize <- function(model, ...) {
  UseMethod("dfpopsize")
}

#' @method dfpopsize singleRStaticCountData
#' @rdname regDiagSingleR
#' @exportS3Method 
dfpopsize.singleRStaticCountData <- function(model, dfbeta = NULL, ...) {
  if (isTRUE(model$call$popVar == "bootstrap")) 
    warning("dfpopsize may (in some cases) not work correctly when bootstrap was chosen as population variance estimate.")
  
  dfb <- if (is.null(dfbeta)) dfbeta(model, ...) else dfbeta
  
  X <- model.frame(model, ...)
  
  N <- model$populationSize$pointEstimate
  range <- 1:NROW(dfb)
  
  res <- vector("numeric", length = NROW(dfb))
  pw <- model$priorWeights
  y <- model.response(model.frame(model))
  
  for (k in range) {
    cf <- model$coefficients - dfb[k, ]
    
    if (isTRUE(model$control$controlModel$weightsAsCounts == FALSE)) {
      res[k] <- model$model$pointEst(
        eta = matrix(
          singleRinternalGetXvlmMatrix(
            X = X[rownames(X) != rownames(X)[k], , drop = FALSE], 
            formulas = model$formula, 
            parNames = model$model$etaNames
          ) %*% cf, 
          ncol = length(model$model$etaNames)
        ),
        y = y[-k],
        pw = pw[-k]
      )
    } else {
      # Here additional conditional is not needed since if weights are zero nothing breaks
      kk <- rep(0, length(pw))
      kk[k] <- 1
      res[k] <- model$model$pointEst(
        eta = matrix(
          model.matrix(model, "vlm") %*% cf, 
          ncol = length(model$model$etaNames)
        ),
        y = y,
        pw = pw - kk
      )
    }
  }
  
  N - res
}

#' @method dfbeta singleRStaticCountData
#' @importFrom stats dfbeta
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @rdname regDiagSingleR
#' @exportS3Method 
dfbeta.singleRStaticCountData <- function(model,
                           maxitNew = 1,
                           trace = FALSE,
                           cores = 1,
                           ...) {
  # formula method removed since it doesn't give good results will reimplement if we find better formula
  X <- model.frame.singleRStaticCountData(model, ...)
  y <- if (is.null(model$y)) stats::model.response(X) else model$y
  X <- X
  y <- y
  cf <- coef(model)
  pw <- model$priorWeights
  offset <- model$offset
  eta <- model$linearPredictors
  
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    #parallel::clusterExport(cl, c("singleRinternalGetXvlmMatrix", "cf", "y", "X", "maxitNew", "model", "pw", "offset", "eta"), envir = environment())
    
    if (isFALSE(model$control$controlModel$weightsAsCounts)) {
      res <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:NROW(X), .combine = rbind),
        ex = {
          c(cf - estimatePopsizeFit(
            control = controlMethod(
              silent = TRUE,
              maxiter = maxitNew + 1,
              ...
            ),
            y = y[-k],
            X = singleRinternalGetXvlmMatrix(
              X        = X[rownames(X) != rownames(X)[k], , drop = FALSE],
              formulas = model$formula,
              parNames = model$model$etaNames
            ),
            coefStart = cf,
            etaStart  = eta[-k, , drop = FALSE] + offset[-k, , drop = FALSE],
            family = model$model,
            priorWeights = pw[-k],
            method = if (is.null(model$call$method)) "IRLS" else model$call$method,
            offset = offset[-k, , drop = FALSE]
          )$beta)
        }
      )
    } else {
      res <- foreach::`%dopar%`(
        obj = foreach::foreach(k = 1:NROW(X), .combine = rbind),
        ex = {
          if (isFALSE(pw[k] - 1 > 0)) {
            c(cf - estimatePopsizeFit(
              control = controlMethod(
                silent = TRUE,
                maxiter = maxitNew + 1,
                ...
              ),
              y = y[-k],
              X = singleRinternalGetXvlmMatrix(
                X        = X[-k, , drop = FALSE],
                formulas = model$formula,
                parNames = model$model$etaNames
              ),
              coefStart = cf,
              etaStart  = eta[-k, , drop = FALSE] + offset[-k, , drop = FALSE],
              family = model$model,
              priorWeights = pw[-k],
              method = if (is.null(model$call$method)) "IRLS" else model$call$method,
              offset = offset[-k, , drop = FALSE]
            )$beta)
          } else {
            kk <- rep(0, length(pw))
            kk[k] <- 1
            
            c(cf - estimatePopsizeFit(
              control = controlMethod(
                silent = TRUE,
                maxiter = maxitNew + 1,
                ...
              ),
              y = y,
              X = model.matrix(model, "vlm"),
              coefStart = cf,
              etaStart  = eta + offset,
              family = model$model,
              priorWeights = pw - kk,
              method = if (is.null(model$call$method)) "IRLS" else model$call$method,
              offset = offset
            )$beta)
          }
        }
      )
    }
  } else {
    res <- matrix(nrow = nrow(X), ncol = length(cf))
    
    for (k in 1:nrow(X)) {
      if (isTRUE(trace)) {
        cat("-----\nRemoving observation number: ", k, "\n", sep = "")
      }
      if (isFALSE(model$control$controlModel$weightsAsCounts)) {
        res[k, ] <- cf - estimatePopsizeFit(
          control = controlMethod(
            silent = TRUE,
            maxiter = maxitNew + 1,
            ...
          ),
          y = y[-k],
          X = singleRinternalGetXvlmMatrix(
            X        = X[rownames(X) != rownames(X)[k], , drop = FALSE],
            formulas = model$formula, 
            parNames = model$model$etaNames
          ),
          coefStart = cf,
          etaStart  = eta[-k, , drop = FALSE] + offset[-k, , drop = FALSE],
          family = model$model,
          priorWeights = pw[-k],
          method = if (is.null(model$call$method)) "IRLS" else model$call$method,
          offset = offset[-k, , drop = FALSE]
        )$beta
      } else {
        if (isFALSE(pw[k] - 1 > 0)) {
          res[k, ] <- cf - estimatePopsizeFit(
            control = controlMethod(
              silent = TRUE,
              maxiter = maxitNew + 1,
              ...
            ),
            y = y[-k],
            X = singleRinternalGetXvlmMatrix(
              X        = X[rownames(X) != rownames(X)[k], , drop = FALSE],
              formulas = model$formula, 
              parNames = model$model$etaNames
            ),
            coefStart = cf,
            etaStart  = eta[-k, , drop = FALSE] + offset[-k, , drop = FALSE],
            family = model$model,
            priorWeights = pw[-k],
            method = if (is.null(model$call$method)) "IRLS" else model$call$method,
            offset = offset[-k, , drop = FALSE]
          )$beta
        } else {
          kk <- rep(0, length(pw))
          kk[k] <- 1
          
          res[k, ] <- cf - estimatePopsizeFit(
            control = controlMethod(
              silent = TRUE,
              maxiter = maxitNew + 1,
              ...
            ),
            y = y,
            X = model.matrix(model, "vlm"),
            coefStart = cf,
            etaStart  = eta + offset,
            family = model$model,
            priorWeights = pw - kk,
            method = if (is.null(model$call$method)) "IRLS" else model$call$method,
            offset = offset
          )$beta
        }
      }
    }
  }
  
  colnames(res) <- names(cf)
  res
}

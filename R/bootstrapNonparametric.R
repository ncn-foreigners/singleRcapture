# These functions are only used internally in the package so there is no need for documenting them
#' @importFrom graphics points
noparBoot <- function(family, formulas, y, X, modelFrame,
                      beta, weights, trcount, numboot,
                      eta, trace, visT, controlBootstrapMethod = NULL,
                      method, N, offset, weightsFlag, ...) {
  strappedStatistic <- vector("numeric", length = numboot)
  n <- length(y)
  if (isTRUE(weightsFlag)) {
    n <- sum(weights)
  }
  famName <- family$family
  
  if (length(weights) == 1) {
    weights <- rep(1, n)
  }
  
  k <- 1
  
  if (visT) {
    plot(
      1, type = "n",
      xlab = "Bootstrap sample", 
      ylab = "Value of population size estimator",
      main = expression(paste("Plot of values of ", hat(N), " obtained from bootstrap samples")),
      sub = "Points will be added in real time",
      xlim = c(0, numboot + 1), ylim = c(0, 2 * N)
    )
  }
  
  #terms <- terms(modelFrame)
  
  while (k <= numboot) {
    if (isTRUE(weightsFlag)) {
      strap        <- sample(x = 1:length(y), 
                             size = n, 
                             prob = weights / n, 
                             replace = TRUE)
      
      ystrap       <- as.numeric(y[strap])
      offsetStrap  <- offset[strap, , drop = FALSE]
      Xstrap       <- modelFrame[strap, , drop = FALSE]
      
      weightsStrap <- rep(1, length(ystrap))
      
    } else {
      strap        <- sample.int(replace = TRUE, n = n)
      
      ystrap       <- as.numeric(y[strap])
      weightsStrap <- as.numeric(weights[strap])
      offsetStrap  <- offset[strap, , drop = FALSE]
      Xstrap       <- modelFrame[strap, , drop = FALSE]
    }
    
    if (!is.data.frame(Xstrap)) {
      Xstrap <- as.data.frame(Xstrap)
      colnames(Xstrap) <- colnames(modelFrame)
    }
    
    Xstrap <- singleRinternalGetXvlmMatrix(X = Xstrap, 
                                           formulas = formulas, 
                                           family$etaNames)
    theta <- NULL
    try(
      theta <- estimatePopsizeFit(
        y = ystrap,
        X = Xstrap,
        family = family,
        control = controlBootstrapMethod,
        method = method,
        priorWeights = weightsStrap,
        coefStart = jitter(beta),
        etaStart = matrix(Xstrap %*% jitter(beta), ncol = NCOL(offsetStrap)) + offsetStrap,
        offset = offsetStrap
      )$beta,
      silent = TRUE
    )
    k <- k + 1
    
    if (is.null(theta)) {
      if (isTRUE(trace)) cat("\n")
      k <- k - 1
    } else {
      theta <- matrix(Xstrap %*% theta, ncol = length(family$etaNames)) + offsetStrap
      
      if (isTRUE(trace)) {print(summary(theta))}
      est <- family$pointEst(pw = weightsStrap, eta = theta, y = ystrap)
      
      if (visT) graphics::points(k - 1, est, pch = 1)
      if (isTRUE(trace)) cat(" Estimated population size: ", est,"\n",sep = "")
      
      strappedStatistic[k - 1] <- est
    }
  }
  
  strappedStatistic
}
# multicore
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
noparBootMultiCore <- function(family, formulas, y, X, modelFrame,
                                beta, weights, trcount, numboot,
                                eta, cores, controlBootstrapMethod = NULL,
                                method, N, offset, weightsFlag, visT, ...) {
  n <- length(y)
  if (isTRUE(weightsFlag)) {
    n <- sum(weights)
  }
  famName <- family$family
  
  if (length(weights) == 1) {
    weights <- rep(1, n)
  }
  
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  
  strappedStatistic <- foreach::`%dopar%`(
    obj = foreach::foreach(k = 1:numboot, .combine = c),
    ex = {
      theta <- NULL
      while (is.null(theta)) {
        if (isTRUE(weightsFlag)) {
          strap        <- sample(x = 1:length(y), 
                                 size = n, 
                                 prob = weights / n, 
                                 replace = TRUE)
          
          ystrap       <- as.numeric(y[strap])
          offsetStrap  <- offset[strap, , drop = FALSE]
          Xstrap       <- modelFrame[strap, , drop = FALSE]
          
          weightsStrap <- rep(1, length(ystrap))
          
        } else {
          strap        <- sample.int(replace = TRUE, n = n)
          
          ystrap       <- as.numeric(y[strap])
          weightsStrap <- as.numeric(weights[strap])
          offsetStrap  <- offset[strap, , drop = FALSE]
          Xstrap       <- modelFrame[strap, , drop = FALSE]
        }
        
        if (!is.data.frame(Xstrap)) {
          Xstrap <- as.data.frame(Xstrap)
          colnames(Xstrap) <- colnames(modelFrame)
        }
        
        Xstrap <- singleRinternalGetXvlmMatrix(X = Xstrap, 
                                               formulas = formulas, 
                                               family$etaNames)
        
        try(
          theta <- estimatePopsizeFit(
            y = ystrap,
            X = Xstrap,
            family = family,
            control = controlBootstrapMethod,
            method = method,
            priorWeights = weightsStrap,
            coefStart = jitter(beta),
            etaStart = matrix(Xstrap %*% jitter(beta), ncol = NCOL(offsetStrap)) + offsetStrap,
            offset = offsetStrap
          )$beta,
          silent = TRUE
        )
      }
      theta <- matrix(Xstrap %*% theta, ncol = length(family$etaNames)) + offsetStrap
      family$pointEst(pw = weightsStrap, eta = theta, y = ystrap)
    }
  )
  
  strappedStatistic
}
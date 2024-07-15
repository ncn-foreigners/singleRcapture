# semi parametric
#' @importFrom graphics points
semparBoot <- function(family, formulas, y, X, beta,
                       weights, trcount, numboot, eta,
                       trace, visT, controlBootstrapMethod = NULL,
                       method, N, modelFrame, offset, weightsFlag, ...) {
  strappedStatistic <- vector("numeric", length = numboot)
  n <- length(y)
  if (isTRUE(weightsFlag)) {
    n <- sum(weights)
  }
  famName <- family$family
  
  if (length(weights) == 1) {
    weights <- rep(1, n)
  }
  
  if (visT) {
    plot(
      1, type = "n",
      xlab = "Bootstrap sample", 
      ylab = "Value of population size estimator",
      main = expression(paste("Plot of values of ", hat(N), 
                              " obtained from bootstrap samples")),
      sub = "Points will be added in real time",
      xlim = c(0, numboot + 1), ylim = c(0, 2 * N)
    )
  }
  
  N <- round(sum(N))
  
  k <- 1
  while (k <= numboot) {
    # get info of how much units will be sampled
    strap <- sum(rbinom(size = 1, n = N, prob = n / N))
    # the rest are sampled uniformly form population
    if (isTRUE(weightsFlag)) {
      strap        <- sample(x = 1:length(y), 
                             size = strap, 
                             prob = weights / n, 
                             replace = TRUE)
      
      ystrap       <- as.numeric(y[strap])
      offsetStrap  <- offset[strap, , drop = FALSE]
      Xstrap       <- modelFrame[strap, , drop = FALSE]
      
      weightsStrap <- rep(1, length(ystrap))
      
    } else {
      strap <- sample.int(n = n, size = strap, replace = TRUE)
      
      ystrap       <- y[as.numeric(strap)]
      weightsStrap <- weights[as.numeric(strap)]
      #etaStrap     <- eta[as.numeric(strap), , drop = FALSE]
      offsetStrap  <- offset[as.numeric(strap), , drop = FALSE]
      Xstrap       <- modelFrame[strap, , drop = FALSE]
    }
    
    if (!is.data.frame(Xstrap)) {
      Xstrap <- as.data.frame(Xstrap)
    }

    theta <- NULL
    if (isTRUE(trace)) cat("Iteration number:", k, 
                           "sample size:", length(ystrap), sep = " ")
    colnames(Xstrap) <- colnames(modelFrame)
    
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
    
    k <- k + 1
    if (is.null(theta)) {
      if (isTRUE(trace)) cat("\n")
      k <- k - 1
    } else {
      theta <- matrix(Xstrap %*% theta, ncol = length(family$etaNames)) + offsetStrap
      est <- family$pointEst(pw = weightsStrap, eta = theta, y = ystrap)
      
      if (visT) graphics::points(k - 1, est, pch = 1)
      if (isTRUE(trace)) cat(" Estimated population size: ", est,"\n",sep = "")
      #if (visT) graphics::points(k - 1, est, pch = 1)
      
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
semparBootMultiCore <- function(family, formulas, y, X, modelFrame,
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
  
  N <- round(sum(N))
  
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  
  strappedStatistic <- foreach::`%dopar%`(
    obj = foreach::foreach(k = 1:numboot, .combine = c),
    ex = {
      theta <- NULL
      while (is.null(theta)) {
        # get info of how much units will be sampled
        strap <- sum(rbinom(size = 1, n = N, prob = n / N))
        # the rest are sampled uniformly form population
        if (isTRUE(weightsFlag)) {
          strap        <- sample(x = 1:length(y), 
                                 size = strap, 
                                 prob = weights / n, 
                                 replace = TRUE)
          
          ystrap       <- as.numeric(y[strap])
          offsetStrap  <- offset[strap, , drop = FALSE]
          Xstrap       <- modelFrame[strap, , drop = FALSE]
          
          weightsStrap <- rep(1, length(ystrap))
          
        } else {
          strap <- sample.int(n = n, size = strap, replace = TRUE)
          
          ystrap       <- y[as.numeric(strap)]
          weightsStrap <- weights[as.numeric(strap)]
          #etaStrap     <- eta[as.numeric(strap), , drop = FALSE]
          offsetStrap  <- offset[as.numeric(strap), , drop = FALSE]
          Xstrap       <- modelFrame[strap, , drop = FALSE]
        }
        
        if (!is.data.frame(Xstrap)) {
          Xstrap <- as.data.frame(Xstrap)
        }
        
        theta <- NULL
        colnames(Xstrap) <- colnames(modelFrame)
        
        
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
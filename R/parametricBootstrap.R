#' @importFrom graphics points
#' @importFrom stats rbinom
parBoot <- function(family, formulas, y, X, beta, weights,
                    trcount, numboot, eta, trace, visT,
                    controlBootstrapMethod = NULL, method,
                    modelFrame, offset, weightsFlag, ...) {
  strappedStatistic <- vector("numeric", length = numboot)
  n <- length(y)
  if (isTRUE(weightsFlag)) {
    n <- sum(weights)
  }
  if (length(weights) == 1) {
    weights <- rep(1, length(y))
  }
  
  if (family$family %in% c("chao", "zelterman")) {
    message(paste(
      "Probability model will be taken as given by poisson distribution",
      "since zelterman and chao models are based on mixture of poisson",
      "distribution semi-parametric bootstrap may be a better choice.", 
      sep = "\n"
    ))
  }
  
  
  contr <- family$pointEst(pw = weights,
                           eta = eta,
                           contr = TRUE,
                           y = y)
  
  N <- sum(contr)
  
  prob <- contr / N
  
  if (visT) {
    plot(
      1, type = "n",
      xlab = "Bootstrap sample", 
      ylab = "Value of population size estimator",
      main = expression(paste("Plot of values of ", hat(N), " obtained from bootstrap samples")),
      sub = "Points will be added in real time",
      xlim = c(0, numboot + 1), ylim = c(0, 2.5 * N)
    )
  }
  
  k <- 1
  while (k <= numboot) {
    nn <- floor(N) + stats::rbinom(n = 1, size = 1, prob = N - floor(N))
    if (isTRUE(weightsFlag)) {
      strap <- sample.int(replace = TRUE, n = length(y), size = nn, prob = prob)
      
      weightsStrap <- as.numeric(weights[strap])
      offsetStrap  <- offset[strap, , drop = FALSE]
      #etaStrap     <- eta[strap, , drop = FALSE]
      Xstrap       <- modelFrame[strap, , drop = FALSE]
      
      colnames(Xstrap) <- colnames(modelFrame)
      
      Xstrap <- singleRinternalGetXvlmMatrix(X = Xstrap, 
                                             formulas = formulas, 
                                             family$etaNames)
      
      ystrap <- family$simulate(
        n = nn,
        eta = matrix(Xstrap %*% beta, ncol = length(family$etaNames)),
        lower = -1, upper = Inf
      )
      offsetStrap  <- offsetStrap[ystrap > 0, , drop = FALSE]
      
      strap <- rep(ystrap > 0, length(family$etaNames))
      hwm <- attr(Xstrap, "hwm")
      
      Xstrap <- Xstrap[strap, , drop = FALSE]
      ystrap <- ystrap[ystrap > 0]
      
      weightsStrap <- rep(1, length(ystrap))
    } else {
      strap <- sample.int(replace = TRUE, n = length(y), size = nn, prob = prob)
      
      weightsStrap <- as.numeric(weights[strap])
      offsetStrap  <- offset[strap, , drop = FALSE]
      #etaStrap     <- eta[strap, , drop = FALSE]
      Xstrap       <- modelFrame[strap, , drop = FALSE]
      
      colnames(Xstrap) <- colnames(modelFrame)
      
      Xstrap <- singleRinternalGetXvlmMatrix(X = Xstrap, 
                                             formulas = formulas, 
                                             family$etaNames)
      
      ystrap <- family$simulate(
        n = nn,
        eta = matrix(Xstrap %*% beta, ncol = length(family$etaNames)),
        lower = -1, upper = Inf
      )
      
      weightsStrap <- weightsStrap[ystrap > 0]
      #etaStrap     <- etaStrap[ystrap > 0, , drop = FALSE]
      offsetStrap  <- offsetStrap[ystrap > 0, , drop = FALSE]
      
      strap <- rep(ystrap > 0, length(family$etaNames))
      hwm <- attr(Xstrap, "hwm")
      
      Xstrap <- Xstrap[strap, , drop = FALSE]
      ystrap <- ystrap[ystrap > 0]
    }
    
    if (isTRUE(trace)) cat("Iteration number:", k, 
                           "sample size:", length(ystrap), sep = " ")
    
    theta <- NULL
    attr(Xstrap, "hwm") <- hwm
    
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
#' @importFrom stats rbinom
parBootMultiCore <- function(family, formulas, y, X, modelFrame,
                             beta, weights, trcount, numboot,
                             eta, cores, controlBootstrapMethod = NULL,
                             method, N, offset, weightsFlag, visT, ...) {
  n <- length(y)
  if (isTRUE(weightsFlag)) {
    n <- sum(weights)
  }
  
  if (length(weights) == 1) {
    weights <- rep(1, length(y))
  }
  
  if (family$family %in% c("chao", "zelterman")) {
    message(paste("Probability model will be taken as given by poisson distribution",
                  "since zelterman and chao models are based on mixture of poisson",
                  "distribution semi-parametric bootstrap may be a better choice.", 
                  sep = "\n"))
  }
  
  
  contr <- family$pointEst(pw = weights,
                           eta = eta,
                           contr = TRUE,
                           y = y)
  N <- sum(contr)
  
  prob <- contr / N
  
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  #doRNG::registerDoRNG()
  
  ### TODO:: This gives a different results for family = "chao" and "zelterman"
  ### when compared to non paralelized version
  strappedStatistic <- foreach::`%dopar%`(
    obj = foreach::foreach(k = 1:numboot, .combine = c),
    ex = {
      theta <- NULL
      while (is.null(theta)) {
        nn <- floor(N) + stats::rbinom(n = 1, size = 1, prob = N - floor(N))
        if (isTRUE(weightsFlag)) {
          strap <- sample.int(replace = TRUE, n = length(y), size = nn, prob = prob)
          
          weightsStrap <- as.numeric(weights[strap])
          offsetStrap  <- offset[strap, , drop = FALSE]
          #etaStrap     <- eta[strap, , drop = FALSE]
          Xstrap       <- modelFrame[strap, , drop = FALSE]
          
          colnames(Xstrap) <- colnames(modelFrame)
          
          Xstrap <- singleRinternalGetXvlmMatrix(X = Xstrap, 
                                                 formulas = formulas, 
                                                 family$etaNames)
          
          ystrap <- family$simulate(
            n = nn,
            eta = matrix(Xstrap %*% beta, ncol = length(family$etaNames)),
            lower = -1, upper = Inf
          )
          offsetStrap  <- offsetStrap[ystrap > 0, , drop = FALSE]
          
          strap <- rep(ystrap > 0, length(family$etaNames))
          hwm <- attr(Xstrap, "hwm")
          
          Xstrap <- Xstrap[strap, , drop = FALSE]
          ystrap <- ystrap[ystrap > 0]
          
          weightsStrap <- rep(1, length(ystrap))
        } else {
          strap <- sample.int(replace = TRUE, n = length(y), size = nn, prob = prob)
          
          weightsStrap <- as.numeric(weights[strap])
          offsetStrap  <- offset[strap, , drop = FALSE]
          #etaStrap     <- eta[strap, , drop = FALSE]
          Xstrap       <- modelFrame[strap, , drop = FALSE]
          
          colnames(Xstrap) <- colnames(modelFrame)
          
          Xstrap <- singleRinternalGetXvlmMatrix(X = Xstrap, 
                                                 formulas = formulas, 
                                                 family$etaNames)
          
          ystrap <- family$simulate(
            n = nn,
            eta = matrix(Xstrap %*% beta, ncol = length(family$etaNames)),
            lower = -1, upper = Inf
          )
          
          weightsStrap <- weightsStrap[ystrap > 0]
          #etaStrap     <- etaStrap[ystrap > 0, , drop = FALSE]
          offsetStrap  <- offsetStrap[ystrap > 0, , drop = FALSE]
          
          strap <- rep(ystrap > 0, length(family$etaNames))
          hwm <- attr(Xstrap, "hwm")
          
          Xstrap <- Xstrap[strap, , drop = FALSE]
          ystrap <- ystrap[ystrap > 0]
        }
        
        attr(Xstrap, "hwm") <- hwm
        
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
      family$pointEst(
        pw = weightsStrap, 
        eta = matrix(Xstrap %*% theta, 
                     ncol = length(family$etaNames)) + offsetStrap, 
        y = ystrap
      )
    }
  )
  
  strappedStatistic
}

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
    # TODO:: since modelframe is needed maybe revisit it and save some memory on response
    if (FALSE) {
      # maybe add this as an option??
      # in tests this is slower than "normal" option 
      # but idk maybe if data is large enough it will be faster??
      # the problem here is that aggregate is slow and I don't want to
      # add another dependency in dplyr/data.table just to make it faster
      strap        <- sample(x = 1:length(y), 
                             size = n, 
                             prob = weights / n, 
                             replace = TRUE)
      
      ystrap       <- as.numeric(y[strap])
      offsetStrap  <- offset[strap, , drop = FALSE]
      Xstrap       <- modelFrame[strap, , drop = FALSE]
      
      # get data into right format
      XXX <- data.frame(Xstrap, offsetStrap, ystrap, singleRcaptureInternalColnameFreq = 0)
      XXX <- aggregate(singleRcaptureInternalColnameFreq ~ ., data = XXX, length)
      
      weightsStrap <- XXX$singleRcaptureInternalColnameFreq
      ystrap <- XXX[,NCOL(XXX)-1, drop = TRUE]
      offsetStrap <- XXX[,(NCOL(Xstrap)+1):(NCOL(XXX) - 2), drop = FALSE]
      offsetStrap <- as.matrix(offsetStrap)
      
      Xstrap <- XXX[, 1:NCOL(Xstrap), drop = FALSE]
      attr(Xstrap, "terms") <- terms
      # free memmory
      XXX <- NULL
    } else if (isTRUE(weightsFlag)) {
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
      #etaStrap     <- eta[as.numeric(strap), , drop = FALSE]
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
noparBootMultiCore <- function(family, formulas, y, X, modelFrame,
                                beta, weights, trcount, numboot,
                                eta, cores, controlBootstrapMethod = NULL,
                                method, N, offset, weightsFlag, ...) {
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
  
  
  ### TODO:: This gives a different results for family = "chao" and "zelterman"
  ### when compared to non paralelized version
  strappedStatistic <- foreach::`%dopar%`(
    obj = foreach::foreach(k = 1:numboot, .combine = c),
    #obj = foreach::foreach(k = 1:numboot, .export = "singleRcaptureinternalIRLSmultipar"),
    ex = {
      theta <- NULL
      while (is.null(theta)) {
        # TODO:: since modelframe is needed maybe revisit it and save some memory on response
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
                               method, N, offset, weightsFlag, ...) {
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
                             method, N, offset, weightsFlag, ...) {
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

# These functions are only used internally in the package so there is no need for documenting them
#' @importFrom graphics points
noparBoot <- function(family, formulas, y, X, modelFrame,
                      beta, weights, trcount, numboot,
                      eta, trace, visT, controlBootstrapMethod = NULL,
                      method, N, offset, ...) {
  strappedStatistic <- vector("numeric", length = numboot)
  n <- length(y)
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
  
  while (k <= numboot) {
    # TODO:: since modelframe is needed maybe revisit it and save some memory on response
    strap        <- sample.int(replace = TRUE, n = n)
    
    ystrap       <- as.numeric(y[strap])
    weightsStrap <- as.numeric(weights[strap])
    #etaStrap     <- eta[as.numeric(strap), , drop = FALSE]
    offsetStrap  <- offset[strap, , drop = FALSE]
    Xstrap       <- modelFrame[strap, , drop = FALSE]
    
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
                                method, N, offset, ...) {
  n <- length(y)
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
    ex = {
      theta <- NULL
      while (is.null(theta)) {
        # TODO:: since modelframe is needed maybe revisit it and save some memory on response
        strap        <- sample.int(replace = TRUE, n = n)
        
        ystrap       <- as.numeric(y[strap])
        weightsStrap <- as.numeric(weights[strap])
        #etaStrap     <- eta[as.numeric(strap), , drop = FALSE]
        offsetStrap  <- offset[strap, , drop = FALSE]
        Xstrap       <- modelFrame[strap, , drop = FALSE]
        
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
                       method, N, modelFrame, offset, ...) {
  strappedStatistic <- vector("numeric", length = numboot)
  n <- length(y)
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
  
  yTab <- table(y)
  yTab <- c("0" = N - sum(yTab), yTab) / N
  prob <- 0:max(as.numeric(names(yTab)))
  names(prob) <- prob
  prob[names(yTab)] <- yTab
  prob[!(names(prob) %in% names(yTab))] <- 0
  yTab <- table(y)
  getX <- function(val, num) {
    sample(which(y == names(strap)[j]), size = num, replace = TRUE)
  }

  k <- 1
  while (k <= numboot) {
    strap1 <- stats::rmultinom(n = 1, size = N, prob = prob)
    strap <- as.numeric(strap1)
    names(strap) <- rownames(strap1)
    strap <- strap[as.numeric(names(strap)) != 0]
    rows <- NULL
    for (j in 1:length(strap)) {
      rows <- c(rows, getX(val = names(strap)[j], num = strap[j]))
    }
    strap <- rows
    
    ystrap       <- y[as.numeric(strap)]
    weightsStrap <- weights[as.numeric(strap)]
    #etaStrap     <- eta[as.numeric(strap), , drop = FALSE]
    offsetStrap  <- offset[as.numeric(strap), , drop = FALSE]
    Xstrap       <- modelFrame[strap, , drop = FALSE]
    
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
                               method, N, offset, ...) {
  n <- length(y)
  famName <- family$family
  
  if (length(weights) == 1) {
    weights <- rep(1, n)
  }
  
  N <- round(sum(N))
  
  yTab <- table(y)
  yTab <- c("0" = N - sum(yTab), yTab) / N
  prob <- 0:max(as.numeric(names(yTab)))
  names(prob) <- prob
  prob[names(yTab)] <- yTab
  prob[!(names(prob) %in% names(yTab))] <- 0
  yTab <- table(y)
  getX <- function(val, num) {
    sample(which(y == names(strap)[j]), size = num, replace = TRUE)
  }
  
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  
  ### TODO:: This gives a different results for family = "chao" and "zelterman"
  ### when compared to non paralelized version
  strappedStatistic <- foreach::`%dopar%`(
    obj = foreach::foreach(k = 1:numboot, .combine = c),
    ex = {
      theta <- NULL
      while (is.null(theta)) {
        strap1 <- stats::rmultinom(n = 1, size = N, prob = prob)
        strap <- as.numeric(strap1)
        names(strap) <- rownames(strap1)
        strap <- strap[as.numeric(names(strap)) != 0]
        rows <- NULL
        for (j in 1:length(strap)) {
          rows <- c(rows, getX(val = names(strap)[j], num = strap[j]))
        }
        strap <- rows
        
        ystrap       <- y[as.numeric(strap)]
        weightsStrap <- weights[as.numeric(strap)]
        #etaStrap     <- eta[as.numeric(strap), , drop = FALSE]
        offsetStrap  <- offset[as.numeric(strap), , drop = FALSE]
        Xstrap       <- modelFrame[strap, , drop = FALSE]
        
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
parBoot <- function(family,
                    formulas,
                    y, X,
                    beta,
                    weights,
                    trcount,
                    numboot,
                    eta,
                    trace,
                    visT,
                    controlBootstrapMethod = NULL,
                    method,
                    modelFrame,
                    offset,
                    ...) {
  strappedStatistic <- vector("numeric", length = numboot)
  n <- length(y)
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
    strap <- sample.int(replace = TRUE, n = n, size = nn, prob = prob)
    
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
                             method, N, offset, ...) {
  n <- length(y)
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
        strap <- sample.int(replace = TRUE, n = n, size = nn, prob = prob)
        
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

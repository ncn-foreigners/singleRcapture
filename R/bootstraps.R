# These functions are only used internally in the package so there is no need for documenting them
#' @importFrom graphics points
noparBoot <- function(family, formulas, y, X, modelFrame,
                      beta, weights, trcount, numboot,
                      eta, trace, visT, controlBootstrapMethod = NULL,
                      method, N, ...) {
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
    # TODO:: since modelframe is needed maybe revisit it and save some memmory on response
    strap <- sample.int(replace = TRUE, n = n)
    ystrap <- as.numeric(y[strap])
    weightsStrap <- as.numeric(weights[strap])
    Xstrap <- modelFrame[strap, , drop = FALSE]
    if (!is.data.frame(Xstrap)) {
      Xstrap <- as.data.frame(Xstrap)
      colnames(Xstrap) <- colnames(modelFrame)
    }
    
    wch <- singleRcaptureinternalDataCleanupSpecialCases(
    family = family, observed = ystrap, popVar = "analytic")
    if (famName == "zelterman") {
      Xstrap1 <- singleRinternalGetXvlmMatrix(
      X = Xstrap, formulas = formulas, family$etaNames)
    }
    Xstrap <- singleRinternalGetXvlmMatrix(
    X = Xstrap, formulas = formulas, family$etaNames)
    theta <- NULL
    try(
      theta <- estimatePopsize.fit(
        y = ystrap[wch$reg],
        X = Xstrap,
        family = family,
        control = controlBootstrapMethod,
        method = method,
        priorWeights = weightsStrap[wch$reg],
        start = jitter(beta)
      )$beta,
      silent = TRUE
    )
    k <- k + 1
    
    if (is.null(theta)) {
      if (isTRUE(trace)) cat("\n")
      k <- k - 1
    } else {
      if (famName == "zelterman") {
        theta <- matrix(Xstrap1 %*% theta, ncol = length(family$etaNames))
      } else {
        theta <- matrix(Xstrap %*% theta, ncol = length(family$etaNames))
      }
      if (isTRUE(trace)) {print(summary(theta))}
      est <- family$pointEst(pw = weightsStrap[wch$est], eta = theta) + wch$trr
      if (visT) graphics::points(k - 1, est, pch = 1)
      if (isTRUE(trace)) cat(" Estimated population size: ", est,"\n",sep = "")
      #if (visT) graphics::points(k - 1, est, pch = 1)
      
      strappedStatistic[k - 1] <- est
    }
  }
  
  strappedStatistic
}
# semi parametric
semparBoot <- function(family, formulas, y, X, beta,
                       weights, trcount, numboot, eta,
                       trace, visT, controlBootstrapMethod = NULL,
                       method, N, modelFrame, ...) {
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
    
    ystrap <- y[as.numeric(strap)]
    weightsStrap <- weights[as.numeric(strap)]
    Xstrap <- modelFrame[strap, , drop = FALSE]
    if (!is.data.frame(Xstrap)) {
      Xstrap <- as.data.frame(Xstrap)
    }

    wch <- singleRcaptureinternalDataCleanupSpecialCases(
    family = family, observed = ystrap, popVar = "analytic")
    theta <- NULL
    if (isTRUE(trace)) cat("Iteration number:", k, 
                           "sample size:", length(ystrap), sep = " ")
    colnames(Xstrap) <- colnames(modelFrame)
    if (famName == "zelterman") {
      Xstrap1 <- singleRinternalGetXvlmMatrix(
      X = Xstrap, formulas = formulas, family$etaNames)
    }
    Xstrap <- singleRinternalGetXvlmMatrix(
    X = Xstrap, formulas = formulas, family$etaNames)
    try(
      theta <- estimatePopsize.fit(
        y = ystrap[wch$reg],
        X = Xstrap,
        family = family,
        control = controlBootstrapMethod,
        method = method,
        priorWeights = weightsStrap[wch$reg],
        start = jitter(beta)
      )$beta,
      silent = TRUE
    )
    
    k <- k + 1
    if (is.null(theta)) {
      if (isTRUE(trace)) cat("\n")
      k <- k - 1
    } else {
      if (famName == "zelterman") {
        theta <- matrix(Xstrap1 %*% theta, ncol = length(family$etaNames))
      } else {
        theta <- matrix(Xstrap %*% theta, ncol = length(family$etaNames))
      }
      est <- family$pointEst(pw = weightsStrap[wch$est], eta = theta) + wch$trr
      if (visT) graphics::points(k - 1, est, pch = 1)
      if (isTRUE(trace)) cat(" Estimated population size: ", est,"\n",sep = "")
      #if (visT) graphics::points(k - 1, est, pch = 1)
      
      strappedStatistic[k - 1] <- est
    }
  }
  
  strappedStatistic
}
#' @importFrom stats rpois
#' @importFrom stats rnbinom
#' @importFrom stats rgeom
#' @importFrom stats optim
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
                    ...) {
  strappedStatistic <- vector("numeric", length = numboot)
  n <- length(y)
  famName <- family$family
  if (length(weights) == 1) {
    weights <- rep(1, length(y))
  }
  
  if (family$family %in% c("chao", "zelterman")) {
    cat("Probability model will be taken as given by poisson distribution",
        "since zelterman and chao models are based on mixture of poisson",
        "distribution semi-parametric bootstrap may be a better choice.", 
        sep = "\n")
  }
  
  dataFunc <- family$simulate
  
  
  if (isFALSE(grepl(x = famName, pattern = "^(zot|cha).*"))) {
    contr <- family$pointEst(pw = weights,
                             eta = eta,
                             contr = TRUE)
  } else {
    contr <- vector(mode = "numeric", length = length(y))
    cond <- switch(
      substr(famName, 1, 3),
      "zot" = (y != 1),
      "cha" = (y < 3)
    )
    contr[cond] <- family$pointEst(pw = weights[cond],
                                   eta = eta,
                                   contr = TRUE)
    contr[!cond] <- 1
  }
  N <- round(sum(contr))
  prob <- contr - floor(contr)
  contr <- floor(contr) + stats::rbinom(n = n, size = 1, prob = prob)
  
  prob <- contr / N
  
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
  
  k <- 1
  while (k <= numboot) {
    strap <- sample.int(replace = TRUE, n = n, size = N, prob = prob)
    weightsStrap <- as.numeric(weights[strap])
    Xstrap <- modelFrame[strap, , drop = FALSE]
    colnames(Xstrap) <- colnames(modelFrame)
    
    Xstrap <- singleRinternalGetXvlmMatrix(
    X = Xstrap, formulas = formulas, family$etaNames)

    ystrap <- dataFunc(
      n = N,
      eta = matrix(Xstrap %*% beta, ncol = length(family$etaNames)),
      lower = -1, upper = Inf
    )
    weightsStrap <- weightsStrap[ystrap > 0]
    strap <- rep(FALSE, length(family$etaNames) * length(ystrap))
    strap[rep(ystrap > 0, length(family$etaNames))] <- TRUE
    hwm <- attr(Xstrap, "hwm")
    Xstrap <- subset(Xstrap, subset = strap)
    ystrap <- ystrap[ystrap > 0]

    if (isTRUE(trace)) cat("Iteration number:", k, 
                           "sample size:", length(ystrap), sep = " ")
    wch <- singleRcaptureinternalDataCleanupSpecialCases(
    family = family, observed = ystrap, popVar = "analytic")
    
    theta <- NULL
    if (famName == "zelterman") {
      Xstrap1 <- Xstrap
    }
    Xstrap <- subset(Xstrap, subset = rep(wch$reg, length(family$etaNames)))
    attr(Xstrap, "hwm") <- hwm
    
    try(
      theta <- estimatePopsize.fit(
        y = ystrap[wch$reg],
        X = Xstrap,
        family = family,
        control = controlBootstrapMethod,
        method = method,
        priorWeights = weightsStrap[wch$reg],
        start = jitter(beta)
      )$beta,
      silent = TRUE
    )
    
    k <- k + 1
    
    if (is.null(theta)) {
      if (isTRUE(trace)) cat("\n")
      k <- k - 1
    } else {
      if (famName != "zelterman") {
        theta <- matrix(Xstrap %*% theta, ncol = length(family$etaNames))
      } else {
        theta <- matrix(Xstrap1 %*% theta, ncol = length(family$etaNames))
      }
      est <- family$pointEst(pw = weightsStrap[wch$est], eta = theta) + wch$trr
      if (visT) graphics::points(k - 1, est, pch = 1)
      if (isTRUE(trace)) cat(" Estimated population size: ", est,"\n",sep = "")
      #if (visT) graphics::points(k - 1, est, pch = 1)
      
      strappedStatistic[k - 1] <- est
    }
    
  }
  
  strappedStatistic
}
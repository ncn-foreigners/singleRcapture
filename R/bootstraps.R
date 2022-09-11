# These functions are only used internally in the package so there is no need for documenting them
noparBoot <- function(family,
                      formulas,
                      y, X,
                      beta,
                      weights,
                      trcount,
                      numboot,
                      eta,
                      trace,
                      visT,
                      hwm,
                      control.bootstrap.method = NULL,
                      method) {
  strappedStatistic <- NULL
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
      main = "Points will be added in real time.",
      xlim = c(0, numboot + 1), ylim = c(0, 5 * n)
    )
  }
  
  while (k <= numboot) {
    strap <- sample.int(replace = TRUE, n = n)
    ystrap <- as.numeric(y[strap])
    weightsStrap <- as.numeric(weights[strap])
    strapXvlm <- NULL
    for (j in 1:family$parNum) {
      strapXvlm <- c(strapXvlm, strap + (j - 1) * n)
    }
    Xstrap <- as.matrix(X[strapXvlm, ]) # Xvlm matrix
    
    wch <- singleRcaptureinternalDataCleanupSpecialCases(family = family, observed = ystrap, pop.var = "analytic")
    theta <- NULL
    try(
      theta <- estimate_popsize.fit(
        y = ystrap[wch$reg],
        X = Xstrap[rep(wch$reg, family$parNum), ],
        family = family,
        control = control.bootstrap.method,
        method = method,
        prior.weights = weightsStrap[wch$reg],
        start = jitter(beta),
        hwm = hwm
      )$beta,
      silent = TRUE
    )
    k <- k + 1
    
    if (is.null(theta)) {
      k <- k - 1
    } else {
      theta <- matrix(Xstrap %*% theta, ncol = family$parNum)
      if (isTRUE(trace)) {print(summary(theta))}
      est <- family$pointEst(pw = weightsStrap[wch$est], eta = theta[wch$est, ]) + wch$trr
      if (visT) points(k - 1, est, col = ((k - 1) %% 4) + 1, pch = 20)
      
      strappedStatistic <- c(strappedStatistic,est)
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
                    hwm,
                    control.bootstrap.method = NULL,
                    method) {
  ###################
  ### TODO update ###
  ###################
  stop("Parametric bootstrap not yet updated")
  strappedStatistic <- NULL
  n <- length(y)
  famName <- family$family
  if (length(weights[[1]]) == 1) {
    weights[[1]] <- rep(1, length(lambda))
  }
  
  if (length(weights[[2]]) == 1) {
    weights[[2]] <- rep(1, n)
  }
  
  if (family$family %in% c("chao", "zelterman")) {
    cat("Probability model will be taken as given by poisson distribution",
    "since zelterman and chao models are based on mixture of poisson",
    "distribution semi-parametric bootstrap may be a better choice.", 
    sep = "\n")
  } else if (family$family %in% c("ztnegbin",
                                  "zotnegbin")) {
    cat("Due to many possible computational problems with truncated binomial models",
    "bootstrap samples will be drawn until specified number of them can be fitted", 
    "this may significantly contribute to increase in runtime.", 
    sep = "\n")
  }
  dataFunc <- switch(famName,
  "ztpoisson"  = function(lambda, n, disp) {stats::rpois(n = n, lambda = lambda)},
  "chao"       = function(lambda, n, disp) {stats::rpois(n = n, lambda = lambda)},
  "zelterman"  = function(lambda, n, disp) {stats::rpois(n = n, lambda = lambda)},
  "zotpoisson" = function(lambda, n, disp) {stats::rpois(n = n, lambda = lambda)},
  "ztnegbin"   = function(lambda, n, disp) {stats::rnbinom(n = n, mu = lambda, size = exp(-disp))},
  "zotnegbin"  = function(lambda, n, disp) {stats::rnbinom(n = n, mu = lambda, size = exp(-disp))},
  "ztgeom"     = function(lambda, n, disp) {stats::rgeom(n = n, prob = (1 / (1 + lambda)))},
  "zotgeom"    = function(lambda, n, disp) {stats::rgeom(n = n, prob = (1 / (1 + lambda)))})

  if (isFALSE(grepl(x = famName, pattern = "^(zot|cha).*"))) {
  contr <- family$pointEst(disp = dispersion, 
                           pw = weights[[2]],
                           lambda = lambda,
                           contr = TRUE)
  } else {
    contr <- vector(mode = "numeric", length = length(y))
    cond <- switch(
      substr(famName, 1, 3),
      "zot" = (y != 1),
      "cha" = (y < 3)
    )
    contr[cond] <- family$pointEst(disp = dispersion, 
                                   pw = weights[[1]], 
                                   lambda = lambda,
                                   contr = TRUE)
    contr[!cond] <- 1
  }
  N <- round(sum(contr))
  prob <- contr - floor(contr)
  contr <- floor(contr) + stats::rbinom(n = n, size = 1, prob = prob)
  
  prob <- contr / N
  weights <- weights[[2]]
  getlambda <- family$linkinv

  k <- 1
  while (k <= numboot) {
    strap <- sample.int(replace = TRUE, n = n, size = N, prob = prob)
    weightsstrap <- as.numeric(weights[strap])
    Xstrap <- as.matrix(X[strap, ])

    #print(summary(getlambda(Xstrap %*% beta)))
    ystrap <- dataFunc(n = N,
                       lambda = getlambda(Xstrap %*% beta),
                       disp = dispersion)
    weightsstraptemp <- weightsstrap
    Xstraptemp <- Xstrap

    df <- data.frame(ystrap,
                     weightsstraptemp,
                     Xstraptemp)
    trcount <- 0
    if (grepl(x = famName, pattern = "^zot.*")) trcount <- nrow(df[df["ystrap"] == 1, ])
    ifelse(grepl(x = famName, pattern = "^zot.*"),
           df <- df[df["ystrap"] > 1, ],
           df <- df[df["ystrap"] > 0, ])

    if (isTRUE(trace)) cat("Iteration number:", k, "sample size:", nrow(df), "\n", sep = " ")
    ystrap <- as.numeric(unlist(df[, 1]))
    weightsstrap <- weightsstraptemp <- as.numeric(unlist(df[, 2]))
    Xstrap <- Xstraptemp <- as.matrix(df[, -c(1, 2)])
    
    if (famName == "zelterman") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      df <- df[df[1] == 1 | df[1] == 2, ]
      
      ystrap <- as.numeric(df[, 1])
      weightsstraptemp <- as.numeric(unlist(df[, 2]))
      Xstraptemp <- as.matrix(df[, -c(1, 2)])
    } else if (famName == "chao") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      trcount <- nrow(df[df[1] > 2, ])

      df <- df[df[1] == 1 | df[1] == 2, ]
      
      ystrap <- as.numeric(df[, 1])
      weightsstraptemp <- as.numeric(unlist(df[, 2]))
      weightsstrap <- weightsstraptemp
      Xstraptemp <- as.matrix(df[, -c(1, 2)])
      Xstrap <- Xstraptemp
    }
    theta <- NULL
    
    try(
      {theta <- estimate_popsize.fit(
        y = ystrap,
        X = Xstraptemp,
        family = family,
        control = control.bootstrap.method,
        method = method,
        prior.weights = weightsstraptemp,
        start = c(dispersion, beta),
        dispersion = dispersion,
        omegaTheta = NULL
      )$beta;
      if (grepl(x = family$family, pattern = "negbin")) {theta <- theta[-1]}},
      silent = TRUE
    )
    
    k <- k + 1
    
    if (is.null(theta)) {
      k <- k - 1
    } else {
      theta <- getlambda(Xstrap %*% theta)

      strappedStatistic <- c(strappedStatistic,
                             family$pointEst(disp = dispersion,
                                             pw = weightsstrap,
                                             lambda = theta) + trcount)
    }

  }
  
  strappedStatistic
}
semparBoot <- function(family,
                       formulas,
                       y, X,
                       beta,
                       weights,
                       trcount,
                       numboot,
                       eta,
                       trace,
                       visT,
                       hwm,
                       control.bootstrap.method = NULL,
                       method) {
  ###################
  ### TODO update ###
  ###################
  stop("Not yet updated")
  strappedStatistic <- NULL
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
      main = "Points will be added in real time.",
      xlim = c(0, numboot + 1), ylim = c(0, 5 * n)
    )
  }
  
  if (isFALSE(grepl(x = famName, pattern = "^(zot|cha).*"))) {
    N <- round(family$pointEst(pw = weights[[2]], 
                               eta = eta) + trcount)
  } else {
    N <- round(family$pointEst(pw = weights[[1]], 
                               eta = eta) + trcount)
  }
  
  
  yTab <- table(y)
  yTab <- c("0" = N - sum(yTab), yTab) / N
  prob <- 0:max(as.numeric(names(yTab)))
  names(prob) <- prob
  prob[names(yTab)] <- yTab
  prob[!(names(prob) %in% names(yTab))] <- 0
  weights <- weights[[2]]
  getlambda <- family$linkinv
  yTab <- table(y)
  dfPerm <- data.frame(y, weights, X)
  getX <- function(val, num) {
    sample(rownames(dfPerm[dfPerm$y == val,]), size = num, replace = TRUE)
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
    
    ystrap <- as.numeric(dfPerm[strap, 1])
    weightsstrap <- as.numeric(dfPerm[strap, 2])
    Xstrap <- as.matrix(dfPerm[strap, -c(1, 2)])
    weightsstraptemp <- weightsstrap
    Xstraptemp <- Xstrap
    
    df <- data.frame(ystrap,
                     weightsstraptemp,
                     Xstraptemp)
    trcount <- 0
    if (grepl(x = famName, pattern = "^zot.*")) trcount <- nrow(df[df["ystrap"] == 1, ])
    ifelse(grepl(x = famName, pattern = "^zot.*"),
           df <- df[df["ystrap"] > 1, ],
           df <- df[df["ystrap"] > 0, ])
    
    if (isTRUE(trace)) cat("Iteration number:", k, "sample size:", nrow(df), "\n", sep = " ")
    ystrap <- as.numeric(unlist(df[, 1]))
    weightsstrap <- weightsstraptemp <- as.numeric(unlist(df[, 2]))
    Xstrap <- Xstraptemp <- as.matrix(df[, -c(1, 2)])
    
    if (famName == "zelterman") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      df <- df[df[1] == 1 | df[1] == 2, ]
      
      ystrap <- as.numeric(df[, 1])
      weightsstraptemp <- as.numeric(unlist(df[, 2]))
      Xstraptemp <- as.matrix(df[, -c(1, 2)])
    } else if (famName == "chao") {
      df <- data.frame(ystrap,
                       weightsstrap,
                       Xstrap)
      
      trcount <- nrow(df[df[1] > 2, ])
      
      df <- df[df[1] == 1 | df[1] == 2, ]
      
      ystrap <- as.numeric(df[, 1])
      weightsstraptemp <- as.numeric(unlist(df[, 2]))
      weightsstrap <- weightsstraptemp
      Xstraptemp <- as.matrix(df[, -c(1, 2)])
      Xstrap <- Xstraptemp
    }
    theta <- NULL
    
    try(
      {theta <- estimate_popsize.fit(
        y = ystrap,
        X = Xstraptemp,
        family = family,
        control = control.bootstrap.method,
        method = method,
        prior.weights = weightsstraptemp,
        start = c(dispersion, beta),
        dispersion = dispersion,
        omegaTheta = NULL
      )$beta;
      if (grepl(x = family$family, pattern = "negbin")) {theta <- theta[-1]}},
      silent = TRUE
    )
    
    k <- k + 1
    
    if (is.null(theta)) {
      k <- k - 1
    } else {
      theta <- getlambda(Xstrap %*% theta)
      
      strappedStatistic <- c(strappedStatistic,
                             family$pointEst(disp = dispersion,
                                             pw = weightsstrap,
                                             lambda = theta) + trcount)
    }
  }
  
  strappedStatistic
}
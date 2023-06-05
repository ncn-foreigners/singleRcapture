#' @rdname singleRmodels
#' @importFrom stats uniroot
#' @importFrom stats dnbinom
#' @importFrom rootSolve multiroot
#' @export
ztnegbin <- function(nSim = 1000, epsSim = 1e-8, eimStep = 6,
                     lambdaLink = c("log", "neglog"), 
                     alphaLink = c("log", "neglog"),
                     ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  if (missing(alphaLink))  alphaLink <- "log"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink, alphaLink)
  
  lambdaLink <- switch(lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  alphaLink <- switch(alphaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  links[1:2] <- c(lambdaLink, alphaLink)
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    if (!deriv) {
      switch (type,
          "nontrunc" = lambda,
          "trunc" = lambda / (1 - (1 + alpha * lambda) ^ (-1 / alpha))
      )
    } else {
      switch (type,
          "nontrunc" = {
            matrix(c(1, 0) * c(
              lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
              alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
            ), ncol = 2)
          },
          "trunc" = {
            matrix(c(
              (alpha * lambda + 1) ^ (1 / alpha - 1) *
              ((alpha * lambda + 1) ^ (1 / alpha + 1) +
              (-alpha - 1) * lambda - 1) / 
              ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2, 
              lambda * (lambda * alpha + 1) ^ (1 / alpha - 1) *
              ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha) /
              (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)
            ) * c(
              lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
              alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
            ), ncol = 2)
          }
      )
    }
  }

  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    switch (type,
    nontrunc = lambda * (1 + alpha * lambda),
    trunc = (lambda + alpha * (lambda ^ 2) - alpha * (lambda ^ 2) * P0) / ((1 - P0) ^ 2)
    )
  }
  
  compdigamma <- function(y, alpha) {
    #temp <- 0:(y-1)
    #sum(-(alpha ^ (-2)) / (temp + 1 / alpha))
    (-digamma(y + 1 / alpha) + digamma(1 / alpha)) / (alpha ^ 2)
  }
  
  
  comptrigamma <- function(y, alpha) {
    #temp <- 0:(y-1)
    #sum(-(temp ^ 2) / (((1 + temp * alpha)) ^ 2))
    (2 * (digamma(y + 1 / alpha) - digamma(1 / alpha)) * alpha +
     trigamma(y + 1 / alpha) - trigamma(1 / alpha)) / (alpha ^ 4)
  }
  
  # Computing the expected value of di/trigamma functions on (y + 1/alpha)
  
  compExpect <- function(eta) {
    lambda <- lambdaLink(eta[1], inverse = TRUE)
    alpha  <-  alphaLink(eta[2], inverse = TRUE)
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    res <- res1 <- 0
    k <- 1
    finished <- c(FALSE, FALSE)
    while ((k < nSim) & !all(finished)) {
      prob <- stats::dnbinom(x = k:(k + eimStep), size = 1 / alpha, mu = lambda) / (1 - P0)
      if (any(!is.finite(prob))) prob <- 0
      toAdd <- cbind(compdigamma(y = k:(k + eimStep), alpha = alpha), 
                     comptrigamma(y = k:(k + eimStep), alpha = alpha)) * prob
      toAdd <- colSums(toAdd)
      k <- k + eimStep + 1
      res <- res + toAdd
      finished <- abs(toAdd) < epsSim
    }
    res
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    Ey <- mu.eta(eta = eta)
    Edig <- apply(X = eta, MARGIN = 1, FUN = function(x) {compExpect(x)})
    Etrig <- Edig[2,]
    Edig <- Edig[1,]
    
    G00 <- ((Etrig + ((lambda * alpha + 1) ^ (1 / alpha) * 
    (lambda ^ 2 * alpha ^ 2 + 2 * lambda * alpha + 1) * log(lambda * alpha + 1) ^ 2 +
    ((lambda * alpha + 1) ^ (1 / alpha) * (2 * lambda ^ 2 * alpha ^ 3 + 
    (4 * lambda - 2 * lambda ^ 2) * alpha ^ 2 + (2 - 2 * lambda) * alpha) + 
    (lambda * alpha + 1) ^ (2 / alpha) * (-2 * lambda ^ 2 * alpha ^ 3 - 4 * lambda * alpha ^ 2 - 2 * alpha)) * 
    log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (2 / alpha) * 
    ((3 * lambda ^ 2 - 2 * Ey * lambda) * alpha ^ 3 +
    (2 * lambda - Ey) * alpha ^ 2) + (lambda * alpha + 1) ^ (1 / alpha) *
    ((4 * Ey * lambda - 3 * lambda ^ 2) * alpha ^ 3 +
    (lambda ^ 2 - 2 * lambda + 2 * Ey) * alpha ^ 2) -
    2 * Ey * lambda * alpha ^ 3 - Ey * alpha ^ 2) /
    (alpha ^ 4 * (lambda * alpha + 1) ^ 2 * 
    ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)) * 
    alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
    (Edig + ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) +
    (Ey - lambda) * alpha * (lambda * alpha + 1) ^ (1 / alpha) - Ey * alpha) /
    (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) - 1))) *
    alphaLink(eta[, 1], inverse = TRUE, deriv = 2))
    
    # mixed derivative
    G01 <- -((alpha * lambda + 1) ^ (1 / alpha + 1) *
    log(alpha * lambda + 1) + (alpha * lambda + 1) ^ (1 / alpha) *
    ((alpha ^ 2 - alpha) * lambda - 2 * alpha ^ 2 * Ey) +
    (alpha * lambda + 1) ^ (2 / alpha) * (alpha ^ 2 * Ey - alpha ^ 2 * lambda) + alpha ^ 2 * Ey) /
    (alpha ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)
    
    G01 <- G01 * prior * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
           lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    # second beta derivative
    
    G11 <- (((alpha * lambda + 1) ^ (2 / alpha) * (alpha * lambda ^ 2 - 2 * alpha * Ey * lambda - Ey) +
    (alpha * lambda + 1) ^ (1 / alpha) * ((1 - alpha) * lambda ^ 2 + 4 * alpha * Ey * lambda + 2 * Ey) - 2 * alpha * Ey * lambda - Ey) /
    (lambda ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 -
    ((lambda - Ey) * (alpha * lambda + 1) ^ (1 / alpha) + Ey) /
    (lambda * (alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) - 1)) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)) * prior
    
    matrix(
      -c(G11, # lambda
         G01, # mixed
         G01, # mixed
         G00  # alpha
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "alpha")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    dig <- compdigamma(y = y, alpha = alpha)

    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    # log(alpha) derivative
    dig <- compdigamma(y = y, alpha = alpha)
    
    G0 <- (dig + ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) +
    (y - lambda) * alpha * (lambda * alpha + 1) ^ (1 / alpha) - y * alpha) /
    (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) - 1))) *
    alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    # Beta derivative
    G1 <- -((lambda - y) * (alpha * lambda + 1) ^ (1 / alpha) + y) /
    (lambda * (alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) - 1)) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- solve(weight[[x]]) #more stable
      xx %*% uMatrix[x, ]
    })
    
    pseudoResid <- t(pseudoResid)
    dimnames(pseudoResid) <- dimnames(eta)
    pseudoResid
  }

  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
        
        -sum(weight * (lgamma(y + 1 / alpha) - lgamma(1 / alpha) - lgamma(y + 1) - 
        (y + 1 / alpha) * log(1 + lambda * alpha) + y * log(lambda * alpha) - 
        log(1 - (1 + lambda * alpha) ^ (-1 / alpha))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
        
        # log(alpha) derivative
        dig <- compdigamma(y = y, alpha = alpha)
        
        G0 <- (dig + ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) +
        (y - lambda) * alpha * (lambda * alpha + 1) ^ (1 / alpha) - y * alpha) /
        (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) - 1))) *
        alphaLink(eta[, 2], inverse = TRUE, deriv = 1) * weight
        
        # Beta derivative
        G1 <- -((lambda - y) * (alpha * lambda + 1) ^ (1 / alpha) + y) /
        (lambda * (alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) - 1)) *
        lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * weight
        
        if (NbyK) {
          XX <- 1:(attr(X, "hwm")[1])
          return(cbind(as.data.frame(X[1:nrow(eta), XX]) * G1, as.data.frame(X[-(1:nrow(eta)), -XX]) * G0))
        }
        if (vectorDer) {
          return(cbind(G1, G0))
        }
        
        as.numeric(c(G1, G0) %*% X)
      },
      function(beta) {
        lambdaPredNumber <- attr(X, "hwm")[1]
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
        res <- matrix(nrow = length(beta), 
                      ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        trig <- comptrigamma(y = y, alpha = alpha)
        dig <-  compdigamma(y = y, alpha = alpha)
        
        # 2nd log(alpha) derivative
        G00 <- ((trig + ((lambda * alpha + 1) ^ (1 / alpha) * 
        (lambda ^ 2 * alpha ^ 2 + 2 * lambda * alpha + 1) * log(lambda * alpha + 1) ^ 2 +
        ((lambda * alpha + 1) ^ (1 / alpha) * (2 * lambda ^ 2 * alpha ^ 3 + 
        (4 * lambda - 2 * lambda ^ 2) * alpha ^ 2 + (2 - 2 * lambda) * alpha) + 
        (lambda * alpha + 1) ^ (2 / alpha) * (-2 * lambda ^ 2 * alpha ^ 3 - 4 * lambda * alpha ^ 2 - 2 * alpha)) * 
        log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (2 / alpha) * 
        ((3 * lambda ^ 2 - 2 * y * lambda) * alpha ^ 3 +
        (2 * lambda - y) * alpha ^ 2) + (lambda * alpha + 1) ^ (1 / alpha) *
        ((4 * y * lambda - 3 * lambda ^ 2) * alpha ^ 3 +
        (lambda ^ 2 - 2 * lambda + 2 * y) * alpha ^ 2) -
        2 * y * lambda * alpha ^ 3 - y * alpha ^ 2) /
        (alpha ^ 4 * (lambda * alpha + 1) ^ 2 * 
        ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)) * 
        alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
        (dig + ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) +
        (y - lambda) * alpha * (lambda * alpha + 1) ^ (1 / alpha) - y * alpha) /
        (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) - 1))) *
        alphaLink(eta[, 2], inverse = TRUE, deriv = 2))
        
        # mixed derivative
        G01 <- -((alpha * lambda + 1) ^ (1 / alpha + 1) *
        log(alpha * lambda + 1) + (alpha * lambda + 1) ^ (1 / alpha) *
        ((alpha ^ 2 - alpha) * lambda - 2 * alpha ^ 2 * y) +
        (alpha * lambda + 1) ^ (2 / alpha) * (alpha ^ 2 * y - alpha ^ 2 * lambda) + alpha ^ 2 * y) /
        (alpha ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)
        
        G01 <- G01 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
               lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
        
        # second beta derivative
        
        G11 <- (((alpha * lambda + 1) ^ (2 / alpha) * (alpha * lambda ^ 2 - 2 * alpha * y * lambda - y) +
        (alpha * lambda + 1) ^ (1 / alpha) * ((1 - alpha) * lambda ^ 2 + 4 * alpha * y * lambda + 2 * y) - 2 * alpha * y * lambda - y) /
        (lambda ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2) *
        lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 -
        ((lambda - y) * (alpha * lambda + 1) ^ (1 / alpha) + y) /
        (lambda * (alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) - 1)) *
        lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)) * weight
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- 
          t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]) *
            G00 * weight) %*% X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
        
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- 
          t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * 
            G11 * weight) %*% X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- 
          t(t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * 
            G01 * weight) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]))
        
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- 
          t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * 
            G01 * weight) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        res
      }
    )
  }

  validmu <- function(mu) {
    all(is.finite(mu)) && all(0 < mu)
  }

  devResids <- function (y, eta, wt, ...) {
    # TODO:: implement theoritical results
    0
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    N <- pw / (1 - (1 + alpha * lambda) ^ (- 1 / alpha))
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    pr <- 1 - (1 + alpha * lambda) ^ (- 1 / alpha)

    bigTheta1 <- pw  *  alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
    ((lambda * alpha + 1) ^ (1 / alpha - 1) * ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha)) /
    (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)# w.r to alpha
    bigTheta2 <- -pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
    ((alpha * lambda + 1) ^ (1 / alpha - 1) / ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)# w.r to lambda

    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)

    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * (1 - pr) / (pr ^ 2))

    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)

    switch (type,
      "trunc" = {
        stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / 
        (1 - (1 + alpha * lambda) ^ (-1 / alpha))
      },
      "nontrunc" = stats::dnbinom(x = x, mu = lambda, size = 1 / alpha)
    )
  }

  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    lb <- stats::pnbinom(lower, mu = lambda, size = 1 / alpha)
    ub <- stats::pnbinom(upper, mu = lambda, size = 1 / alpha)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qnbinom(p_u, mu = lambda, size = 1 / alpha)
    sims
  }
  
  getStart <- expression(
    start <- stats::glm.fit(
      x = variables[wch$reg, 1:attr(Xvlm, "hwm")[1]],
      y = observed[wch$reg],
      family = stats::poisson(),
      weights = priorWeights[wch$reg],
      ...
    )$coefficients,
    if (attr(family$links, "linkNames")[1] == "neglog") start <- -start,
    if (controlModel$alphaFormula == ~ 1) {
      start <- c(start, log(abs(mean(observed[wch$reg] ^ 2) - mean(observed[wch$reg])) / (mean(observed[wch$reg]) ^ 2 + .25)))
    } else {
      cc <- colnames(Xvlm)
      cc <- cc[grepl(x = cc, pattern = "alpha$")]
      cc <- unlist(strsplit(x = cc, ":alpha"))
      cc <- sapply(cc, FUN = function(x) {
        ifelse(x %in% names(start), start[x], 0) # TODO: gosh this is terrible pick a better method
      })
      if (attr(family$links, "linkNames")[1] == attr(family$links, "linkNames")[2])
        start <- c(start, cc)
      else
        start <- c(start, -cc)
    }
  )
  
  structure(
    list(
      makeMinusLogLike = minusLogLike,
      densityFunction  = dFun,
      links     = links,
      mu.eta    = mu.eta,
      valideta  = function (eta) {TRUE},
      variance  = variance,
      Wfun      = Wfun,
      funcZ     = funcZ,
      devResids = devResids,
      validmu   = validmu,
      pointEst  = pointEst,
      popVar    = popVar,
      family    = "ztnegbin",
      etaNames  = c("lambda", "alpha"),
      simulate  = simulate,
      getStart  = getStart,
      extraInfo = c(
        mean       = "lambda",
        variance   = "lambda * (1 + alpha * lambda)",
        popSizeEst = "1 / (1 - (1 + alpha * lambda) ^ (- 1 / alpha))",
        meanTr     = "lambda / (1 - (1 + alpha * lambda) ^ (-1 / alpha))",
        varianceTr =
          "(lambda + alpha * (lambda ^ 2) - alpha * (lambda ^ 2) * (1 + alpha * lambda) ^ (-1 / alpha)) / ((1 - (1 + alpha * lambda) ^ (-1 / alpha)) ^ 2)"
      )
    ),
    class = c("singleRfamily", "family")
  )
}

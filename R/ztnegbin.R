#' @rdname singleRmodels
#' @importFrom stats uniroot
#' @importFrom stats dnbinom
#' @importFrom rootSolve multiroot
#' @export
ztnegbin <- function(nSim = 1000, epsSim = 1e-8, ...) {
  # Fist for lambda second for alpha
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2])), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),exp(x[,2])), ncol = 2, dimnames = dimnames(x))}

  mu.eta <- function(eta, type = "trunc", ...) {
    A <- invlink(eta)
    lambda <- A[, 1]
    A <- A[, 2]
    switch (type,
    nontrunc = lambda,
    trunc = lambda / (1 - (1 + A * lambda) ^ (-1 / A))
    )
  }

  variance <- function(eta, type = "nontrunc", ...) {
    A <- invlink(eta)
    lambda <- A[, 1]
    A <- A[, 2]
    P0 <- (1 + A * lambda) ^ (-1 / A)
    switch (type,
    nontrunc = lambda * (1 + A * lambda),
    trunc = (lambda + A * (lambda ^ 2) - A * (lambda ^ 2) * P0) / ((1 - P0) ^ 2)
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
    (alpha ^ 2 * (1 - y) + 2 * alpha * digamma(y + 1 / alpha) + trigamma(y + 1 / alpha) - 2 * alpha * digamma(1 + 1 / alpha) - trigamma(1 + 1 / alpha)) / (alpha ^ 4)
  }
  
  # Computing the expected value of di/trigamma functions on (y + 1/alpha)
  
  compExpect <- function(eta) {
    alpha <- invlink(matrix(eta, ncol = 2))
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    res <- res1 <- 0
    k <- 0
    finished <- c(FALSE, FALSE)
    while ((k < nSim) & !all(finished)) {
      k <- k + 1 # 1 is the first possible y value for 0 truncated distribution
      prob <- stats::dnbinom(x = k, size = 1 / alpha, mu = lambda) / (1 - P0)
      toAdd <- c(compdigamma(y = k, alpha = alpha) * prob, comptrigamma(y = k, alpha = alpha) * prob)
      res <- res + toAdd
      finished <- abs(toAdd) < epsSim
    }
    res
  }
  
  Wfun <- function(prior, eta, ...) {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    z <- 1 / alpha
    M <- ((1 + lambda / z) ^ z) - 1
    S <- 1 / (1 + lambda / z)
    G <- 1 / (1 - (1 + lambda / z) ^ (-z))
    Ey <- mu.eta(eta = eta)
    T3 <- (log(S) * (z ^ 2) + z * lambda * S) * (S ^ z) * (S ^ (1 + z)) * (G ^ 2) + G * (S ^ (1 + z)) * (lambda * (1 + z) * S + log(S) * (z ^ 2))
    Edig <- apply(X = eta, MARGIN = 1, FUN = function(x) {compExpect(x)})
    Etrig <- Edig[2,]
    Edig <- Edig[1,]
    matrix(
      c(-lambda * (((1/S ^ z) * (lambda - 1) + 1) * (S ^ 2) / ((1/S ^ z - 1) ^ 2) - (1 + Ey / z) * (S ^ 2)) * prior, # lambda predictor derivative without X matrix,
        -(-lambda * (Ey - lambda) * (S ^ 2) + lambda * T3) * prior * alpha, # mixed derivative without X matrix
        -(-lambda * (Ey - lambda) * (S ^ 2) + lambda * T3) * prior * alpha, # mixed derivative without X matrix
        -((2 * lambda * S * (z ^ 2) + 2 * log(S) * (z ^ 3) + Etrig + # alpha predictor derivative without X matrix
        (Ey + z) * (lambda ^ 2) * (S ^ 2) + (z ^ 3) * 2 * (lambda / z + log(S) / S) * S / M +
        (z ^ 2) * (S ^ (1 - z)) * (z * lambda * S + log(S) * (z ^ 2)) *
        (lambda / z + log(S) / S) / (M ^ 2) + (z ^ 2) * lambda * log(1/S) * S / M +
        (z ^ 2) * lambda * (S ^ 2) * (lambda / z + log(S) / S) / M) * (alpha ^ 2)+
        (log(1/S) * (z ^ 2) + Ey * z - (Ey + z) * lambda * S + Edig +
        G * (S ^ z) * (log(1/S) * (z ^ 2) - z * lambda * S)) * alpha) * prior
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "alpha")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    z <- 1 / alpha
    S <- 1 / (1 + lambda / z)
    G <- 1 / (1 - (1 + lambda / z) ^ (-z))
    dig <- compdigamma(y = y, alpha = alpha)

    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    uMatrix <- matrix(
      c((y + (lambda - y) * S ^ (-z)) * S / (1 - S ^ (-z)),
        (log(1 / S) * (z ^ 2) + y * z - (y + z) * lambda * S + dig +
        G * (S ^ z) * (log(1 / S) * (z ^ 2) - z * lambda * S)) * alpha),
      ncol = 2)
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      #xx <- solve(weight[[x]]) #more stable
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
    deriv <- deriv + 1 # to make it comfort to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        alpha <- invlink(eta)
        lambda <- alpha[, 1]
        alpha <- alpha[, 2]
        z <- 1 / alpha
        
        -sum(weight * (lgamma(y + z) - lgamma(z) - log(factorial(y)) - 
        (y + z) * log(1 + lambda * alpha) + y * log(lambda * alpha) - 
        log(1 - (1 + lambda * alpha) ^ (-z))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        alpha <- invlink(eta)
        lambda <- alpha[, 1]
        alpha <- alpha[, 2]
        # z is inverse of alpha in documentation
        z <- 1 / alpha
        S <- 1 / (1 + lambda / z)
        
        # log(alpha) derivative
        dig <- compdigamma(y = y, alpha = alpha)
        G0 <- t((log(1 / S) * (z ^ 2) + y * z - (y + z) * lambda * S + dig + 
        (1 / (1 - (1 + lambda / z) ^ (-z))) * (S ^ z) * (log(1 / S) * (z ^ 2) - 
        z * lambda * S)) * alpha)
        
        # Beta derivative
        G1 <- t(((y + (lambda - y) / (S ^ z)) * S / (1 - S ^ (-z)))  * weight)
        
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
        alpha <- invlink(eta)
        lambda <- alpha[, 1]
        alpha <- alpha[, 2]
        # z is inverse of alpha in documentation
        z <- 1 / alpha
        M <- ((1 + lambda / z) ^ z) - 1
        S <- 1 / (1 + lambda / z)
        G <- 1 / (1 - (1 + lambda / z) ^ (-z))
        res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
        
        trig <- comptrigamma(y = y, alpha = alpha)
        dig <- compdigamma(y = y, alpha = alpha)
        
        # 2nd log(alpha) derivative
        G00 <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]) *
        ((2 * lambda * S * (z ^ 2) + 2 * log(S) * (z ^ 3) + trig +
        (y + z) * (lambda ^ 2) * (S ^ 2) + (z ^ 3) * 2 * (lambda / z + log(S) / S) * S / M +
        (z ^ 2) * (S ^ (1 - z)) * (z * lambda * S + log(S) * (z ^ 2)) *
        (lambda / z + log(S) / S) / (M ^ 2) + (z ^ 2) * lambda * log(1/S) * S / M +
        (z ^ 2) * lambda * (S^2) * (lambda / z + log(S) / S) / M) * (alpha ^ 2) +
        (log(1/S) * (z ^ 2) + y * z - (y + z) * lambda * S + dig +
        G * (S ^ z) * (log(1/S) * (z ^ 2) - z * lambda * S)) * alpha) * weight) %*% X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
        # mixed derivative
        T1 <- lambda * (y - lambda) * (S^2)
        
        T3 <- (G ^ 2) * (S ^ z) * (log(S) * (z ^ 2) + z * lambda * S) * (S ^ (1 + z)) + (S ^ (1 + z)) * (lambda * (1 + z) * S + log(S) * (z ^ 2)) * G
        
        G01 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * as.numeric(-T1 + lambda * T3) * alpha * weight) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        # second beta derivative
        
        G11 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * lambda * (((1/S ^ z) * (lambda - 1) + 1) * (S ^ 2) / ((1/S ^ z - 1) ^ 2) - (1 + y / z) * (S^2)) * weight) %*% X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(G01)
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- G01
        
        res
      }
    )
  }

  validmu <- function(mu) {
    all(is.finite(mu)) && all(0 < mu)
  }

  devResids <- function (y, eta, wt, ...) {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    
    logLikFit <- wt * (lgamma(y + 1/alpha) - lgamma(1/alpha) -
    log(factorial(y)) - (y + 1/alpha) * log(1+alpha * lambda) +
    y * log(lambda * alpha) - log(1 - (1+alpha * lambda) ^ (-1/alpha)))
    
    yUnq <- unique(y) # see comments in zotpoisson
    findL <- function(yNow) {
      rootSolve::multiroot(
        start = c(.5, log(yNow), 1),# maybe pick better starting points
        f = function(x) { # TODO:: provide analytic jacobian matrix will make it faster and more reliable
          s <- x[1] # this is the lagrange multiplier and has no constraints of positivity
          l <- exp(x[2])
          a <- exp(x[3]) # including constraints
          prob <- 1 - (1+a*l)^(-1/a)
          prob <- 1 / prob
          c(l*prob - yNow,# s der
          yNow/l+(1+yNow*a)/(1+l*a)+s*yNow*l-(yNow-l)/(l*(1+l*a))-s*yNow*(yNow-l)/(l*(1+l*a)),# lambda der
          (1+yNow*a)*l+(yNow-l)*(1+s*yNow)*((1+a*l)*log(1+a*l)-a*l)+l*(1+a*l)*(digamma(1/a)+log(a)+log(yNow+1/a))-(digamma(yNow+1/a)+1)*(1+a*l)*l)#alpha der
        }
      )$f.root
    }
    logLikIdeal <- sapply(yUnq, FUN = function(x) {
      ifelse(y[x] == 1, 0, findL(x))
    })
    
    logLikIdeal <- sapply(y, FUN = function(x) logLikIdeal[yUnq == x])
    
    sign(y - mu.eta(eta = eta)) * sqrt(-2 * wt * (logLikFit - logLikIdeal))
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    disp <- invlink(eta)
    lambda <- disp[, 1]
    disp <- disp[, 2]
    pr <- 1 - (1 + disp * lambda) ^ (- 1 / disp)
    N <- pw / pr
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    z <- invlink(eta)
    lambda <- z[, 1]
    z <- z[, 2]
    pr <- 1 - (1 + z * lambda) ^ (- 1 / z)
    S <- 1 / (1 + z * lambda)

    bigTheta1 <- pw * z * ((S ^ (1 - 1 / z)) * (log(1 / S) / S - z * lambda) / ((z ^ 2) * ((1 - (1 / S) ^ (1 / z)) ^ 2))) # w.r to alpha
    bigTheta2 <- -(pw * as.numeric(lambda * (S ^ (1 - 1 / z)) / ((1 - (1 / S) ^ (1 / z)) ^ 2))) # w.r to lambda

    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)

    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * (1 - pr) / (pr ^ 2))

    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / (1 - P0)
  }

  simulate <- function(n, eta, lower = 0, upper = Inf) {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
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
    if (!is.null(controlMethod$alphaStart)) {
      start <- c(start, controlMethod$alphaStart)
    } else {
      if (controlModel$alphaFormula == ~ 1) {
        start <- c(start, log(abs(mean(observed[wch$reg] ^ 2) - mean(observed[wch$reg])) / (mean(observed[wch$reg]) ^ 2 + .25)))
      } else {
        cc <- colnames(Xvlm)
        cc <- cc[grepl(x = cc, pattern = "alpha$")]
        cc <- unlist(strsplit(x = cc, ":alpha"))
        cc <- sapply(cc, FUN = function(x) {
          ifelse(x %in% names(start), start[x], 0) # TODO: gosh this is terrible pick a better method
        })
        start <- c(start, cc)
      }
    }
  )
  
  structure(
    list(
      makeMinusLogLike = minusLogLike,
      linkfun = link,
      linkinv = invlink,
      mu.eta = mu.eta,
      link = c("log", "log"),
      valideta = function (eta) {TRUE},
      variance = variance,
      Wfun = Wfun,
      funcZ = funcZ,
      devResids = devResids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      simulate = simulate,
      family = "ztnegbin",
      parNum = 2,
      etaNames = c("lambda", "alpha"),
      densityFunction = dFun,
      getStart = getStart
    ),
    class = "family"
  )
}

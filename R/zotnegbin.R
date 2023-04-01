#' @rdname singleRmodels
#' @export
zotnegbin <- function(nSim = 1000, epsSim = 1e-8, ...) {
  # Fist for lambda second for alpha
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2])), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),exp(x[,2])), ncol = 2, dimnames = dimnames(x))}

  mu.eta <- function(eta, type = "trunc", ...) {
    A <- exp(eta[, 2])
    lambda <- exp(eta[, 1])
    switch (type,
    "nontrunc" = lambda,
    "trunc" = (lambda - lambda * ((1 + A * lambda) ^ (-1-1/A))) / (1 - (1 + A * lambda) ^ (-1/A) - lambda * ((1 + A * lambda) ^ (-1-1/A)))
    )
  }

  variance <- function(eta, type = "nontrunc", ...) {
    A <- exp(eta[, 2])
    lambda <- exp(eta[, 1])
    switch (type,
            "nontrunc" = lambda * (1 + A * lambda),
            "trunc" = {
              EX2 <- lambda * (1 + A * lambda) + lambda ^ 2;
              prob <- (1 - (1 + A * lambda) ^ (-1 / A) - lambda * ((1 + A * lambda) ^ (-1 - 1 / A)));
              term <- lambda * ((1 + A * lambda) ^ (-1 - 1 / A));
              EX <- (lambda - lambda * ((1 + A * lambda) ^ (-1 - 1 / A))) / (1 - (1 + A * lambda) ^ (-1 / A) - lambda * ((1 + A * lambda) ^ (-1 - 1 / A)));
              (EX2 - term) / prob - EX ^ 2
            }
    )
  }
  
  compExpect <- function (eta) {
    alpha <- invlink(matrix(eta, ncol = 2))
    lambda <- alpha[,1]
    alpha <- alpha[,2]
    z <- 1 / alpha
    res <- res1 <- 0
    k <- 1
    repeat{
      k <- k + 1 # 2 is the first possible y value for 0 truncated distribution
      prob <- stats::dnbinom(x = k, size = 1 / alpha, mu = lambda) / (1 - stats::dnbinom(x = 0, size = 1 / alpha, mu = lambda) - stats::dnbinom(x = 1, size = 1 / alpha, mu = lambda))
      toAdd <- (digamma(k + z) - digamma(z)) * z * prob
      toAdd1 <- (trigamma(k + z) - trigamma(z)) * (z ^ 2) * prob
      res <- res + toAdd
      res1 <- res1 + toAdd1
      if ((k == nSim) | ((abs(toAdd) < 1e-8) & (abs(toAdd1) < 1e-8))) {
        break
      }
    }
    c(res, res1)
  }
  
  Wfun <- function(prior, eta, ...) {
    alpha <- exp(eta[, 2])
    lambda <- exp(eta[, 1])
    z <- 1 / alpha
    M <- 1 + lambda * alpha
    S <- 1 / M
    Ey <- mu.eta(eta = eta)
    prob <- 1 - S ^ z - lambda * (S ^ (1 + z))
    # The following vectors are meant to decrease compitation time
    Edig <- apply(X = eta, MARGIN = 1, FUN = function(x) {compExpect(x)})
    Etrig <- Edig[2,]
    Edig <- Edig[1,]
    G11 <- lambda * ((-(alpha ^ 3) * Ey * (lambda ^ 2) * ((M^z - 1) ^ 2) + (alpha ^ 2) * lambda * ((lambda ^ 2) * (M ^ z) - lambda * ((M ^ z) - 1) * (1 - 2 * Ey + M ^ z) - 2 * Ey * ((M^z - 1) ^ 2)) + (lambda ^ 2) * (M ^ z) + alpha * ((lambda ^ 3) * (M ^ z) + (lambda ^ 2) * (1 - Ey + M^z) - 2 * lambda * (M^z - 1) * (M^z - Ey) - Ey * ((M^z - 1) ^ 2)) - ((M^z - 1) ^ 2)) / ((M ^ 2) * ((M^z + lambda * (alpha * (M^z - 1) - 1) - 1) ^ 2))) * prior
    G01 <- lambda * (-alpha * (Ey + z) * S + lambda * (alpha ^ 2) * (Ey + z) * (M ^ -2) + ((-1 - z) * alpha * lambda / (M ^ (2 + z))) / prob + (lambda / (M ^ (2 + z))) / prob + S - alpha * lambda * (1 + z) * (alpha * lambda * (-2 - z) * S + log(M) * z) / (prob * (M ^ (2 + z))) - (alpha * (-1 - z) * lambda * (1 / (M ^ (2 + z))) * ((M^-z) * (z * log(S) + lambda * S) - lambda * (1 / (M ^ (1 + z))) * (z * log(M) + S * lambda * alpha * (-1 - z)))) / (prob ^ 2)) * prior
    
    G00 <- (Etrig + Edig + z * (((M ^ (2 + z)) * log(M) * (log(M) * (-z) + lambda * S * (2 * alpha + 1))) +
    (Ey - lambda) * alpha * (M ^ (1 + z)) * (log(M) * (1 - z * (1 + alpha)) + lambda * S * (1 + alpha)) - 2 * lambda * (Ey + lambda) * (alpha ^ 2) + lambda * alpha * (M ^ (z+1)) +
    (Ey - lambda) * alpha * (M ^ (1 + z)) - (1 + lambda) * Ey * alpha) / (M * ((M ^ (1 + z)) - alpha * lambda - lambda - 1)) - z * (log(M) * (M ^ (2 + z)) - lambda * (Ey + lambda) * (alpha ^ 2) + (Ey - lambda) * alpha * (M ^ (1 + z)) - (1 + lambda) * Ey * alpha) * ((M ^ (1 + z)) * ((1 - z * (1 + alpha)) * log(M) + lambda * S * (1 + alpha)) - lambda * alpha) /
    (M * (((M ^ (1 + z)) - alpha * lambda - lambda - 1) ^ 2)) - z * (log(M) * (M ^ (2 + z)) - lambda * (Ey + lambda) * (alpha ^ 2) + (Ey - lambda) * alpha * (M ^ (1 + z)) - (1 + lambda) * Ey * alpha) /
    (M * ((M ^ (1 + z)) - alpha * lambda - lambda - 1)) - lambda * (log(M) * (M ^ (2 + z)) - lambda * (Ey + lambda) * (alpha ^ 2) + (Ey - lambda) * alpha * (M ^ (1 + z)) - (1 + lambda) * Ey * alpha) / ((M ^ 2) * ((M ^ (1 + z)) - alpha * lambda - lambda - 1))) * prior
    
    matrix(
      -c(G11, # lambda predictor derivative without X matrix,
        G01,  # mixed derivative without X matrix
        G01,  # mixed derivative without X matrix
        G00   # alpha predictor derivative without X matrix
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "alpha")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    alpha <- exp(eta[, 2])
    lambda <- exp(eta[, 1])
    z <- 1 / alpha
    M <- 1 + lambda * alpha
    S <- 1 / M
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    uMatrix <- matrix(
      c((y - alpha * (y + z) * lambda * S + lambda * alpha * (-1 - z) * lambda * (S ^ (2 + z)) / (1 - S ^ z - lambda * (S ^ (1 + z)))),
        (-(digamma(y + z) - digamma(z)) * (z ^ 2) - 
        z * ((S ^ z) * (z * log(S) + lambda * S) - lambda *
        (S ^ (1 + z)) * (z * log(M) - alpha * (1 + z) * lambda * S)) / (1 - S ^ z - lambda * (S ^ (1 + z))) +
        z * log(M) * z - lambda * S * (y + z) + y * z) * alpha),
      ncol = 2)
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- solve(weight[[x]]) # more stable
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
    deriv <- deriv + 1 # to make it comfort to how swith in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        alpha <- exp(eta[, 2])
        lambda <- exp(eta[, 1])
        z <- 1 / alpha
        S <- 1 / (1 + lambda * alpha)

        -sum(weight * (lgamma(y + z) - lgamma(z) -
        log(factorial(y)) + (y + z) * log(S) +
        y * log(lambda * alpha) - log(1 - S ^ z - lambda * (S ^ (1 + z)))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        alpha <- exp(eta[, 2])
        lambda <- exp(eta[, 1])
        # z is inverse of alpha in documentation
        z <- 1 / alpha
        S <- (1 + lambda * alpha) ^ -1
        
        # log(alpha) derivative
        G0 <- t((-(digamma(y + z) - digamma(z)) * (z ^ 2)  - 
        z * ((S ^ z) * (z * log(S) + lambda * S) - lambda *
        (S ^ (1 + z)) * (-z * log(S) - alpha * (1 + z) * lambda * S)) / (1 - S ^ z - lambda * (S ^ (1 + z))) -
        z * log(S) * z - lambda * S * (y + z) + y * z) * alpha * weight)
        # Beta derivative
        G1 <- (weight * (y - alpha * (y + z) * lambda * S - lambda * alpha * (1 + z) * lambda * (S ^ (2 + z)) / (1 - S ^ z - lambda * (S ^ (1 + z)))))
        
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
        alpha <- exp(eta[, 2])
        lambda <- exp(eta[, 1])
        # z is inverse of alpha in documentation
        z <- 1 / alpha
        M <- 1 + lambda * alpha
        S <- 1 / M
        prob <- 1 - S ^ z - lambda * (S ^ (1 + z))
        res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
        
        # 2nd log(alpha) derivative
        
        G00 <- weight * ((trigamma(y + z) - trigamma(z)) * (z ^ 2) + 
        (digamma(y + z) - digamma(z)) * z + z * (((M ^ (z * (2 * alpha + 1))) * log(M) * (-log(M) * z + lambda * S * (2 + z))) +
        (y - lambda) * alpha * (M ^ (1 + z)) * (log(M) * (1 - z * (1 + alpha)) + 
        lambda * S * (1 + alpha)) - 2 * lambda * (y + lambda) * (alpha ^ 2) + lambda * alpha * (M ^ (z + 1)) +
        (y - lambda) * alpha * (M ^ (1 + z)) - (1 + lambda) * y * alpha) / (M * (M ^ (1 + z) - alpha * lambda - lambda - 1)) - 
        z * (log(M) * (M ^ (2 + z)) - lambda * (y + lambda) * (alpha ^ 2) + (y - lambda) * alpha * (M ^ (1 + z)) - 
        (1 + lambda) * y * alpha) * ((M ^ (1 + z)) * ((1 - z * (1 + alpha)) * log(M) + lambda * S * (1 + alpha)) - lambda * alpha) /
        (M * ((M ^ (1 + z) - alpha * lambda - lambda - 1) ^ 2)) - z * (log(M) * (M ^ (2 + z)) - lambda * (y + lambda) * (alpha ^ 2) + 
        (y - lambda) * alpha * (M ^ (1 + z)) - (1 + lambda) * y * alpha) /
        (M * (M ^ (1 + z) - alpha * lambda - lambda - 1)) - lambda * (log(M) * (M ^ (2 + z)) - 
        lambda * (y + lambda) * (alpha ^ 2) + (y - lambda) * alpha * (M ^ (1 + z)) - 
        (1 + lambda) * y * alpha) / ((M ^ 2) * (M ^ (1 + z) - alpha * lambda - lambda - 1)))
          
        
        # mixed derivative
        G01 <- lambda * weight * (-alpha * (y + z) * S + lambda * (alpha ^ 2) * (y + z) * (M ^ -2) -
        ((1 + z) * alpha * lambda / (M ^ (2 + z))) / prob + (lambda / (M ^ (2 + z))) / prob + S -
        alpha * lambda * (1 + z) * (alpha * lambda * (-2 - z) * S + log(M) * z) / (prob * (M ^ (2 + z))) -
        (alpha * (-1 - z) * lambda *
        ((M ^ (-z)) * (z * log(S) + lambda * S) - lambda * (M ^ (-1 - z)) *
        (z * log(M) - S * lambda * alpha * (1 + z)))) / ((prob ^ 2) * (M ^ (2 + z))))
        
        # second beta derivative
        G11 <- lambda * weight * ((-(alpha ^ 3) * y * (lambda ^ 2) * ((M ^ z - 1) ^ 2) +
        (alpha ^ 2) * lambda * ((lambda ^ 2) * (M ^ z) - lambda * (M ^ z - 1) *
        (1 - 2 * y + M ^ z) - 2 * y * ((M ^ z - 1) ^ 2)) +
        (lambda ^ 2) * (M ^ z) + alpha * ((lambda ^ 3) * (M ^ z) +
        (lambda ^ 2) * (1 - y + M ^ z) - 2 * lambda * (M ^ z - 1) *
        (M ^ z - y) - y * ((M ^ z - 1) ^ 2)) - (M ^ z - 1) ^ 2) /
        ((M ^ 2) * ((M ^ z + lambda * (alpha * (M ^ z - 1) - 1) - 1) ^ 2)))
        
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]) * G00) %*% X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * G11) %*% X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * G01) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]))
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * G01) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        res
      }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  devResids <- function (y, eta, wt, ...) {
    # lambda <- exp(eta[, 1])
    # alpha <- exp(eta[, 2])
    
    # I won't be able to do this before eurostat conference it has to be delayed
    # 
    # logLikFit <- (lgamma(y + 1/alpha) - lgamma(1/alpha) -
    # log(factorial(y)) - (y + 1/alpha) * log(1+alpha * lambda) +
    # y * log(lambda * alpha) - log(1 - (1+alpha * lambda) ^ (-1/alpha) -
    # lambda * ((1+alpha * lambda) ^ (-1-1/alpha))))
    # 
    # 
    # # TODO:: alpha goes to 0 check if this is not analytic
    # yUnq <- unique(y) # see comments in zotpoisson
    # findL <- function(yNow) {
    #   root <- rootSolve::multiroot(
    #     start = c(1, mean(lambda), mean(alpha)),# maybe pick better starting points
    #     f = function(x) {# TODO:: provide analytic jacobian matrix will make it faster and more reliable
    #       s <- log(x[1]) # this is the lagrange multiplier and has no constraints of positivity
    #       l <- x[2]
    #       a <- x[3] # including constraints
    #       
    #       # musi być gdzieś błąd w pochodnej tutaj
    #       # może daj jakobian i będzie stabilniej
    #       # c(log(l-l*((1+a*l)^(-1-1/a)))-log(1-(1+a*l)^(-1/a)-l*((1+a*l)^(-1-1/a)))-log(yNow),#sder
    #       #   yNow/l-(1+a*yNow)/(1+a*l)+s/l+s*(1+a)*((1+a*l)^(-2-1/a))/(1-(1+a*l)^(-1-1/a))-(1+s)*l*(1+a)*((1+a*l)^(-2-1/a))/(1-(1+a*l)^(-1/a)-l*((1+a*l)^(-1-1/a))),#lambda der
    #       #   (digamma(1/a)-digamma(yNow+1/a))/(a^2)-l*(yNow+1/a)/(1+a*l)+log(1+a*l)/(a^2)+
    #       #   yNow/a-s*((1+a*l)^(-1-1/a))*(log(1+a*l)/(a^2)-l*(1+1/a)/(1+a*l))/(1-(1+a*l)^(-1-1/a))+
    #       #   (1+s)*((log(1+a*l)/(a^2)-l/(a*(1+a*l)))/((1+a*l)^(1/a))+l*((1+a*l)^(-1-1/a))*
    #       #   (log(1+a*l)/(a^2)-l*(1+1/a)/(1+a*l)))/(1-(1+a*l)^(-1/a)-l*((1+a*l)^(-1-1/a))))#alpha der
    #       c(log(l-l*(a*l+1)^(-1/a-1))-log(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)-log(yNow),#sder
    #         -((-1/a-1)*a*(-s-1)*l*(a*l+1)^(-1/a-2))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)+(s*(-(a*l+1)^(-1/a-1)-(-1/a-1)*a*l*(a*l+1)^(-1/a-2)+1))/(l-l*(a*l+1)^(-1/a-1))+(a*(-yNow-1/a))/(a*l+1)+yNow/l,#lambda der
    #         ((-s-1)*(-(log(l*a+1)/a^2-l/(a*(l*a+1)))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1))))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)-(l*s*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1)))/(l-l*(l*a+1)^(-1/a-1))+(log(l*a+1)-digamma(1/a+yNow)+digamma(1/a))/a^2+(l*(-1/a-yNow))/(l*a+1)+yNow/a)
    #     }, maxiter = 10000, positive = TRUE
    #   )$root
    #   print(root)
    #   root <- log(root)
    #   
    #   (lgamma(yNow + exp(-root[3])) - lgamma(exp(-root[3])) -
    #   log(factorial(yNow)) - (yNow + exp(-root[3])) * log(1+exp(root[2] + root[3])) +
    #   yNow * (root[2]+root[3]) - log(1 - (1+exp(root[2] + root[3])) ^ (-exp(-root[3])) -
    #   exp(root[2]) * ((1+exp(root[2] + root[3])) ^ (-1-exp(-root[3])))))
    # }
    # logLikIdeal <- sapply(yUnq, FUN = function(x) {
    #   ifelse(x == 2, 0, findL(x))
    # })
    # 
    # logLikIdeal <- sapply(y, FUN = function(x) logLikIdeal[yUnq == x])
    # 
    # sign(y - mu.eta(eta = eta)) * sqrt(-2 * wt * (logLikFit - logLikIdeal))
    NULL
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    z <- exp(-eta[, 2])
    lambda <- exp(eta[, 1])
    S <- 1 / (1 + lambda / z)
    prob <- 1 - S ^ z - lambda * (S ^ (1 + z))
    N <- (pw * (1 - lambda * (S ^ (1 + z))) / prob)
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    alpha <- exp(eta[, 2])
    lambda <- exp(eta[, 1])
    z <- 1 / alpha
    M <- 1 + lambda / z
    S <- 1 / M
    prob <- 1 - S ^ z - lambda * (S ^ (1 + z))
    
    bigTheta1 <- (pw * alpha *  as.numeric(
                  (prob * lambda * (S ^ (2 + z)) * 
                  (alpha * (alpha + 1) * lambda -M * log(M)) -
                  (1 - lambda * (S ^ (1 + z))) *(-lambda * (S ^ (1 + z)) *
                  (log(M) * (z ^ 2) - (1 + z) * lambda * S * z) -
                  (S ^ z) * (log(M) * (z ^ 2) - lambda * z * S))) /
                  (prob ^ 2)))
    
    bigTheta2 <- (pw * as.numeric(lambda *
                  (prob * (lambda - 1) * (S ^ (2 + z)) -
                  (1 + alpha) * lambda * (S ^ (2 + z)) *
                  (1 - lambda * (S ^ (1 + z)))) / (prob ^ 2)))
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * ((1 - lambda * (S ^ (1 + z))) ^ 2) *
              (1 - prob) / (prob ^ 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    alpha <- exp(eta[, 2])
    lambda <- exp(eta[, 1])
    stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / (1 - stats::dnbinom(x = 0, mu = lambda, size = 1 / alpha) - stats::dnbinom(x = 1, mu = lambda, size = 1 / alpha))
  }

  simulate <- function(n, eta, lower = 0, upper = Inf) {
    alpha <- exp(eta[, 2])
    lambda <- exp(eta[, 1])
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
      popVar = popVar,
      simulate = simulate,
      family = "zotnegbin",
      etaNames = c("lambda", "alpha"),
      densityFunction = dFun,
      getStart = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

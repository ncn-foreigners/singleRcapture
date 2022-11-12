#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @export
ztoigeom <- function(...) {
  # Fist for lambda second for omega
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2]/ (1 - x[,2]))), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),1/(exp(-x[,2]) + 1)), ncol = 2, dimnames = dimnames(x))}
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
            "nontrunc" = omega * lambda / (1 + lambda) + (1 - omega) * lambda / ((1 + lambda) ** 2) + (1 - omega) * (2 + lambda) * lambda * lambda / ((1 + lambda) ** 2),
            "trunc" = omega + (1 - omega) / (1 + lambda) + lambda * (2 + lambda) * (1 - omega) / (1 + lambda)
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
            "nontrunc" = omega * lambda / (1 + lambda) + (1 - omega) * lambda / ((1 + lambda) ** 2) + (1 - omega) * (2 * (lambda ** 2) +  5 * lambda + 4) * lambda * lambda / ((1 + lambda) ** 2),
            "trunc" = omega + (1 - omega) / (1 + lambda) + lambda * (2 * (lambda ** 2) + 5 * lambda + 4) * (1 - omega) / (1 + lambda)
    ) - mu.eta(eta = eta) ** 2
  }
  
  Wfun <- function(prior, y, eta, ...) {
    # TODO
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- omega + (1 - omega) / (1 + lambda) # expected for I's
    #z <- ifelse(y == 1, y, 0)
    Ey <- mu.eta(eta = eta)
    
    G1 <- z * (omega - 1) / ((1 + lambda) * (1 + lambda * omega))
    G1 <- G1 - (1 - z) * (lambda + 1 - Ey) / (lambda ** 2 + lambda) # this one
    #G1 <- G1 - (lambda + 1 - Ey) / (lambda ** 2 + lambda) + z * lambda / (lambda ** 2 + lambda)
    G0 <- z * lambda / (1 + lambda * omega) - (1 - z) / (1 - omega)
    
    G00 <- -z * (lambda ** 2) / ((1 + lambda * omega) ** 2) - (1 - z) / ((1 - omega) ** 2)
    G00 <- prior * (G0 * (omega * (1 - omega) * (1 - 2 * omega)) + G00 * ((omega * (1 - omega)) ** 2))  # second derivative of inverse logistic link
    
    # mixed derivative
    G01 <- z / ((1 + lambda * omega) ** 2)
    G01 <- G01 * lambda * omega * (1 - omega) * prior
    
    # Beta^2 derivative
    G11 <- -z * (omega - 1) * (2 * lambda * omega + omega + 1) / (((1 + lambda * omega) ** 2) * ((1 + lambda) ** 2))
    G11 <- G11 + (1 - z) * (lambda ** 2 - 2 * lambda * (Ey - 1) - Ey + 1) / ((lambda ** 2) * ((1 + lambda) ** 2)) # This one
    #G11 <- G11 + (lambda ** 2 - 2 * lambda * (Ey - 1) - Ey + 1) / ((lambda ** 2) * ((1 + lambda) ** 2)) - z * (lambda ** 2) / ((lambda ** 2) * ((1 + lambda) ** 2))
    G11 <- (G11 * lambda * lambda + G1 * lambda) * prior # second derivative of log link
    matrix(
      -c(G11, # lambda
         G01, # mixed
         G01, # mixed
         G00  # omega
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "omega")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- ifelse(y == 1, y, 0)
    G1 <- z * (omega - 1) / ((1 + lambda) * (1 + lambda * omega))
    G1 <- G1 - (1 - z) * (lambda + 1 - y) / (lambda ** 2 + lambda)
    G1 <- G1 * lambda
    G0 <- z * lambda / (1 + lambda * omega) - (1 - z) / (1 - omega) # omega derivative
    G0 <- G0 * omega * (1 - omega)
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      #xx <- solve(weight[[x]])
      xx %*% uMatrix[x, ]
    })
    pseudoResid <- t(pseudoResid)
    dimnames(pseudoResid) <- dimnames(eta)
    pseudoResid
  }
  
  minusLogLike <- function(y, X, weight = 1, ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    z <- as.numeric(y == 1)
    function(beta) {
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      lambda <- invlink(eta)
      omega <- lambda[, 2]
      lambda <- lambda[, 1]
      theta <- exp(eta[, 2])
      -sum(weight * (z * log(omega + (1 - omega) / (1 + lambda)) + (1 - z) * 
      (log(1 - omega) - log(1 + lambda) + (y - 1) * log(lambda) - (y - 1) * log(1 + lambda))))
    }
  }
  
  gradient <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    z <- as.numeric(y == 1)
    function(beta) {
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      lambda <- invlink(eta)
      omega <- lambda[, 2]
      lambda <- lambda[, 1]
      G1 <- z * (omega - 1) / ((1 + lambda) * (1 + lambda * omega))
      G1 <- G1 - (1 - z) * (lambda + 1 - y) / (lambda ** 2 + lambda)
      G1 <- G1 * weight * lambda
      G0 <- z * lambda / (1 + lambda * omega) - (1 - z) / (1 - omega) # omega derivative
      G0 <- G0 * weight * omega * (1 - omega)
      if (NbyK) {
        XX <- sapply(as.data.frame(X[1:nrow(eta), ]), FUN = function(x) {all(x == 0)})
        return(cbind(as.data.frame(X[1:nrow(eta), !(XX)]) * G1, as.data.frame(X[-(1:nrow(eta)), XX]) * G0))
      }
      if (vectorDer) {
        return(cbind(G1, G0))
      }
      as.numeric(c(G1, G0) %*% X)
    }
  }
  
  hessian <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    z <- as.numeric(y == 1)
    function (beta) {
      lambdaPredNumber <- attr(X, "hwm")[1]
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      lambda <- invlink(eta)
      omega <- lambda[, 2]
      lambda <- lambda[, 1]
      Xlambda <- X[1:(nrow(X) / 2), 1:lambdaPredNumber]
      Xomega <- X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
      res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
      
      
      G1 <- z * (omega - 1) / ((1 + lambda) * (1 + lambda * omega))
      G1 <- G1 - (1 - z) * (lambda + 1 - y) / (lambda ** 2 + lambda)
      G0 <- z * lambda / (1 + lambda * omega) - (1 - z) / (1 - omega) # omega derivative
      
      # omega^2 derivative
      G00 <- -z * (lambda ** 2) / ((1 + lambda * omega) ** 2) - (1 - z) / ((1 - omega) ** 2)
      G00 <- weight * (G0 * (omega * (1 - omega) * (1 - 2 * omega)) + G00 * ((omega * (1 - omega)) ** 2))  # second derivative of inverse logistic link
      G00 <- t(as.data.frame(Xomega * G00)) %*% as.matrix(Xomega)
      
      # mixed derivative
      G01 <- z / ((1 + lambda * omega) ** 2)
      G01 <- G01 * lambda * omega * (1 - omega) * weight
      G01 <- t(as.data.frame(Xlambda) * G01) %*% as.matrix(Xomega)
      
      # Beta^2 derivative
      G11 <- -z * (omega - 1) * (2 * lambda * omega + omega + 1) / (((1 + lambda * omega) ** 2) * ((1 + lambda) ** 2))
      G11 <- G11 + (1 - z) * (lambda ** 2 - 2 * lambda * (y - 1) - y + 1) / ((lambda ** 2) * ((1 + lambda) ** 2))
      G11 <- (G11 * lambda * lambda + G1 * lambda) * weight # second derivative of log link
      G11 <- t(as.data.frame(Xlambda * G11)) %*% Xlambda
      res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
      res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
      res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(G01)
      res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- G01
      
      res
    }
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  dev.resids <- function(y, mu, wt, theta, ...) {
    #TODO
    0
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
    lambda <- lambda[, 1]
    N <- pw * (1 + 1 / lambda)
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    
    bigTheta1 <- rep(0, nrow(eta)) # w.r to omega
    bigTheta2 <- -lambda * pw / (lambda ** 2) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * (1 + lambda) / (lambda ** 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    ifelse(x == 1, omega + (1 - omega) / (1 + lambda), 
          (1 - omega) * (lambda ** (x - 1)) / ((1 + lambda) ** x))
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    CDF <- function(x) {
      p <- lambda / (1 + lambda)
      ifelse(x == Inf, 1, ifelse(x < 0, 0, ifelse(x < 1, 1 - p, 1 - p + omega * p +  (1 - omega) * (lambda - lambda * (p ** x)) / (1 + lambda))))
    }
    lb <- CDF(lower)
    ub <- CDF(upper)
    p_u <- stats::runif(n, lb, ub)
    sims <- NULL
    for (k in 1:n) {
      m <- 0
      while(CDF(m) < p_u[k]) {
        m <- m + 1
      }
      sims <- c(sims, m)
    }
    sims
  }
  
  structure(
    list(
      makeMinusLogLike = minusLogLike,
      makeGradient = gradient,
      makeHessian = hessian,
      linkfun = link,
      linkinv = invlink,
      mu.eta = mu.eta,
      link = c("log", "logit"),
      valideta = function (eta) {TRUE},
      variance = variance,
      Wfun = Wfun,
      funcZ = funcZ,
      dev.resids = dev.resids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      family = "ztoigeom",
      parNum = 2,
      etaNames = c("lambda", "omega"),
      densityFunction = dFun,
      simulate = simulate
    ),
    class = "family"
  )
}

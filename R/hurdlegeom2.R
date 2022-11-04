#' @rdname ztoipoisson
#' @importFrom lamW lambertW0
#' @export
Hurdleztgeom <- function(...) {
  # Fist for lambda second for PI
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2]/ (1 - x[,2]))), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),1/(exp(-x[,2]) + 1)), ncol = 2, dimnames = dimnames(x))}
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
            "nontrunc" = PI + (1 - PI) * lambda * lambda * (2 + lambda) / (lambda ** 2 + lambda + 1),
            "trunc" = PI * (lambda ** 2 + lambda + 1) / (lambda ** 2 + PI * (lambda + 1)) + (1 - PI) * lambda * lambda * (2 + lambda) / (lambda ** 2 + PI * (lambda + 1))
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    switch(type,
           "nontrunc" = PI + (1 - PI) * lambda * lambda * (2 * lambda * lambda + 5 * lambda + 4) / (lambda ** 2 + lambda + 1),
           "trunc" = PI * (lambda ** 2 + lambda + 1) / (lambda ** 2 + PI * (lambda + 1)) + (1 - PI) * lambda * lambda * (2 * lambda * lambda + 5 * lambda + 4) / (lambda ** 2 + PI * (lambda + 1))
    ) - (mu.eta(eta = eta, type = type) ** 2)
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- PI * (lambda ** 2 + lambda + 1) / (lambda ** 2 + PI * (lambda + 1))
    Ey <- mu.eta(eta)
    
    G1 <- z * (2 * lambda + 1) / (lambda ** 2 + lambda + 1)
    #G1 <- G1 + (1 - z) * (y / lambda - (y - 1) / (1 + lambda)) #this part
    G1 <- G1 + Ey / lambda - (Ey - 1) / (1 + lambda) - z / lambda
    G1 <- G1 - (2 * lambda + PI) / (lambda ** 2 + PI * (lambda + 1)) # lambda derivative
    G0 <- z / PI - (1 - z) / (1 - PI) - (1 + lambda) / (lambda ** 2 + PI * (1 + lambda)) # PI derivative
    
    # PI^2 derivative
    G00 <- -z / (PI ** 2) - (1 - z) / ((1 - PI) ** 2) + ((lambda + 1) ** 2) / ((lambda ** 2 + PI * (lambda + 1)) ** 2)
    G00 <- prior * (G0 * (PI * (1 - PI) * (1 - 2 * PI)) + G00 * ((PI * (1 - PI)) ** 2))
    
    # mixed
    
    G01 <- lambda * (lambda + 2) / ((lambda ** 2 + PI * (lambda + 1)) ** 2)
    G01 <- G01 * lambda * PI * (1 - PI) * prior
    
    # Beta^2 derivative
    #G11 <- (1 - z) * ((y - 1) / ((1 + lambda) ** 2) - y / (lambda ** 2)) # This part
    G11 <- (Ey - 1) / ((1 + lambda) ** 2) - Ey / (lambda ** 2) + z / (lambda ** 2)
    G11 <- G11 + z * (2 / (lambda ** 2 + lambda + 1) - ((2 * lambda + 1) ** 2) / ((lambda ** 2 + lambda + 1) ** 2))
    G11 <- G11 + ((2 * lambda + PI) ** 2) / ((lambda ** 2 + PI * (lambda + 1)) ** 2) - 2 / (lambda ** 2 + PI * (lambda + 1))
    G11 <- (G11 * lambda * lambda + G1 * lambda) * prior
    
    matrix(
      -c(G11, # lambda
         G01, # mixed
         G01, # mixed
         G00  # pi
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "pi")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- ifelse(y == 1, y, 0)
    G1 <- z * (2 * lambda + 1) / (lambda ** 2 + lambda + 1)
    G1 <- G1 + (1 - z) * (y / lambda - (y - 1) / (1 + lambda))
    G1 <- G1 - (2 * lambda + PI) / (lambda ** 2 + PI * (lambda + 1))
    G1 <- G1 * lambda # lambda derivative
    G0 <- z / PI - (1 - z) / (1 - PI) - (1 + lambda) / (lambda ** 2 + PI * (1 + lambda))
    G0 <- G0 * PI * (1 - PI) # PI derivative
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      #xx <- solve(weight[[x]]) # less computationally demanding
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
      PI <- lambda[, 2]
      lambda <- lambda[, 1]
      -sum(weight * (z * (log(PI) + log(lambda ** 2 + lambda + 1)) + (1 - z) * (log(1 - PI) + y * log(lambda) - (y - 1) * log(1 + lambda)) - log(lambda ** 2 + PI * (lambda + 1))))
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
      PI <- lambda[, 2]
      lambda <- lambda[, 1]
      G1 <- z * (2 * lambda + 1) / (lambda ** 2 + lambda + 1)
      G1 <- G1 + (1 - z) * (y / lambda - (y - 1) / (1 + lambda))
      G1 <- G1 - (2 * lambda + PI) / (lambda ** 2 + PI * (lambda + 1))
      G1 <- G1 * weight * lambda # lambda derivative
      G0 <- z / PI - (1 - z) / (1 - PI) - (1 + lambda) / (lambda ** 2 + PI * (1 + lambda))
      G0 <- G0 * weight * PI * (1 - PI) # PI derivative
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
      PI <- lambda[, 2]
      lambda <- lambda[, 1]
      Xlambda <- X[1:(nrow(X) / 2), 1:lambdaPredNumber]
      XPI <- X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
      res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
      
      G1 <- z * (2 * lambda + 1) / (lambda ** 2 + lambda + 1)
      G1 <- G1 + (1 - z) * (y / lambda - (y - 1) / (1 + lambda))
      G1 <- G1 - (2 * lambda + PI) / (lambda ** 2 + PI * (lambda + 1)) # lambda derivative
      G0 <- z / PI - (1 - z) / (1 - PI) - (1 + lambda) / (lambda ** 2 + PI * (1 + lambda)) # PI derivative
      
      # PI^2 derivative
      G00 <- -z / (PI ** 2) - (1 - z) / ((1 - PI) ** 2) + ((lambda + 1) ** 2) / ((lambda ** 2 + PI * (lambda + 1)) ** 2)
      G00 <- weight * (G0 * (PI * (1 - PI) * (1 - 2 * PI)) + G00 * ((PI * (1 - PI)) ** 2))  # second derivative of inverse logistic link
      G00 <- t(as.data.frame(XPI * G00)) %*% as.matrix(XPI)
      
      # mixed
      
      G01 <- lambda * (lambda + 2) / ((lambda ** 2 + PI * (lambda + 1)) ** 2)
      G01 <- G01 * lambda * PI * (1 - PI) * weight
      G01 <- t(as.data.frame(Xlambda * G01)) %*% as.matrix(XPI)
      
      # Beta^2 derivative
      G11 <- (1 - z) * ((y - 1) / ((1 + lambda) ** 2) - y / (lambda ** 2))
      G11 <- G11 + z * (2 / (lambda ** 2 + lambda + 1) - ((2 * lambda + 1) ** 2) / ((lambda ** 2 + lambda + 1) ** 2))
      G11 <- G11 + ((2 * lambda + PI) ** 2) / ((lambda ** 2 + PI * (lambda + 1)) ** 2) - 2 / (lambda ** 2 + PI * (lambda + 1))
      G11 <- (G11 * lambda * lambda + G1 * lambda) * weight # second derivative of log link
      G11 <- t(as.data.frame(Xlambda * G11)) %*% Xlambda
      res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
      res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
      res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- G01
      res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- t(G01)
      
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
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    N <- pw * (lambda ** 2 + lambda + 1) / (lambda ** 2 + PI * (lambda + 1))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    
    bigTheta1 <- -pw *  PI * (1 - PI) * (lambda ** 2 + lambda + 1) * (lambda + 1) / (lambda ** 2 + PI * (lambda + 1)) # w.r to PI
    bigTheta2 <- pw * lambda * ((2 * lambda + 1) * (lambda ** 2 + PI * (lambda + 1)) - (2 * lambda + PI) * (lambda ** 2 + lambda + 1)) / ((lambda ** 2 + PI * (lambda + 1)) ** 2) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * (lambda ** 2 + lambda + 1) * (1 - PI) * (1 + lambda) / ((lambda ** 2 + PI * (lambda + 1)) ** 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    ifelse(x == 1, PI * (lambda ** 2 + lambda + 1), (1 - PI) * (lambda ** x) / ((1 + lambda) ** (x - 1))) / (lambda ** 2 + PI * (lambda + 1))
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    CDF <- function(x) {
      p <- lambda / (1 + lambda)
      const <- -lambda * (p ** x + lambda * (p ** x - 1))
      polly <- lambda ** 2 + lambda + 1
      ifelse(x == Inf, 1, ifelse(x < 0, 0, ifelse(x < 1, (1 - PI) * (1 + lambda) / polly, (1 - PI) * (1 + lambda) / polly + PI + (1 - PI) * const / polly)))
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
      family = "Hurdleztgeom",
      parNum = 2,
      etaNames = c("lambda", "pi"),
      densityFunction = dFun,
      simulate = simulate
    ),
    class = "family"
  )
}

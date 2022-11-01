#' hurdle Zero truncated poisson model
#'
#' @return A object of class "family" containing objects \cr
#' makeMinusLogLike(y,X) - for creating negative likelihood function \cr
#' makeGradient(y,X) - for creating gradient function \cr
#' makeHessian(X) - for creating hessian \cr
#' linkfun - a link function to connect between linear predictor and model parameter in regression and a name of link function\cr
#' linkinv - an inverse function of link \cr
#' Dlink - a 1st derivative of link function \cr
#' mu.eta,Variance - Expected Value and Variance \cr
#' valedmu, valideta - for checking if regression arguments and valid\cr
#' family - family name\cr
#' Where: \cr
#' y is a vector of observed values \cr
#' X is a matrix / data frame of covariates
#' @importFrom lamW lambertW0
#' @export
Hurdleztpoisson <- function() {
  # Fist for lambda second for PI
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2]/ (1 - x[,2]))), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),1/(exp(-x[,2]) + 1)), ncol = 2, dimnames = dimnames(x))}
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
            "nontrunc" = PI + (1 - PI) * exp(-lambda) * (lambda * exp(lambda) - lambda) / (1 - lambda * exp(-lambda)),
            "trunc" = PI * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) + (1 - PI) * (lambda * exp(lambda) - 1) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    switch(type,
            "nontrunc" = PI + (1 - PI) * exp(-lambda) * lambda * (exp(lambda) * (1 + lambda) - 1) / (1 - lambda * exp(-lambda)),
            "trunc" = PI * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) + (1 - PI) * exp(-lambda) * lambda * (exp(lambda) * (1 + lambda) - 1) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    ) - (mu.eta(eta = eta, type = type) ** 2)
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- (1 - PI) * exp(-lambda) / (1 - lambda * exp(-lambda))
    
    G00 <- -1 + PI * (1 - PI) * exp(-2 * lambda) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) ** 2) - (1 - 2 * PI) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    G00 <- G00 * PI * (1 - PI) * prior
    
    # mixed
    
    G01 <- -PI * (1 - PI) * (-exp(-lambda) * (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) - exp(-2 * lambda) * (lambda - PI)) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) ** 2)
    G01 <- G01 * lambda * prior
    
    # Beta^2 derivative
    G11 <- -(PI * (exp(lambda) * (lambda + 1) - 1) + exp(lambda) * (-lambda ** 2 - 2) + exp(2 * lambda) + 1) / ((PI + exp(lambda) - lambda - 1) ** 2) - (exp(lambda) * (lambda ** 2 - lambda - exp(lambda) + 1) * z) / ((lambda - exp(lambda)) ** 2)
    G11 <- G11 * lambda * prior
    
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
    G1 <- z * lambda * (lambda - 1) / (exp(lambda) - lambda)
    G1 <- G1 + (1 - z) * (y - lambda)
    G1 <- G1 - lambda * exp(-lambda) * (lambda - PI) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    G0 <- (z - PI- PI * (1 - PI) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))) # PI derivative
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- solve(weight[[x]]) # less computationally demanding
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
      logistic <- -sum(weight * (z * log(PI) + (1 - z) * log(1 - PI)))
      rest <- -sum(z * log(1 - lambda * exp(-lambda)) +
      (1 - z) * (y * log(lambda) - lambda - log(factorial(y))) -
      log(1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)))
      rest + logistic
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
      G1 <- z * lambda * (lambda - 1) / (exp(lambda) - lambda)
      G1 <- G1 + (1 - z) * (y - lambda)
      G1 <- G1 - lambda * exp(-lambda) * (lambda - PI) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
      G1 <- G1 * weight
      G0 <- (z - PI- PI * (1 - PI) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))) # PI derivative
      G0 <- G0 * weight
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
      
      # PI^2 derivative
      G00 <- -1 + PI * (1 - PI) * exp(-2 * lambda) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) ** 2) - (1 - 2 * PI) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
      G00 <- G00 * PI * (1 - PI)
      G00 <- t(as.data.frame(XPI * G00 * weight)) %*% as.matrix(XPI)
      
      # mixed
      
      G01 <- -PI * (1 - PI) * (-exp(-lambda) * (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) - exp(-2 * lambda) * (lambda - PI)) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) ** 2)
      G01 <- G01 * lambda
      G01 <- t(as.data.frame(Xlambda * G01 * weight)) %*% as.matrix(XPI)
      
      # Beta^2 derivative
      G11 <- -(PI * (exp(lambda) * (lambda + 1) - 1) + exp(lambda) * (-lambda ** 2 - 2) + exp(2 * lambda) + 1) / ((PI + exp(lambda) - lambda - 1) ** 2) - (exp(lambda) * (lambda ** 2 - lambda - exp(lambda) + 1) * z) / ((lambda - exp(lambda)) ** 2)
      G11 <- G11 * lambda
      G11 <- t(as.data.frame(Xlambda * G11 * weight)) %*% Xlambda
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
    N <- pw * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    prob <- 1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)
    
    bigTheta1 <- -pw * (1 - lambda * exp(-lambda)) * exp(-lambda) / (prob ** 2) # w.r to PI
    bigTheta2 <- pw * (prob * (lambda * exp(-lambda) - exp(-lambda)) - ((1 - PI) * exp(-lambda) + lambda * exp(-lambda) - exp(-lambda))) / (prob ** 2) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    #f2 <- sum(pw * (1 - PI) * exp(-lambda) * (1 - lambda * exp(-lambda)) / (prob ** 2))
    f2 <- sum(pw * (1 - lambda * exp(-lambda)) * (1 - PI) * exp(-lambda) / (prob ** 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    ifelse(x == 1, PI * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)),
           (1 - PI) * (lambda ** x) * exp(-lambda) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) * factorial(x)))
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    CDF <- function(x) {
      ifelse(x == Inf, 1, ifelse(x < 0, 0, ifelse(x < 1, (1 - PI) * exp(-lambda) / (1 - lambda * exp(-lambda)), PI + (1 - PI) * (stats::ppois(x, lambda) - lambda * exp(-lambda)) / (1 - lambda * exp(-lambda)))))
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
      family = "Hurdleztpoisson",
      parNum = 2,
      etaNames = c("lambda", "pi"),
      densityFunction = dFun,
      simulate = simulate
    ),
    class = "family"
  )
}

#' One Inflated Zero truncated Poisson Model
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
oiztpoisson <- function() {
  # Fist for lambda second for omega
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2]/ (1 - x[,2]))), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),1/(exp(-x[,2]) + 1)), ncol = 2, dimnames = dimnames(x))}
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
            "nontrunc" = (omega + (1 - omega) * lambda / (exp(lambda) - 1) + (1 - omega) * lambda) * (1 - exp(-lambda)),
            "trunc" = omega + (1 - omega) * lambda / (exp(lambda) - 1) + (1 - omega) * lambda
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
            "nontrunc" = ((omega + (1 - omega) * lambda / (exp(lambda) - 1) + (1 - omega) * lambda + (lambda ** 2) * (1 - omega) / (1 - exp(-lambda))) * (1 - exp(-lambda)) - 
                            ((omega + (1 - omega) * lambda / (exp(lambda) - 1) + (1 - omega) * lambda) * (1 - exp(-lambda))) ** 2),
            "trunc" = ((omega + (1 - omega) * lambda / (exp(lambda) - 1) + (1 - omega) * lambda + (lambda ** 2) * (1 - omega) / (1 - exp(-lambda))) - 
                         ((omega + (1 - omega) * lambda / (exp(lambda) - 1) + (1 - omega) * lambda) ** 2))
    )
  }
  
  Wfun <- function(prior, y, eta, ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    theta <- exp(eta[, 2])
    temp <- exp(lambda) * omega + lambda * (1 - omega)
    temp1 <- exp(lambda) * omega + 1 - omega
    temp2 <- exp(lambda) + 1 - omega
    z <- temp / temp2 # expected for I's
    #z <- ifelse(y == 1, y, 0)
    G00 <- z * (-(1 - omega) * omega * ((exp(lambda) - lambda) ** 2) / ((omega * (exp(lambda) - lambda) + lambda) ** 2) +
                  (1 - omega) * (exp(lambda) - lambda) / (omega * (exp(lambda) - lambda) + lambda) - omega * (exp(lambda) - lambda) / (omega * (exp(lambda) - lambda) + lambda))
    G00 <- G00 - (1 - z)
    G00 <- G00 - ((1 - 2 * omega) * (exp(lambda) - 1 + omega) - omega * (1 - omega)) / ((exp(lambda) - 1 + omega) ** 2)
    G00 <- G00 * (1 - omega) * omega * prior
    
    # mixed derivative
    G01 <- lambda * exp(lambda) / ((exp(lambda) - 1 + omega) ** 2)
    G01 <- G01 + z * lambda * exp(lambda) * (lambda - 1) / ((omega * (exp(lambda) - lambda) + lambda) ** 2)
    G01 <- G01 * (1 - omega) * omega
    G01 <- G01 * prior
    
    G11 <- z * omega * exp(lambda) * (omega * (lambda + exp(lambda) - 1 - (lambda ** 2)) + lambda ** 2 - lambda + 1) / ((omega * (exp(lambda) - lambda) + lambda) ** 2)
    G11 <- G11 - exp(lambda) * (omega * lambda + omega + exp(lambda) - lambda - 1) / ((omega + exp(lambda) - 1) ** 2)
    G11 <- G11 * prior * lambda
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
    theta <- exp(eta[, 2])
    z <- ifelse(y == 1, y, 0)
    G1 <- z  * (exp(lambda) * omega + 1 - omega) / (exp(lambda) * omega + lambda - omega * lambda)
    G1 <- G1 + (1 - z) * y / lambda
    G1 <- G1 - exp(lambda) / (exp(lambda) - 1 + omega)
    G1 <- G1 * lambda # log link
    G0 <- z * (exp(lambda) - lambda) / (lambda + omega * (exp(lambda) - lambda))
    G0 <- G0 - (1 - z) / (1 - omega)
    G0 <- G0 - 1 / (exp(lambda) - 1 + omega)
    G0 <- G0 * omega * (1 - omega) # logit link
    
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
      -sum(weight * (z * (log(exp(lambda) * omega + (1 - omega) * lambda)) +
      (1 - z) * (log(1 - omega) + y * log(lambda) - log(factorial(y))) - 
      log(exp(lambda) - 1 + omega)))
    }
  }
  
  gradient <- function(y, X, weight = 1, ...) {
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
      G1 <- z  * (exp(lambda) * omega + 1 - omega) / (exp(lambda) * omega + lambda - omega * lambda)
      G1 <- G1 + (1 - z) * y / lambda
      G1 <- G1 - exp(lambda) / (exp(lambda) - 1 + omega)
      G1 <- G1 * weight * lambda # log link
      G0 <- z * (exp(lambda) - lambda) / (lambda + omega * (exp(lambda) - lambda))
      G0 <- G0 - (1 - z) / (1 - omega)
      G0 <- G0 - 1 / (exp(lambda) - 1 + omega)
      G0 <- G0 * weight * omega * (1 - omega) # logit link
      as.numeric(c(G1, G0) %*% X)
    }
  }
  
  hessian <- function(y, X, weight = 1, lambdaPredNumber, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    z <- as.numeric(y == 1)
    function (beta) {
      # TODO:
      # something wrong with derivatives (all)
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      lambda <- invlink(eta)
      omega <- lambda[, 2]
      lambda <- lambda[, 1]
      Xlambda <- X[1:(nrow(X) / 2), 1:lambdaPredNumber]
      Xomega <- X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
      res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
      
      temp <- exp(lambda) * omega + lambda * (1 - omega)
      temp1 <- exp(lambda) * omega + 1 - omega
      temp2 <- exp(lambda) + 1 - omega
      
      # omega^2 derivative
      G00 <- z * (-(1 - omega) * omega * ((exp(lambda) - lambda) ** 2) / ((omega * (exp(lambda) - lambda) + lambda) ** 2) +
      (1 - 2 * omega) * (exp(lambda) - lambda) / (omega * (exp(lambda) - lambda) + lambda))
      G00 <- G00 - (1 - z)
      G00 <- G00 - ((1 - 2 * omega) * (exp(lambda) - 1 + omega) - omega * (1 - omega)) / ((exp(lambda) - 1 + omega) ** 2)
      G00 <- G00 * (1 - omega) * omega
      G00 <- t(as.data.frame(Xomega * G00 * weight)) %*% as.matrix(Xomega)
      
      # mixed derivative
      G01 <- lambda * exp(lambda) / ((exp(lambda) - 1 + omega) ** 2)
      G01 <- G01 + z * lambda * exp(lambda) * (lambda - 1) / ((omega * (exp(lambda) - lambda) + lambda) ** 2)
      G01 <- G01 * (1 - omega) * omega
      G01 <- t(as.data.frame(Xlambda) * G01 * weight) %*% as.matrix(Xomega)
      
      # Beta^2 derivative
      G11 <- z * omega * exp(lambda) * (omega * (lambda + exp(lambda) - 1 - lambda ** 2) + lambda ** 2 - lambda + 1) / ((omega * (exp(lambda) - lambda) + lambda) ** 2)
      G11 <- G11 - exp(lambda) * (omega * lambda + omega + exp(lambda) - lambda - 1) / ((omega + exp(lambda) - 1) ** 2)
      G11 <- G11 * weight * lambda
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
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    N <- pw / (1 - (1 - omega) * exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    
    coefficient <- (1 - (1 - omega) * exp(-lambda))
    bigTheta1 <- -pw * exp(-lambda) / (coefficient ** 2)# w.r to omega
    bigTheta1 <- bigTheta1 * omega * (1 - omega)
    bigTheta2 <- -pw * lambda * exp(-lambda) * (1 - omega) / (coefficient ** 2) # w.r to lambda
    bigTheta2 <- bigTheta2 * lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * (1 - omega) * exp(-lambda) / (coefficient ** 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    ifelse(x == 1, exp(lambda) * omega + (1 - omega) * lambda,
           (1 - omega) * (lambda ** x) / factorial(x)) / (exp(lambda) - 1 + omega)
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    CDF <- function(x) {
      ifelse(x == Inf, 1, ifelse(x < 0, 0, ifelse(x < 1, (1 - omega) * exp(-lambda), omega +  (1 - omega) * stats::ppois(x, lambda))))
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
      family = "oiztpoisson",
      parNum = 2,
      etaNames = c("lambda", "omega"),
      densityFunction = dFun,
      simulate = simulate
    ),
    class = "family"
  )
}

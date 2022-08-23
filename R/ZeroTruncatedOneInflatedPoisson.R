#' Zero truncated One Inflated Poisson Model
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
ztoipoisson <- function() {
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
    z <- omega + (1 - omega) * lambda / (exp(lambda) - 1) # expected for I's
    #z <- ifelse(y == 1, y, 0)
    eml <- exp(-lambda)
    term <- lambda / (exp(lambda) - 1)
    G00 <- -omega * (1 - omega) + z * (theta / (theta + term) - (theta / (theta + term)) ** 2)
    G00 <- G00 * prior
    
    # mixed derivative
    G01 <- -z * theta / ((theta + term) ** 2)
    G01 <- G01 * (term - (term ** 2) * exp(lambda)) * prior
    
    G11 <- (1 - z) * (lambda * eml / ((1 - eml) ** 2) - 1 / (1 - eml))
    wwd <- 1 / (exp(lambda) - 1)
    G11 <- G11 + z * ((wwd - wwd * term * exp(lambda) - 2 * lambda * wwd * (wwd - wwd * term * exp(lambda))) *
                        (theta + term) - (term - term ** 2) * (wwd - wwd * term * exp(lambda))) / (term + theta)
    G11 <- lambda * prior * G11
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
    mm <- exp(lambda) - 1
    mm <- 1 / mm
    G1 <- z * (mm - mm * mm * lambda * exp(lambda)) / (theta + lambda * mm)
    G1 <- G1 + (1 - z) * (y / lambda - mm * exp(lambda))
    G1 <- G1 * lambda
    G0 <- -omega + z * (theta / (theta + lambda / (exp(lambda) - 1)))
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- solve(weight[[x]])
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
      -sum(weight * (-log(1 + theta) + z * log(theta + lambda / (exp(lambda) - 1)) +
                    (1 - z) * (y * log(lambda) - log(exp(lambda) - 1) - log(factorial(y)))))
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
      theta <- exp(eta[, 2])
      mm <- exp(lambda) - 1
      mm <- 1 / mm
      G1 <- z * (mm - mm * mm * lambda * exp(lambda)) / (theta + lambda * mm)
      G1 <- G1 + (1 - z) * (y / lambda - mm * exp(lambda))
      G1 <- G1 * weight * lambda
      G0 <- -omega + z * (theta / (theta + lambda / (exp(lambda) - 1))) # omega derivative
      G0 <- G0 * weight
      as.numeric(c(G1, G0) %*% X)
    }
  }
  
  hessian <- function(y, X, weight = 1, lambdaPredNumber, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    z <- as.numeric(y == 1)
    function (beta) {
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      lambda <- invlink(eta)
      omega <- lambda[, 2]
      lambda <- lambda[, 1]
      theta <- exp(eta[, 2])
      Xlambda <- X[1:(nrow(X) / 2), 1:lambdaPredNumber]
      Xomega <- X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
      res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
      eml <- exp(-lambda)
      coefficient <- (1 / (1 - eml) - lambda * eml / ((1 - eml) ** 2))
      
      # omega^2 derivative
      term <- lambda / (exp(lambda) - 1)
      G00 <- t(as.data.frame(Xomega * (-omega * (1 - omega) + z * (theta / (theta + term) - (theta / (theta + term)) ** 2)) * weight)) %*% as.matrix(Xomega)
      
      # mixed derivative
      xx <- -z * theta / ((theta + term) ** 2)
      xx <- xx * (term - (term ** 2) * exp(lambda))
      G01 <- t(as.data.frame(Xlambda) * xx * weight) %*% as.matrix(Xomega)
      
      # Beta^2 derivative
      G11 <- (1 - z) * (lambda * eml / ((1 - eml) ** 2) - 1 / (1 - eml))
      wwd <- 1 / (exp(lambda) - 1)
      G11 <- G11 + z * ((wwd - wwd * term * exp(lambda) - 2 * lambda * wwd * (wwd - wwd * term * exp(lambda))) *
                        (theta + term) - (term - term ** 2) * (wwd - wwd * term * exp(lambda))) / (term + theta)
      G11 <- lambda * weight * G11
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
    NULL
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
    lambda <- lambda[, 1]
    N <- pw / (1 - exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    theta <- exp(eta[, 2])
    ml <- (1 - exp(-lambda)) ** 2
    
    bigTheta1 <- rep(0, nrow(eta)) # w.r to omega
    bigTheta2 <- pw * (exp(log(lambda) - lambda) / ml) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * exp(-lambda) / ml)
    
    f1 + f2
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
      family = "ztoipoisson",
      parNum = 2,
      etaNames = c("lambda", "omega")
    ),
    class = "family"
  )
}

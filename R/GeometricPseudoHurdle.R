#' Zero truncated hurdle geometric model
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
#' @export
ztHurdlegeom <- function() {
  # Fist for lambda second for PI
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2]/ (1 - x[,2]))), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),1/(exp(-x[,2]) + 1)), ncol = 2, dimnames = dimnames(x))}
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
            "nontrunc" = (1 - exp(-lambda)) * (PI + lambda - PI * lambda),
            "trunc" = PI + (1 - PI) * (2 + lambda)
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
            "nontrunc" = PI * (1 - exp(-lambda)) + (1 - PI) * (lambda ** 2 + lambda - lambda * exp(-lambda)),
            "trunc" = (PI + (1 - PI) * (lambda + lambda ** 2 - lambda * exp(-lambda)) / (1 - exp(-lambda))) - (PI + (1 - PI) * lambda) ** 2
    )
  }
  
  Wfun <- function(prior, y, eta, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- PI
    Ey <- mu.eta(eta = eta)
    #z <- ifelse(y == 1, y, 0)
    term <- -(PI * (1 - PI))
    G00 <- term * prior
    S <- 1 / (1 + lambda)
    term <- (1 - z) * lambda * (Ey - 1) * (S ** 2)
    G11 <- -term * prior
    G01 <- rep(0, nrow (eta))
    
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
    S <- 1 / (1 + lambda)
    G1 <- ifelse(z, 0, ((y - 1) * S - 1)  * weight)
    G0 <- (z - PI) # PI derivative
    
    pseudoResid <- matrix(c(G1 / weight[, 1], G0 / weight[, 4]), ncol = 2)
  
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
      zot <- -sum(ifelse(z, 0, weight * ((y - 2) * log(lambda) - (y - 1) * log(1 + lambda))))
      zot + logistic
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
      PI <- lambda[, 2]
      lambda <- lambda[, 1]
      S <- 1 / (1 + lambda)
      G1 <- weight * ifelse(z, 0, ((y - 1) * S - 1))
      G0 <- (z - PI) # PI derivative
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
      PI <- lambda[, 2]
      lambda <- lambda[, 1]
      Xlambda <- X[1:(nrow(X) / 2), 1:lambdaPredNumber]
      XPI <- X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
      res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
      
      # PI^2 derivative
      term <- -(PI * (1 - PI))
      G00 <- t(as.data.frame(XPI * term * weight)) %*% as.matrix(XPI)
      
      # Beta^2 derivative
      S <- 1 / (1 + lambda)
      term <- ifelse(z, 0, lambda * (y - 1) * (S ** 2))
      G11 <- -term * weight
      G11 <- t(as.data.frame(Xlambda * G11)) %*% Xlambda
      res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
      res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
      res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- 0
      res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- 0
      
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
    S <- 1 / (1 + lambda)
    prob <- 1 - S - lambda * (S ** 2)
    N <- pw * (1 - lambda * (S ** 2)) / prob
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    M <- 1 + lambda
    S <- 1 / M
    prob <- 1 - S - lambda * (S ** 2)
    
    bigTheta1 <- rep(0, nrow(eta)) # w.r to PI
    bigTheta2 <- pw * as.numeric(lambda * (prob * (lambda - 1) * (S ** 3) - 
    2 * lambda * (S ** 3) * (1 - lambda * (S ** 2))) / (prob ** 2)) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * (exp(-lambda) / ((1 - exp(-lambda) - lambda * exp(-lambda)) ** 2)))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    ifelse(x == 1, PI, (1 - PI) * (lambda ** (x - 2)) / ((1 + lambda) ** (x - 1)))
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
      family = "ztHurdlegeom",
      parNum = 2,
      etaNames = c("lambda", "pi"),
      densityFunction = dFun
    ),
    class = "family"
  )
}

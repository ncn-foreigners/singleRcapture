#' Zero truncated hurdle poisson model
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
ztHurdlepoisson <- function() {
  # Fist for lambda second for PI
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2]/ (1 - x[,2]))), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),1/(exp(-x[,2]) + 1)), ncol = 2, dimnames = dimnames(x))}
  
  mu.eta <- function(eta, type = "trunc", ...) {
    # TODO
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
            "nontrunc" = (1 - exp(-lambda)) * (PI + lambda - PI * lambda),
            "trunc" = PI + (1 - PI) * (lambda - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda))
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    # TODO
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
    #z <- ifelse(y == 1, y, 0)
    term <- -(PI * (1 - PI))
    G00 <- term * prior
    term <- (1 - z) * ((2 + lambda ** 2) * exp(lambda) - exp(2 * lambda) - 1) / ((exp(lambda) - lambda - 1) ** 2)
    G11 <- lambda * term * prior
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
    G1 <- ifelse(z, 0 , (y - lambda - lambda * lambda / (exp(lambda) - lambda - 1)))
    G0 <- (z - PI) # PI derivative
    
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
      PI <- lambda[, 2]
      lambda <- lambda[, 1]
      logistic <- -sum(weight * (z * log(PI) + (1 - z) * log(1 - PI)))
      zot <- -sum(ifelse(z, 0, weight * (y * log(lambda) - lambda - log(factorial(y)) - log(1 - exp(-lambda) - lambda * exp(-lambda)))))
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
      G1 <- weight * ifelse(z, 0 , (y - lambda - lambda * lambda / (exp(lambda) - lambda - 1)))
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
      term <- ifelse(z, 0, ((2 + lambda ** 2) * exp(lambda) - exp(2 * lambda) - 1) / ((exp(lambda) - lambda - 1) ** 2))
      G11 <- lambda * term * weight
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
    NULL
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
    lambda <- lambda[, 1]
    N <- pw * (1 - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    
    bigTheta1 <- rep(0, nrow(eta)) # w.r to PI
    bigTheta2 <-as.numeric(pw * lambda * (1 - exp(lambda)) / ((1 + lambda - exp(lambda)) ** 2)) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * ((1 - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda))))
    
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
      family = "ztHurdlepoisson",
      parNum = 2,
      etaNames = c("lambda", "pi")
    ),
    class = "family"
  )
}

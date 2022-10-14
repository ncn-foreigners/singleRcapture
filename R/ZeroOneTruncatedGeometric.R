#' Zero-one truncated geometric model
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
zotgeom <- function() {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    switch (type,
      "nontrunc" = lambda,
      "trunc" = 2 + lambda
    )
  }
  
  variance <- function(eta, disp, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    switch (type,
      "nontrunc" = lambda * (lambda - 1),
      "trunc" = lambda * (lambda + 1)
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- exp(eta)
    lambda / (1 + lambda)
  }
  
  funcZ <- function(eta, weight, y, mu, ...) {
    ((y - 1) / (1 + exp(eta)) - 1) / weight
  }
  
  minusLogLike <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)
    
    function(beta) {
      eta <- as.matrix(X) %*% beta
      lambda <- exp(eta)
      
      -sum(weight * ((y - 2) * log(lambda) - (y - 1) * log(1 + lambda)))
    }
  }
  
  
  gradient <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)
    
    function(beta) {
      eta <- X %*% beta
      lambda <- exp(eta)
      S <- 1 / (1 + lambda)
      
      # Beta derivative
      if (NbyK) {
        return(((y - 1) * S - 1)  * weight * as.data.frame(X))
      }
      if (vectorDer) {
        return(matrix(((y - 1) * S - 1)  * weight, ncol = 1))
      }
      G1 <- t(((y - 1) * S - 1)  * weight) %*% X
      
      G1
    }
  }
  
  hessian <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)
    
    function(beta) {
      eta <- X %*% beta
      lambda <- exp(eta)
      S <- 1 / (1 + lambda)
      
      # second beta derivative
      
      G11 <- -t(as.data.frame(X) * lambda * (y - 1) * (S ** 2) * weight) %*% X
      
      G11
    }
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  dev.resids <- function (y, eta, wt, ...) {
    mu <- invlink(eta)
    mu1 <- mu.eta(eta = eta)
    loghm1y <- ifelse(y > 2, log(y - 2), 0)
    sign(y - mu1) * sqrt(-2 * wt * ((y - 2) * eta - (y - 1) * log(1 + mu1) - (y - 2) * loghm1y + (y - 1) * log(y - 1)))
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
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
    M <- 1 + lambda
    S <- 1 / M
    prob <- 1 - S - lambda * (S ** 2)
    
    bigTheta <- t(Xvlm) %*% (pw * as.numeric(lambda *
                            (prob * (lambda - 1) * (S ** 3) -
                            2 * lambda * (S ** 3) *
                            (1 - lambda * (S ** 2))) / (prob ** 2)))
    
    bigTheta <- as.vector(bigTheta)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * (1 - lambda * S * S) / prob)
    
    f1 + f2
  }
  
simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- invlink(eta)
    lb <- stats::pnbinom(lower, mu=lambda, size = 1)
    ub <- stats::pnbinom(upper, mu=lambda, size = 1)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qnbinom(p_u, mu=lambda, size = 1)
    sims
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    stats::dgeom(x = x, prob = (1 / (1 + lambda))) / (1 - stats::dgeom(x = 0, prob = (1 / (1 + lambda))) - stats::dgeom(x = 1, prob = (1 / (1 + lambda))))
  }
  
  structure(
    list(
      makeMinusLogLike = minusLogLike,
      makeGradient = gradient,
      makeHessian = hessian,
      linkfun = link,
      linkinv = invlink,
      dlink = dlink,
      mu.eta = mu.eta,
      link = "log",
      valideta = function (eta) {TRUE},
      variance = variance,
      Wfun = Wfun,
      funcZ = funcZ,
      dev.resids = dev.resids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      simulate=simulate,
      family = "zotgeom",
      parNum = 1,
      etaNames = "lambda",
      densityFunction = dFun
    ),
    class = "family"
  )
}
#' Zero truncated geometric model
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
ztgeom <- function() {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    switch (type,
      "nontrunc" = lambda,
      "trunc" = 1 + lambda
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    mu <- mu.eta(eta)
    switch (type,
      nontrunc = mu ** 2 - mu - 1 / mu + 2,
      trunc = (mu + 1) / mu
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- exp(eta)
    lambda * (1 + lambda) / ((1 + lambda) ** 2)
  }
  
  funcZ <- function(eta, weight, y, ...) {
    mu <- 1 + exp(eta)
    (y  / mu - 1) / weight
  }
  
  minusLogLike <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)
    
    function(arg) {
      beta <- arg
      eta <- as.matrix(X) %*% beta
      lambda <- exp(eta)
      
      -sum(weight * ((y - 1) * log(lambda) - y * log(1 + lambda)))
    }
  }
  
  
  gradient <- function(y, X, weight = 1, NbyK, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)
    
    function(arg) {
      beta <- arg
      eta <- X %*% beta
      lambda <- exp(eta)
      S <- 1 / (1 + lambda)
      
      # Beta derivative
      if (NbyK) {
        return(as.data.frame(X) * (y * S - 1)  * weight)
      }
      G1 <- t((y * S - 1)  * weight) %*% X
      
      G1
    }
  }
  
  hessian <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)
    
    function(arg) {
      beta <- arg
      eta <- X %*% beta
      lambda <- exp(eta)
      S <- 1 / (1 + lambda)
      
      # second beta derivative
      
      G11 <- -t(as.data.frame(X) * lambda * y * (S ** 2) * weight) %*% X
      
      G11
    }
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  dev.resids <- function (y, eta, wt, ...) {
    mu <- invlink(eta)
    mu1 <- mu.eta(eta = eta)
    hm1y <- y - 1 # thats an analytic inverse for geometric
    #log1mexphm1y <- ifelse(y > 1, log(1 - exp(-hm1y)), 0)
    loghm1y <- ifelse(y > 1, log(hm1y), 0)
    sign(y - mu1) * sqrt(-2 * wt * ((y - 1) * eta - y * log(mu1) - (y - 1) * loghm1y + y * log(y)))
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
    pr <- 1 - 1 / (1 + lambda)
    N <- pw / pr
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    pr <- 1 - 1 / (1 + lambda)
    
    bigTheta <- -(pw * as.numeric(lambda / 
                 ((1 - (1 + lambda)) ** 2))) %*% as.matrix(Xvlm)
    bigTheta <- as.vector(bigTheta)
    
    f1 <- t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <- sum(pw * (1 - pr) / (pr ** 2))
    
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
    stats::dgeom(x = x, prob = (1 / (1 + lambda))) / (1 - stats::dgeom(x = 0, prob = (1 / (1 + lambda))))
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
      simulate = simulate,
      family = "ztgeom",
      parNum = 1,
      etaNames = "lambda",
      densityFunction = dFun
    ),
    class = "family"
  )
}
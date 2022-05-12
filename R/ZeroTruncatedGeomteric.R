#' Zero truncated geometric model
#'
#' @return A object of class "family" containing objects \cr
#' make_minusloglike(y,X) - for creating negative likelihood function \cr
#' make_gradient(y,X) - for creating gradient function \cr
#' make_hessian(X) - for creating hessian \cr
#' linkfun - a link function to connect between linear predictor and model parameter in regression and a name of link function\cr
#' linkinv - an inverse function of link \cr
#' Dlink - a 1st derivative of link function \cr
#' mu.eta,Variance - Expected Value and Variance \cr
#' aic - for aic computation\cr
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
  
  mu.eta <- function(disp = NULL, eta) {
    (1 + exp(eta))
  }
  
  variance <- function(disp = NULL, mu) {
    #(((mu - 1) ** 3) / mu + 2 * mu - 1)
    (mu ** 2 - mu - 1 / mu + 2)
  }
  
  minusLogLike <- function(y, X, weight = 1) {
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
  
  
  gradient <- function(y, X, weight = 1) {
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
      G1 <- t((y * S - 1)  * weight) %*% X
      
      G1
    }
  }
  
  hessian <- function(y, X, weight = 1) {
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
  
  dev.resids <- function (y, mu, wt, disp = NULL) {
    NULL
  }
  
  aic <- function(y, mu, wt, dev) {
    -2 * sum(wt * ((y - 1) * log(mu) - y * log(1 + mu)))
  }
  
  pointEst <- function (disp, pw, lambda, contr = FALSE) {
    pr <- 1 - 1 / (1 + lambda)
    N <- pw / pr
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (beta, pw, lambda, disp, hess, X) {
    pr <- 1 - 1 / (1 + lambda)
    I <- as.matrix(-hess(beta))
    
    bigTheta <- -(pw * as.numeric(lambda / 
                 ((1 - (1 + lambda)) ** 2))) %*% as.matrix(X)
    bigTheta <- as.vector(bigTheta)
    
    f1 <- t(bigTheta) %*% solve(I) %*% bigTheta
    f2 <- sum(pw * (1 - pr) / (pr ** 2))
    
    f1 + f2
  }
  
  
  R <- list(make_minusloglike = minusLogLike,
            make_gradient = gradient,
            make_hessian = hessian,
            linkfun = link,
            linkinv = invlink,
            dlink = dlink,
            mu.eta = mu.eta,
            aic = aic,
            link = "log",
            valideta = function (eta) {TRUE},
            variance = variance,
            dev.resids = dev.resids,
            validmu = validmu,
            pointEst = pointEst,
            popVar= popVar,
            family = "ztgeom")
  class(R) <- "family"
  R
}
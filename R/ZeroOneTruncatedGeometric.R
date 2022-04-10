#' Zero One Truncated geometric model
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
zotgeom <- function() {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }
  
  mu.eta <- function(eta, disp) {
    lambda <- invlink(eta)
    Pr <- 1 / (1 + lambda) + lambda / ((1 + lambda) ** 2)
    G <- (lambda - lambda * ((1 + lambda) ** (-2))) / (1 - Pr)
    G
  }
  
  variance <- function(mu, disp) {
    mu * (1 + mu)
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
      
      -sum(weight * ((y - 2) * log(lambda) - (y - 1) * log(1 + lambda)))
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
      G1 <- t(((y - 1) * S - 1)  * weight) %*% X
      
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
      
      G11 <- -t(as.data.frame(X) * lambda * (y - 1) * (S ** 2) * weight) %*% X
      
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
    -2 * sum(wt * ((y - 2) * log(mu) - (y - 1) * log(1 + mu)))
  }
  
  pointEst <- function (disp, pw, lambda) {
    z <- 1
    S <- 1 / (1 + lambda)
    prob <- 1 - S - lambda * (S ** 2)
    N <- sum(pw * (1 - lambda * (S ** 2)) / prob)
    N
  }
  
  popVar <- function (beta, pw, lambda, disp, hess, X) {
    alpha <- 1
    z <- 1 / alpha
    M <- 1 + lambda
    S <- 1 / M
    prob <- 1 - S - lambda * (S ** 2)
    I <- as.matrix(-hess(beta))
    
    bigTheta <- t(as.matrix(X)) %*% (pw * as.numeric(lambda *
                                    (prob * (lambda - 1) * (S ** 3) -
                                    2 * lambda * (S ** 3) *
                                    (1 - lambda * (S ** 2))) / (prob ** 2)))
    
    bigTheta <- as.vector(bigTheta)
    
    f1 <-  t(bigTheta) %*% solve(I) %*% bigTheta
    f2 <-  sum(pw * ((1 - lambda * (S ** 2)) ** 2) *
               (1 - prob) / (prob ** 2))
    
    variation <- f1 + f2
    variation
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
            family = "zotgeom")
  class(R) <- "family"
  R
}
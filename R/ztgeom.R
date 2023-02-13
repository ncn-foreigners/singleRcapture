#' @rdname singleRmodels
#' @export
ztgeom <- function(...) {
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
  
  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it comfort to how swith in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(arg) {
        beta <- arg
        eta <- as.matrix(X) %*% beta
        lambda <- exp(eta)
        
        -sum(weight * ((y - 1) * log(lambda) - y * log(1 + lambda)))
      },
      function(beta) {
        eta <- X %*% beta
        lambda <- exp(eta)
        S <- 1 / (1 + lambda)
        
        # Beta derivative
        if (NbyK) {
          return(as.data.frame(X) * (y * S - 1)  * weight)
        }
        if (vectorDer) {
          return(matrix((y * S - 1) * weight, ncol = 1))
        }
        t((y * S - 1)  * weight) %*% X
      },
      function(beta) {
        eta <- X %*% beta
        lambda <- exp(eta)
        S <- 1 / (1 + lambda)
        
        # second beta derivative
        
        -t(as.data.frame(X) * lambda * y * (S ** 2) * weight) %*% X
      },
    )
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  dev.resids <- function (y, eta, wt, ...) {
    mu <- invlink(eta)
    mu1 <- mu.eta(eta = eta)
    hm1y <- y - 1 # that's an analytic inverse for geometric
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
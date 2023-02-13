#' @rdname singleRmodels
#' @export
zotpoisson <- function(...) {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }

  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = (lambda - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda))
    )
  }

  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    switch (type,
            "nontrunc" = lambda,
            "trunc" = (lambda - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda))
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- invlink(eta)
    matrix(-lambda * (((2 + lambda ** 2) * exp(lambda) - exp(2 * lambda) - 1) /
                     ((exp(lambda) - lambda - 1) ** 2)), 
           ncol = 1, dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- invlink(eta)
    (y - lambda - lambda * lambda / (exp(lambda) - lambda - 1)) / weight
  }

  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it comfort to how swith in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- exp(eta)
        -sum(weight * (y * eta - lambda - log(factorial(y)) -
        log(1 - exp(-lambda) - lambda * exp(-lambda))))
      },
      function(beta) {
        lambda <- exp(as.matrix(X) %*% beta)
        if (NbyK) {
          return(as.data.frame(X) * weight * (y - lambda - lambda * lambda / (exp(lambda) - lambda - 1)))
        }
        if (vectorDer) {
          return(matrix(weight * (y - lambda - lambda * lambda / (exp(lambda) - lambda - 1)), ncol = 1))
        }
        
        t(as.matrix(X)) %*% (weight * (y - lambda - lambda * lambda / (exp(lambda) - lambda - 1)))
      },
      function(beta) {
        lambda <- exp(as.matrix(X) %*% beta)
        
        term <- ((2 + lambda ** 2) * exp(lambda) - exp(2 * lambda) - 1) / ((exp(lambda) - lambda - 1) ** 2)
        
        t(X) %*% as.matrix(t(t(as.data.frame(X) * lambda * term * weight)))
      }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  devResids <- function(y, eta, wt, ...) {
    lambda <- invlink(eta)
    mu <- mu.eta(eta = eta)
    
    inverseFunction <- function(y) {stats::uniroot(
      f = function(x) {mu.eta(x) - y}, 
      lower = -log(y), upper = y * 10, 
      tol = .Machine$double.eps
    )$root}
    
    etaSat <- vector("numeric", length = length(y))
    
    etaSat[y == 2] <- -Inf
    etaSat[y > 2] <- sapply(y[y > 2], FUN = inverseFunction)
    lambdaSat <- exp(etaSat)
    print(lambdaSat)
    
    lFit <- y * eta - lambda - log(1 - exp(-lambda) - lambda * exp(-lambda))
    lSat <- ifelse(y == 2, log(2), # log(2) is the limit as lambda->0^+
    y * etaSat - lambdaSat - log(1 - exp(-lambdaSat) - lambdaSat * exp(-lambdaSat)))
    
    sign(y - mu) * sqrt(-2 * wt * (lFit - lSat))
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
    N <- (pw * (1 - lambda * exp(-lambda)) /
         (1 - exp(-lambda) - lambda * exp(-lambda)))
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    Xvlm <- as.data.frame(Xvlm)
    prob <- (1 - exp(-lambda) - lambda * exp(-lambda))
    term <- (1 - lambda * exp(-lambda)) ** 2
    
    f1 <- t(Xvlm) %*% (as.numeric(pw * lambda * (1 - exp(lambda)) /
                                 ((1 + lambda - exp(lambda)) ** 2)))
    
    f1 <- t(f1) %*% as.matrix(cov) %*% f1
    
    f2 <- sum(pw * term * (1 - prob) / (prob ** 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    stats::dpois(x = x, lambda = lambda) / (1 - stats::dpois(x = 0, lambda = lambda) - stats::dpois(x = 1, lambda = lambda))
  }

  simulate <- function(n, eta, lower = 1, upper = Inf) {
    lambda <- invlink(eta)
    lb <- stats::ppois(lower, lambda)
    ub <- stats::ppois(upper, lambda)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qpois(p_u, lambda)
    sims
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
      devResids = devResids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      simulate = simulate,
      family = "zotpoisson",
      parNum = 1,
      etaNames = "lambda",
      densityFunction = dFun
    ),
    class = "family"
  )
}

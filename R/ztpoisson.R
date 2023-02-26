#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @importFrom stats glm.fit
#' @importFrom stats poisson
#' @export
ztpoisson <- function(...) {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }

  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = lambda / (1 - exp(-lambda))
    )
  }

  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = mu.eta(eta = eta) * (1 + lambda - mu.eta(eta = eta))
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- invlink(eta)
    matrix(-lambda * ((exp(-lambda) + lambda * exp(- lambda) - 1) / ((1 - exp(-lambda)) ^ 2)), 
           ncol = 1, dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, ...) {
    (y - mu.eta(eta)) / weight
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
        -sum(weight * (y * eta - log(exp(exp(eta)) - 1) - log(factorial(y))))
      },
      function(beta) {
        lambda <- exp(as.matrix(X) %*% beta)
        if (NbyK) {
          return(as.data.frame(X) * weight * (y - lambda / (1 - exp(-lambda))))
        }
        if (vectorDer) {
          return(matrix(weight * (y - lambda / (1 - exp(-lambda))), ncol = 1))
        }
        t(as.matrix(X)) %*% (weight * (y - lambda / (1 - exp(-lambda))))
      },
      function(beta) {
        lambda <- exp(as.matrix(X) %*% beta)
        eml <- exp(-lambda)
        coefficient <- 1 / (1 - eml) - lambda * eml / ((1 - eml) ^ 2)
        
        -(t(as.matrix(X) * as.numeric(weight * (1 / (1 - eml) - lambda * eml / ((1 - eml) ^ 2)) * lambda))) %*% as.matrix(X)
      }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  devResids <- function(y, eta, wt, ...) {
    lambda <- invlink(eta)
    #hm1y <- ifelse(y > 1, VGAM::lambertW(-y * exp(-y)) + y, 0)
    hm1y <- ifelse(y > 1, lamW::lambertW0(-y * exp(-y)) + y, 0)
    sign(y - mu.eta(eta = eta)) * sqrt(-2 * wt * (y * eta - lambda - log(1 - exp(-lambda)) - y * ifelse(y > 1, log(hm1y), 0) + hm1y + ifelse(y > 1, log(1 - exp(-hm1y)), 0)))
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
    N <- pw / (1 - exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    Xvlm <- as.data.frame(Xvlm)
    lambda <- invlink(eta)

    f1 <- colSums(-Xvlm * pw * (exp(log(lambda) - lambda) / ((1 - exp(-lambda)) ^ 2)))
    f1 <- t(f1) %*% as.matrix(cov) %*% f1

    f2 <- sum(pw * exp(-lambda) / ((1 - exp(-lambda)) ^ 2))

    f1 + f2
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- invlink(eta)
    lb <- stats::ppois(lower, lambda)
    ub <- stats::ppois(upper, lambda)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qpois(p_u, lambda)
    sims
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    stats::dpois(x = x, lambda = lambda) / (1 - stats::dpois(x = 0, lambda = lambda))
  }
  
  getStart <- expression(
    start <- stats::glm.fit(
      x = variables[wch$reg, ],
      y = observed[wch$reg],
      family = stats::poisson(),
      weights = priorWeights[wch$reg],
      ...
    )$coefficients
  )
  
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
      family = "ztpoisson",
      parNum = 1,
      etaNames = "lambda",
      densityFunction = dFun,
      getStart = getStart
    ),
    class = "family"
  )
}

#' @rdname singleRmodels
#' @importFrom lamW lambertW0
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
    matrix(-lambda * ((exp(-lambda) + lambda * exp(- lambda) - 1) / ((1 - exp(-lambda)) ** 2)), 
           ncol = 1, dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, ...) {
    (y - mu.eta(eta)) / weight
  }

  minusLogLike <- function(y, X, weight = 1, ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    function(beta) {
      eta <- as.matrix(X) %*% beta
      lambda <- exp(eta)
      -sum(weight * (y * eta - log(exp(lambda) - 1) - log(factorial(y))))
    }
  }

  gradient <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      lambda <- exp(as.matrix(X) %*% beta)
      mu <- lambda / (1 - exp(-lambda))
      if (NbyK) {
        return(as.data.frame(X) * weight * (y - mu))
      }
      if (vectorDer) {
        return(matrix(weight * (y - mu), ncol = 1))
      }
      t(as.matrix(X)) %*% (weight * (y - mu))
    }
  }

  hessian <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      lambda <- exp(as.matrix(X) %*% beta)
      eml <- exp(-lambda)
      coefficient <- 1 / (1 - eml) - lambda * eml / ((1 - eml) ** 2)

      dmu <- weight * as.numeric(coefficient) # This was probably the dumbest mistake I've made in last 8 months
      dlam <- as.matrix(X * as.numeric(lambda))

      -((t(as.matrix(X) * dmu)) %*% dlam)
    }
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  dev.resids <- function(y, eta, wt, ...) {
    mu <- invlink(eta)
    mu1 <- mu.eta(eta = eta)
    #hm1y <- ifelse(y > 1, VGAM::lambertW(-y * exp(-y)) + y, 0)
    hm1y <- ifelse(y > 1, lamW::lambertW0(-y * exp(-y)) + y, 0)
    log1mexphm1y <- ifelse(y > 1, log(1 - exp(-hm1y)), 0)
    loghm1y <- ifelse(y > 1, log(hm1y), 0)
    sign(y - mu1) * sqrt(-2 * wt * (y * eta - mu - log(1 - exp(-mu)) - y * loghm1y + hm1y + log1mexphm1y))
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
    ml <- (1 - exp(-lambda)) ** 2

    f1 <- colSums(-Xvlm * pw * (exp(log(lambda) - lambda) / ml))
    f1 <- t(f1) %*% as.matrix(cov) %*% f1

    f2 <- sum(pw * exp(-lambda) / ml)

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
      family = "ztpoisson",
      parNum = 1,
      etaNames = "lambda",
      densityFunction = dFun
    ),
    class = "family"
  )
}

#' @rdname singleRmodels
#' @importFrom stats glm.fit
#' @importFrom stats poisson
#' @export
chao <- function(...) {
  link <- function(x) {log(x / 2)}
  invlink <- function (x) {2 * exp(x)}
  dlink <- function(lambda) {1 / lambda}
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = 1 / (1 + exp(-eta))
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = (1 / (1 + exp(-eta))) * (1 / (1 + exp(eta)))
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- invlink(eta)
    L1 <- lambda / 2
    (L1 / ((1 + L1) ^ 2))
  }
  
  funcZ <- function(eta, weight, y, mu, ...) {
    lambda <- invlink(eta)
    L1 <- lambda / 2
    (L1 * (y - 1) + y) / (L1 + 1) / weight
  }

  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    y <- as.numeric(y)
    z <- y - 1
    if (is.null(weight)) {
      weight <- 1
    }
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it comfort to how swith in R works, i.e. indexing begins with 1

    switch (deriv,
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- invlink(eta)
        L1 <- lambda / 2
        -sum(weight * (z * log(L1 / (1 + L1)) + (1 - z) * log(1 / (1 + L1))))
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- invlink(eta)
        L1 <- lambda / 2
        if (NbyK) {
          return(as.data.frame(X) * (z - L1 / (1 + L1)) * weight)
        }
        if (vectorDer) {
          return(matrix((z - L1 / (1 + L1)) * weight, ncol = 1))
        }
        t(X) %*% ((z - L1 / (1 + L1)) * weight)
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- invlink(eta)
        L1 <- lambda / 2
        -t(as.data.frame(X) * weight * (L1 / ((1 + L1) ^ 2))) %*% as.matrix(X)
      }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(1 > mu)
  }

  devResids <- function(y, eta, wt, ...) {
    z <- y - 1
    mu <- invlink(eta)
    mu1 <- mu.eta(eta = eta)
    ((-1) ^ y) * sqrt(-2 * wt * (z * log(mu1) + (1 - z) * log(1 - mu1)))
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
    N <- ((1 + 1 / (lambda + (lambda ^ 2) / 2)) * pw)
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    Xvlm <- as.data.frame(Xvlm)
    prob <- lambda * exp(-lambda) + (lambda ^ 2) * exp(-lambda) / 2

    f1 <- colSums(-Xvlm * pw * ((lambda + (lambda ^ 2)) /
                  ((lambda + (lambda ^ 2) / 2) ^ 2)))
    f1 <- t(f1) %*% as.matrix(cov) %*% f1

    f2 <- sum(pw * (1 - prob) * ((1 + exp(-lambda) / prob) ^ 2))

    f1 + f2
  }
  
  simulate <- function(n, eta, lower = 0, upper = 2) {
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
      link = "2 * log",
      valideta = function (eta) {TRUE},
      variance = variance,
      Wfun = Wfun,
      funcZ = funcZ,
      devResids = devResids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      simulate = simulate,
      family = "chao",
      parNum = 1,
      etaNames = "lambda",
      densityFunction = dFun,
      getStart = getStart
    ),
    class = "family"
  )
}

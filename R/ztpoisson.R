#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @importFrom stats glm.fit
#' @importFrom stats poisson
#' @export
ztpoisson <- function(lambdaLink = c("log", "neglog"),
                      ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink)
  
  lambdaLink <- switch (lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  links[1] <- c(lambdaLink)
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = lambda / (1 - exp(-lambda))
    )
  }

  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = mu.eta(eta = eta) * (1 + lambda - mu.eta(eta = eta))
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    Ey <- mu.eta(eta)
    
    G1 <- (Ey / lambda - 1 / (1 - exp(-lambda))) * 
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
    
    G11 <- (exp(2 * lambda) / (exp(lambda) - 1) ^ 2 -
    exp(lambda) / (exp(lambda) - 1) - Ey / lambda ^ 2) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
    matrix(-prior * (G11 + G1), ncol = 1, 
           dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    ((y / lambda - 1 / (1-exp(-lambda))) * 
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / weight)
  }

  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        lambda <- lambdaLink((as.matrix(X) %*% beta)[, 1], inverse = TRUE)
        
        -sum(weight * (y * log(lambda) - log(exp(lambda) - 1) - log(factorial(y))))
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- (y / lambda - 1 / (1-exp(-lambda))) * weight * 
               lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
        
        if (NbyK) {
          return(as.data.frame(X) * G1)
        }
        if (vectorDer) {
          return(matrix(G1, ncol = 1))
        }
        t(as.matrix(X)) %*% G1
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- (y / lambda - 1 / (1 - exp(-lambda))) * 
               lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
        
        G11 <- (exp(2*lambda) / (exp(lambda) - 1) ^ 2 -
                exp(lambda) / (exp(lambda) - 1) - y / lambda ^ 2) *
                lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
        
        t(as.matrix(X) * weight * (G1 + G11)) %*% as.matrix(X)
      }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  devResids <- function(y, eta, wt, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    #hm1y <- ifelse(y > 1, VGAM::lambertW(-y * exp(-y)) + y, 0)
    hm1y <- ifelse(y > 1, lamW::lambertW0(-y * exp(-y)) + y, 0)
    sign(y - mu.eta(eta = eta)) * sqrt(-2 * wt * (y * log(lambda) - lambda - log(1 - exp(-lambda)) - y * ifelse(y > 1, log(hm1y), 0) + hm1y + ifelse(y > 1, log(1 - exp(-hm1y)), 0)))
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    N <- pw / (1 - exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    Xvlm <- as.data.frame(Xvlm)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)

    bigTheta <- -(lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * pw * 
                  exp(lambda) / (exp(lambda) - 1) ^ 2) %*% as.matrix(Xvlm)
    bigTheta <- as.vector(bigTheta)
    
    f1 <- t(bigTheta) %*% as.matrix(cov) %*% bigTheta

    f2 <- sum(pw * exp(-lambda) / ((1 - exp(-lambda)) ^ 2))

    f1 + f2
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    lb <- stats::ppois(lower, lambda)
    ub <- stats::ppois(upper, lambda)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qpois(p_u, lambda)
    sims
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    stats::dpois(x = x, lambda = lambda) / (1 - stats::dpois(x = 0, lambda = lambda))
  }
  
  getStart <- expression(
    start <- stats::glm.fit(
      x = variables[wch$reg, ],
      y = observed[wch$reg],
      family = stats::poisson(),
      weights = priorWeights[wch$reg],
      ...
    )$coefficients,
    if (attr(family$links, "linkNames")[1] == "neglog") start <- -start
  )
  
  structure(
    list(
      makeMinusLogLike = minusLogLike,
      densityFunction  = dFun,
      links     = links,
      mu.eta    = mu.eta,
      valideta  = function (eta) {TRUE},
      variance  = variance,
      Wfun      = Wfun,
      funcZ     = funcZ,
      devResids = devResids,
      validmu   = validmu,
      pointEst  = pointEst,
      popVar    = popVar,
      family    = "ztpoisson",
      etaNames  = c("lambda"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

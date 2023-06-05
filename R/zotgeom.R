#' @rdname singleRmodels
#' @export
zotgeom <- function(lambdaLink = c("log", "neglog"),
                    ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink)
  
  lambdaLink <- switch (lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  links[1] <- c(lambdaLink)
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    if (!deriv) {
      switch (type,
        "nontrunc" = lambda,
        "trunc" = 2 + lambda
      )
    } else {
      switch (type,
        "nontrunc" = lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
        "trunc" = lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
      )
    }
  }
  
  variance <- function(eta, disp, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    switch (type,
    "nontrunc" = lambda * (lambda - 1),
    "trunc" = lambda * (lambda + 1)
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    Ey <- mu.eta(eta)
    
    G1 <- -lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
          (lambda - Ey + 2) / (lambda ^ 2 + lambda)
    
    G11 <- lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 * 
           (lambda ^ 2 + (4 - 2 * Ey) * lambda - Ey + 2) / 
           (lambda ^ 2 * (lambda + 1) ^ 2)
    
    matrix(- prior * (G11 + G1), ncol = 1, 
    dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, mu, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    -lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
    (lambda - y + 2) / (lambda ^ 2 + lambda) / weight
  }
  
  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        -sum(weight * ((y - 2) * log(lambda) - (y - 1) * log(1 + lambda)))
      },
      function(beta) {
        eta <- X %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        G1 <- -lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
              (lambda - y + 2) / (lambda ^ 2 + lambda)
        if (NbyK) {
          return(G1  * weight * as.data.frame(X))
        }
        if (vectorDer) {
          return(matrix(G1  * weight, ncol = 1))
        }
        t(G1  * weight) %*% X
      },
      function(beta) {
        eta <- X %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- -lambdaLink(eta[, 1], inverse = TRUE, deriv = 2) *
              (lambda - y + 2) / (lambda ^ 2 + lambda)
        
        G11 <- lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 * 
               (lambda ^ 2 + (4 - 2 * y) * lambda - y + 2) / 
               (lambda ^ 2 * (lambda + 1) ^ 2)
        
        t(as.data.frame(X) * (G1 + G11) * weight) %*% X
      }
    )
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  devResids <- function (y, eta, wt, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    loghm1y <- ifelse(y > 2, log(y - 2), 0)
    sign(y - mu.eta(eta = eta)) * 
    sqrt(-2 * wt * ((y - 2) * log(lambda) - (y - 1) * log(1 + lambda) - 
                    (y - 2) * loghm1y + (y - 1) * log(y - 1)))
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    N <- pw * (lambda ^ 2 + lambda + 1) / lambda ^ 2
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    bigTheta <- t(Xvlm) %*% (pw * as.numeric(-(lambda + 2) / lambda ^ 3) *
                             lambdaLink(eta[, 1], inverse = TRUE, deriv = 1))
    
    bigTheta <- as.vector(bigTheta)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * (lambda ^ 2 + lambda + 1) / lambda ^ 2)
    
    f1 + f2
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    lb <- stats::pnbinom(lower, mu = lambda, size = 1)
    ub <- stats::pnbinom(upper, mu = lambda, size = 1)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qnbinom(p_u, mu = lambda, size = 1)
    sims
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        stats::dgeom(x = x, prob = 1 / (1 + lambda)) / 
        (1 - stats::dgeom(x = 0, prob = 1 / (1 + lambda)) - 
        stats::dgeom(x = 1, prob = 1 / (1 + lambda)))
      },
      "nontrunc" = stats::dgeom(x = x, prob = 1 / (1 + lambda))
    )
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
      family    = "zotgeom",
      etaNames  = c("lambda"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}
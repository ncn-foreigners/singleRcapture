#' @rdname singleRmodels
#' @export
ztgeom <- function(lambdaLink = c("log", "neglog"),
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
      switch (
        type,
        "nontrunc" = lambda,
        "trunc" = 1 + lambda
      )
    } else {
      switch (type,
        "nontrunc" = lambdaLink(eta, inverse = TRUE, deriv = 1),
        "trunc" = lambdaLink(eta, inverse = TRUE, deriv = 1)
      )
    }
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    mu <- mu.eta(eta)
    switch (type,
      nontrunc = mu ^ 2 - mu - 1 / mu + 2,
      trunc    = (mu + 1) / mu
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    matrix( # returning as matrix to keep type consistency
      1 / (lambda * (lambda + 1)) * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2,
      ncol = 1, dimnames = list(rownames(eta), c("lambda"))
      # 0 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)) first 
      # derivative always zero after we take expected value
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    ((y - 1) / lambda - y / (1 + lambda)) * 
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / weight
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
        
        -sum(weight * ((y - 1) * log(lambda) - y * log(1 + lambda)))
      },
      function(beta) {
        eta <- X %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        G1 <- (y - 1) / lambda - y / (1 + lambda)
        G1 <- weight * G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
        # Beta derivative
        if (NbyK) {
          return(as.data.frame(X) * G1)
        }
        if (vectorDer) {
          return(matrix(G1, ncol = 1))
        }
        t(G1) %*% X
      },
      function(beta) {
        eta <- X %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- (y - 1) / lambda - y / (1 + lambda)
        G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
        
        G11 <- y / (lambda + 1) ^ 2 - (y - 1) / lambda ^ 2
        G11 <- G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
        
        # second beta derivative
        
        t(as.data.frame(X) * (G1 + G11) * weight) %*% X
      },
    )
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  devResids <- function (y, eta, wt, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    hm1y <- y - 1 # that's an analytic inverse for geometric
    loghm1y <- ifelse(y > 1, log(hm1y), 0)
    sign(y - 1 - lambda) * sqrt(-2 * wt * ((y - 1) * log(lambda) - y * log(1 + lambda) - 
                                           (y - 1) * loghm1y + y * log(y)))
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    N <- pw / (1 - 1 / (1 + lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    pr <- 1 - 1 / (1 + lambda)
    
    bigTheta <- -(lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * pw / lambda ^ 2) %*% as.matrix(Xvlm)
    bigTheta <- as.vector(bigTheta)
    
    f1 <- t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <- sum(pw * (1 - pr) / (pr ^ 2))
    
    f1 + f2
  }
  
simulate <- function(n, eta, lower = 0, upper = Inf) {
  lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    lb <- stats::pnbinom(lower, mu=lambda, size = 1)
    ub <- stats::pnbinom(upper, mu=lambda, size = 1)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qnbinom(p_u, mu=lambda, size = 1)
    sims
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    switch (type,
      "trunc" = {
        stats::dgeom(x = x, prob = 1 / (1 + lambda)) / 
        (1 - stats::dgeom(x = 0, prob = 1 / (1 + lambda)))
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
      family    = "ztgeom",
      etaNames  = c("lambda"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}
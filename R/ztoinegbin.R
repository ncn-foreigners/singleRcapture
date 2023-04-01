#' @rdname singleRmodels
#' @importFrom stats uniroot
#' @importFrom stats dnbinom
#' @importFrom rootSolve multiroot
#' @export
ztoinegbin <- function(nSim = 1000, epsSim = 1e-8, 
                       lambdaLink = c("log", "neglog"), 
                       alphaLink = c("log", "neglog"),
                       omegaLink = "logit", ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  if (missing(alphaLink))  alphaLink <- "log"
  if (missing(omegaLink))  omegaLink <- "logit"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink, alphaLink, omegaLink)
  
  lambdaLink <- switch(lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  alphaLink <- switch(alphaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  omegaLink <- switch(omegaLink,
    "logit" = singleRinternallogitLink
  )
  
  links[1:3] <- c(lambdaLink, alphaLink, omegaLink)
  
  
  mu.eta <- function(eta, type = "trunc", ...) {
    #TODO
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    switch (
      type,
      nontrunc = 0,
      trunc = 0
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    #TODO
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    switch (
      type,
      nontrunc = 0,
      trunc = 0
    ) - mu.eta(eta = eta, type = type) ^ 2
  }
  
  compdigamma <- function(y, alpha) {
    #temp <- 0:(y-1)
    #sum(-(alpha ^ (-2)) / (temp + 1 / alpha))
    (-digamma(y + 1 / alpha) + digamma(1 / alpha)) / (alpha ^ 2)
  }
  
  
  comptrigamma <- function(y, alpha) {
    #temp <- 0:(y-1)
    #sum(-(temp ^ 2) / (((1 + temp * alpha)) ^ 2))
    (alpha ^ 2 * (1 - y) + 2 * alpha * digamma(y + 1 / alpha) + trigamma(y + 1 / alpha) - 2 * alpha * digamma(1 + 1 / alpha) - trigamma(1 + 1 / alpha)) / (alpha ^ 4)
  }
  
  # Computing the expected value of di/trigamma functions on (y + 1/alpha)
  
  compExpect <- function(eta) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    res <- res1 <- 0
    k <- 0
    finished <- c(FALSE, FALSE)
    while ((k < nSim) & !all(finished)) {
      k <- k + 1 # 1 is the first possible y value for 0 truncated distribution
      prob <- stats::dnbinom(x = k, size = 1 / alpha, mu = lambda) / (1 - P0)
      if (!is.finite(prob)) prob <- 0
      toAdd <- c(compdigamma(y = k, alpha = alpha), comptrigamma(y = k, alpha = alpha)) * prob
      res <- res + toAdd
      finished <- abs(toAdd) < epsSim
    }
    res
  }
  
  Wfun <- function(prior, eta, ...) {
    #TODO
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    matrix(0,
      dimnames = list(rownames(eta), c("lambda", "lambda:alpha", "lambda:omega",
                                       "lambda:alpha", "alpha", "alpha:omega",
                                       "lambda:omega", "alpha:omega", "omega")),
      ncol = 9
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    #TODO:: the whole thing
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    # weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
    #   matrix(as.numeric(weight[x, ]), ncol = 3)
    # })
    # 
    # pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
    #   #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
    #   xx <- solve(weight[[x]]) #more stable
    #   xx %*% uMatrix[x, ]
    # })
    # pseudoResid <- t(pseudoResid)
    # dimnames(pseudoResid) <- dimnames(eta)
    # pseudoResid
    
    NULL
  }
  
  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    indY1 <- as.numeric(y == 1)
    X <- as.matrix(X)
    
    if (!(deriv %in% c(0, 1, 2))) 
      stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it comfort to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
            function(beta) {
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3)
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              omega  <-  omegaLink(eta[, 3], inverse = TRUE)
              z <- 1 / alpha
              
              -sum(weight * (indY1 * log(omega + (1 - omega) *
              (lgamma(y + z) - lgamma(z) - log(factorial(y)) - 
              (y + z) * log(1 + lambda * alpha) + y * log(lambda * alpha) - 
              log(1 - (1 + lambda * alpha) ^ (-z)))) + 
              (1 - indY1) * (lgamma(y + z) - lgamma(z) - log(factorial(y)) - 
              (y + z) * log(1 + lambda * alpha) + y * log(lambda * alpha) - 
              log(1 - (1 + lambda * alpha) ^ (-z)))))
            },
            function(beta) {
              ## TODO
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3)
              
              G2 <- "lambda derivative"
              G1 <- "alpha derivative"
              G0 <- "omega derivative"
              
              if (NbyK) {
                return(cbind(
                  G2 * as.data.frame(
                    X[1:nrow(eta), 
                      1:(attr(X, "hwm")[1])]
                  ), 
                  G1 * as.data.frame(
                    X[(nrow(eta)+1):(2*nrow(eta)), 
                      (attr(X, "hwm")[1]+1):(attr(X, "hwm")[2])]
                  ), 
                  G0 * as.data.frame(
                    X[(2*nrow(eta)+1):(3*nrow(eta)), 
                      (attr(X, "hwm")[2]+1):(attr(X, "hwm")[3])]
                  )
                ))
              }
              if (vectorDer) {
                return(cbind(G2, G1, G0))
              }
              
              as.numeric(c(G2, G1, G0) %*% X)
              NULL
            },
            function(beta) {
              ## TODO
              predNumbers <- attr(X, "hwm")
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3)
              NULL
            }
    )
  }
  
  validmu <- function(mu) {
    all(is.finite(mu)) && all(0 < mu)
  }
  
  devResids <- function (y, eta, wt, ...) {
    #TODO
    # AAAAAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaa
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    #TODO
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    #TODO
  }
  
  dFun <- function (x, eta, type = "trunc") {
    #TODO
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    #TODO:: uzupełnić CDF
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    CDF <- function(x) {
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, exp(-lambda), 
      exp(-lambda) + omega * (1 - exp(-lambda)) + 
      (1 - omega) * (stats::ppois(x, lambda) - exp(-lambda)))))
    }
    lb <- CDF(lower)
    ub <- CDF(upper)
    p_u <- stats::runif(n, lb, ub)
    sims <- rep(0, n)
    cond <- CDF(sims) <= p_u
    while (any(cond)) {
      sims[cond] <- sims[cond] + 1
      cond <- CDF(sims) <= p_u
    }
    sims
  }
  
  getStart <- expression(
    #TODO
    start <- stats::glm.fit(
      x = variables[wch$reg, 1:attr(Xvlm, "hwm")[1]],
      y = observed[wch$reg],
      family = stats::poisson(),
      weights = priorWeights[wch$reg],
      ...
    )$coefficients,
    if (!is.null(controlMethod$alphaStart)) {
      start <- c(start, controlMethod$alphaStart)
    } else {
      if (controlModel$alphaFormula == ~ 1) {
        start <- c(start, log(abs(mean(observed[wch$reg] ^ 2) - mean(observed[wch$reg])) / (mean(observed[wch$reg]) ^ 2 + .25)))
      } else {
        cc <- colnames(Xvlm)
        cc <- cc[grepl(x = cc, pattern = "alpha$")]
        cc <- unlist(strsplit(x = cc, ":alpha"))
        cc <- sapply(cc, FUN = function(x) {
          ifelse(x %in% names(start), start[x], 0) # TODO: gosh this is terrible pick a better method
        })
        start <- c(start, cc)
      }
    },
    if (is.null(controlMethod$omegaStart)) {
      if (controlModel$omegaFormula == ~ 1) {
        omg <- (length(observed[wch$reg]) - sum(observed == 1)) / (sum(observed[wch$reg]) - length(observed[wch$reg]))
        start <- c(start, log(omg))
      } else {
        cc <- colnames(Xvlm)
        cc <- cc[grepl(x = cc, pattern = "omega$")]
        cc <- unlist(strsplit(x = cc, ":omega"))
        cc <- sapply(cc, FUN = function(x) {
          ifelse(x %in% names(start), start[x], 0) # TODO: gosh this is terrible pick a better method
        })
        start <- c(start, cc)
      }
    } else {
      start <- c(start, controlMethod$omegaStart)
    }
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
      family    = "ztoinegbin",
      etaNames  = c("lambda", "alpha", "omega"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

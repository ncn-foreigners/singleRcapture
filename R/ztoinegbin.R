#' @rdname singleRmodels
#' @importFrom stats uniroot
#' @importFrom stats dnbinom
#' @importFrom rootSolve multiroot
#' @export
ztoinegbin <- function(nSim = 1000, epsSim = 1e-8, 
                       lambdaLink = c("log", "neglog"), 
                       alphaLink = c("log", "neglog"),
                       omegaLink = c("logit", "cloglog"), ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  if (missing(alphaLink))  alphaLink <- "log"
  if (missing(omegaLink))  omegaLink <- "logit"
  ## TODO:: there is an error in W and variance/devaince are not yet completed
  ## -- apart from that this is completed
  
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
    switch (type,
      nontrunc = omega * (1 + alpha * lambda) ^ (-1 / alpha) + (1 - omega) * lambda,
      trunc = omega + (1 - omega) * lambda / (1 - (1 + alpha * lambda) ^ (-1 / alpha))
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    #TODO
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    switch (type,
      nontrunc = 2 * mu.eta(eta = eta, type = type) ^ 2,
      trunc = 2 * mu.eta(eta = eta, type = type) ^ 2
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
    (2 * (digamma(y + 1 / alpha) - digamma(1 / alpha)) * alpha +
       trigamma(y + 1 / alpha) - trigamma(1 / alpha)) / (alpha ^ 4)
  }
  
  # Computing the expected value of di/trigamma functions on (y + 1/alpha)
  
  compExpectG1 <- function(eta) {
    lambda <- lambdaLink(eta[1], inverse = TRUE)
    alpha  <-  alphaLink(eta[2], inverse = TRUE)
    omega  <-  omegaLink(eta[3], inverse = TRUE)
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    res <- res1 <- 0
    k <- 1 # 1 is the first possible y value for 0 truncated distribution 1 inflated
    finished <- c(FALSE, FALSE)
    while ((k < nSim) & !all(finished)) {
      k <- k + 1 # but here we compute the (1 - z) * psi function which takes 0 at y = 1
      prob <- (1 - omega) * stats::dnbinom(x = k, size = 1 / alpha, mu = lambda) / (1 - P0)
      if (!is.finite(prob)) prob <- 0
      toAdd <- c( compdigamma(y = k, alpha = alpha), 
                 comptrigamma(y = k, alpha = alpha)) * prob
      res <- res + toAdd
      finished <- abs(toAdd) < epsSim
    }
    res
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    
    Ey <- mu.eta(eta)
    z  <- omega + (1 - omega) *
    (log(1 / alpha) - (1 + 1 / alpha) * log(1 + lambda * alpha) + 
    1 * log(lambda * alpha) - log(1 - (1 + lambda * alpha) ^ (-1 / alpha)))
    XXX <- Ey - z
    
    Edig  <- apply(X = eta, MARGIN = 1, FUN = function(x) {compExpectG1(x)})
    Etrig <- Edig[2, ]
    Edig  <- Edig[1, ]
    
    
    # omega
    G0 <- (z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1))) /
    (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) -
    (1 - z) / (1 - omega)
    
    G00 <- -(z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1)) ^ 2) /
    (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) ^ 2 -
    (1 - z) / (1 - omega) ^ 2
    
    G00 <- G00 * omegaLink(eta[, 3], inverse = TRUE, deriv = 1) ^ 2 +
            G0 * omegaLink(eta[, 3], inverse = TRUE, deriv = 2)
    
    # omega alpha
    G01 <- -(z * lambda * (lambda * alpha + 1) ^ (1 / alpha) * 
    ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha ^ 2 -
    lambda * alpha)) / (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) *
    (omega * lambda * alpha + omega) + (1 - omega) * lambda) ^ 2)
    
    G01 <- G01 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
                 omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
    
    # omega lambda
    G02 <- (z * (lambda - 1) * (alpha * lambda + 1) ^ (1 / alpha)) /
    ((alpha * lambda + 1) ^ (1 / alpha) * (alpha * omega * lambda + omega) +
    (1 - omega) * lambda) ^ 2
    
    G02 <- G02 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                  omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
    
    # alpha
    
    G1 <- ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1))) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) +
    (1 - z) * ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
    ((lambda * alpha + 1) ^ (1 / alpha) * (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) +
    log(lambda * alpha + 1) / alpha ^ 2) + (lambda * (-(1 - z) / alpha - XXX)) /
    (lambda * alpha + 1) + XXX / alpha + Edig
    
    G11 <- (1 - z) * ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / 
    (alpha * (lambda * alpha + 1))) ^ 2 / ((lambda * alpha + 1) ^ (1 / alpha) *
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) + (log(lambda * alpha + 1) / alpha ^ 2 -
    lambda / (alpha * (lambda * alpha + 1))) ^ 2 / ((lambda * alpha + 1) ^ (2 / alpha) *
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2) + (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 +
    (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) + lambda ^ 2 / 
    (alpha * (lambda * alpha + 1) ^ 2)) / ((lambda * alpha + 1) ^ (1 / alpha) * 
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) - (2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
    (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1))) - 
    (lambda ^ 2 * (-(1 - z) / alpha - XXX)) / 
    (lambda * alpha + 1) ^ 2 - XXX / alpha ^ 2 + Etrig + 
    ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) ^ 2) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
    ((1 - omega) ^ 2 * z * lambda ^ 2 * (lambda * alpha + 1) ^ (-2 / alpha - 2) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) ^ 2) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) ^ 2 +
    ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * 
    (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 *(lambda * alpha + 1)) -
    (lambda ^ 2 * (-1 / alpha - 1)) / (lambda * alpha + 1) ^ 2)) /
    ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega)
    
    G11 <- G11 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
            G1 * alphaLink(eta[, 2], inverse = TRUE, deriv = 2)
    
    # alpha lambda
    G12 <- (z * ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) + (1 - omega) * lambda * (-1 / alpha - 1) *
    alpha * (lambda * alpha + 1) ^ (-1 / alpha - 2) * 
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 2)) /
    (lambda * alpha + 1)) + ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 2)) / alpha + 
    (1 - omega) * lambda * (-1 / alpha - 1) * (lambda * alpha + 1) ^ (-1 / alpha - 2))) /
    ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
    ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 
    (1 - omega) * lambda * (-1 / alpha - 1) * alpha * (lambda * alpha + 1) ^ (-1 / alpha - 2)) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1))) /((1 - omega) * lambda * 
    (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) ^ 2 +
    (1 - z) * (-((lambda * alpha + 1) ^ (-2 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1)))) /
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2 - ((lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) / 
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha)) + 1 / (alpha * (lambda * alpha + 1))) + 
    (-(1 - z) / alpha - XXX) / (lambda * alpha + 1) - 
    (lambda * (-(1 - z) / alpha - XXX) * alpha) / (lambda * alpha + 1) ^ 2
    
    G12 <- G12 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
      alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    #lambda
    G2 <- (1 - z) * (-(alpha * lambda + 1) ^ (-1 / alpha - 1) / 
    (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha))) + 
    (alpha * (-XXX - (1 - z) / alpha)) / (alpha * lambda + 1) + XXX / lambda + 
    (z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) +
    (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2))) /
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega)
    
    G22 <- (1 - z) * (-((-1 / alpha - 1) * alpha * (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
    (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) + (alpha * lambda + 1) ^ (-2 / alpha - 2) /
    (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2) - XXX / lambda ^ 2 - 
    (alpha ^ 2 * (-XXX - (1 - z) / alpha)) / (alpha * lambda + 1) ^ 2 + 
    (z * (2 * (-1 / alpha - 1) * alpha *
    (1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 2) +
    (-1 / alpha - 2) * (-1 / alpha - 1) * alpha ^ 2 * (1 - omega) * lambda * 
    (alpha * lambda + 1) ^ (-1 / alpha - 3))) / 
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) -
    (z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) + (-1 / alpha - 1) * 
    alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) ^ 2) /
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) ^ 2
    
    
    G22 <- G22 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 + 
            G2 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
    
    matrix(
      -c(G22 * prior, G12 * prior, G02 * prior,
         G12 * prior, G11 * prior, G01 * prior,
         G02 * prior, G01 * prior, G00 * prior),
      dimnames = list(
        rownames(eta), 
        c("lambda", "lambda:alpha", "lambda:omega",
          "lambda:alpha", "alpha", "alpha:omega",
          "lambda:omega", "alpha:omega", "omega")
      ),
      ncol = 9
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    z <- as.numeric(y == 1)
    
    dig <-  compdigamma(y = y, alpha = alpha)
    
    G0 <- (z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1))) /
    (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) -
    (1 - z) / (1 - omega)
    
    G1 <- ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1))) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) +
    (1 - z) * ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
    ((lambda * alpha + 1) ^ (1 / alpha) * (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) +
    log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) /
    (lambda * alpha + 1) + y / alpha + dig)
    
    G2 <- (1 - z) * (-(alpha * lambda + 1) ^ (-1 / alpha - 1) / 
    (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) + (alpha * (-y - 1 / alpha)) /
    (alpha * lambda + 1) + y / lambda) + 
    (z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) +
    (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2))) /
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega)
    
    G0 <- G0 *  omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
    G1 <- G1 *  alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    G2 <- G2 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    uMatrix <- matrix(c(G2, G1, G0), ncol = 3)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 3)
    })

    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- solve(weight[[x]]) #more stable
      xx %*% uMatrix[x, ]
    })
    
    pseudoResid <- t(pseudoResid)
    dimnames(pseudoResid) <- dimnames(eta)
    pseudoResid
  }
  
  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    z <- as.numeric(y == 1)
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
              
              -sum(weight * (z * log(omega + (1 - omega) *
              lambda / (1 + alpha * lambda) ^ (1 / alpha + 1)) + 
              (1 - z) * (log(1 - omega) + 
              lgamma(y + 1 / alpha) - lgamma(1 / alpha) - log(factorial(y)) - 
              (y + 1 / alpha) * log(1 + lambda * alpha) + y * log(lambda * alpha) - 
              log(1 - (1 + lambda * alpha) ^ (-1 / alpha)))))
            },
            function(beta) {
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3)
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              omega  <-  omegaLink(eta[, 3], inverse = TRUE)
              
              dig <-  compdigamma(y = y, alpha = alpha)
              
              G0 <- (z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1))) /
              (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) -
              (1 - z) / (1 - omega)
              
              G1 <- ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1))) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) +
              (1 - z) * ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
              ((lambda * alpha + 1) ^ (1 / alpha) * (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) +
              log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) /
              (lambda * alpha + 1) + y / alpha + dig)
              
              G2 <- (1 - z) * (-(alpha * lambda + 1) ^ (-1 / alpha - 1) / 
              (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) + (alpha * (-y - 1 / alpha)) /
              (alpha * lambda + 1) + y / lambda) + 
              (z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) +
              (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2))) /
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega)
              
              
              G0 <- G0 * weight *  omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
              G1 <- G1 * weight *  alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
              G2 <- G2 * weight * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
              
              if (NbyK) {
                return(cbind(
                  G2 * as.data.frame(
                    X[1:nrow(eta), 
                      1:(attr(X, "hwm")[1])]
                  ), 
                  G1 * as.data.frame(
                    X[(nrow(eta) + 1):(2 * nrow(eta)), 
                      (attr(X, "hwm")[1] + 1):(attr(X, "hwm")[2])]
                  ), 
                  G0 * as.data.frame(
                    X[(2 * nrow(eta) + 1):(3 * nrow(eta)), 
                      (attr(X, "hwm")[2] + 1):(attr(X, "hwm")[3])]
                  )
                ))
              }
              if (vectorDer) {
                return(cbind(G2, G1, G0))
              }
              
              as.numeric(c(G2, G1, G0) %*% X)
            },
            function(beta) {
              predNumbers <- attr(X, "hwm")
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3)
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              omega  <-  omegaLink(eta[, 3], inverse = TRUE)
              res <- matrix(nrow = length(beta), 
                            ncol = length(beta), 
                            dimnames = list(names(beta), names(beta)))
              
              trig <- comptrigamma(y = y, alpha = alpha)
              dig <-  compdigamma(y = y, alpha = alpha)
              
              # omega
              G0 <- (z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1))) /
              (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) -
              (1 - z) / (1 - omega)
              
              G00 <- -(z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1)) ^ 2) /
              (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) ^ 2 -
              (1 - z) / (1 - omega) ^ 2
              
              G00 <- G00 * omegaLink(eta[, 3], inverse = TRUE, deriv = 1) ^ 2 +
                      G0 * omegaLink(eta[, 3], inverse = TRUE, deriv = 2)
              
              # omega alpha
              G01 <- -(z * lambda * (lambda * alpha + 1) ^ (1 / alpha) * 
              ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha ^ 2 -
              lambda * alpha)) / (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) *
              (omega * lambda * alpha + omega) + (1 - omega) * lambda) ^ 2)
              
              G01 <- G01 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
                           omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
              
              # omega lambda
              G02 <- (z * (lambda - 1) * (alpha * lambda + 1) ^ (1 / alpha)) /
              ((alpha * lambda + 1) ^ (1 / alpha) * (alpha * omega * lambda + omega) +
              (1 - omega) * lambda) ^ 2
              
              G02 <- G02 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                            omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
              # alpha
              G1 <- ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1))) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) +
              (1 - z) * ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
              ((lambda * alpha + 1) ^ (1 / alpha) * (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) +
              log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) /
              (lambda * alpha + 1) + y / alpha + dig)
              
              G11 <- (1 - z) * ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / 
              (alpha * (lambda * alpha + 1))) ^ 2 / ((lambda * alpha + 1) ^ (1 / alpha) *
              (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) + (log(lambda * alpha + 1) / alpha ^ 2 -
              lambda / (alpha * (lambda * alpha + 1))) ^ 2 / ((lambda * alpha + 1) ^ (2 / alpha) *
              (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2) + (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 +
              (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) + lambda ^ 2 / 
              (alpha * (lambda * alpha + 1) ^ 2)) / ((lambda * alpha + 1) ^ (1 / alpha) * 
              (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) - (2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
              (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) - (lambda ^ 2 * (-1 / alpha - y)) / 
              (lambda * alpha + 1) ^ 2 - y / alpha ^ 2 + trig) + 
              ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) ^ 2) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
              ((1 - omega) ^ 2 * z * lambda ^ 2 * (lambda * alpha + 1) ^ (-2 / alpha - 2) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) ^ 2) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) ^ 2 +
              ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * 
              (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 *(lambda * alpha + 1)) -
              (lambda ^ 2 * (-1 / alpha - 1)) / (lambda * alpha + 1) ^ 2)) /
              ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega)
              
              G11 <- G11 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
                      G1 * alphaLink(eta[, 2], inverse = TRUE, deriv = 2)
              
              # alpha lambda
              G12 <- (z * ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) + (1 - omega) * lambda * (-1 / alpha - 1) *
              alpha * (lambda * alpha + 1) ^ (-1 / alpha - 2) * 
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 2)) /
              (lambda * alpha + 1)) + ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 2)) / alpha + 
              (1 - omega) * lambda * (-1 / alpha - 1) * (lambda * alpha + 1) ^ (-1 / alpha - 2))) /
              ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
              ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 
              (1 - omega) * lambda * (-1 / alpha - 1) * alpha * (lambda * alpha + 1) ^ (-1 / alpha - 2)) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1))) /((1 - omega) * lambda * 
              (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) ^ 2 +
              (1 - z) * (-((lambda * alpha + 1) ^ (-2 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1)))) /
              (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2 - ((lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) / 
              (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha)) + 1 / (alpha * (lambda * alpha + 1)) + 
              (-1 / alpha - y) / (lambda * alpha + 1) - 
              (lambda * (-1 / alpha - y) * alpha) / (lambda * alpha + 1) ^ 2)
              
              G12 <- G12 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                            alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
              
              #lambda
              G2 <- (1 - z) * (-(alpha * lambda + 1) ^ (-1 / alpha - 1) / 
              (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) + (alpha * (-y - 1 / alpha)) /
              (alpha * lambda + 1) + y / lambda) + 
              (z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) +
              (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2))) /
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega)
              
              G22 <- (1 - z) * (-((-1 / alpha - 1) * alpha * (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
              (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) + (alpha * lambda + 1) ^ (-2 / alpha - 2) /
              (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2 - (alpha ^ 2 * (-y - 1 / alpha)) /
              (alpha * lambda + 1) ^ 2 - y / lambda ^ 2) + (z * (2 * (-1 / alpha - 1) * alpha *
              (1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 2) +
              (-1 / alpha - 2) * (-1 / alpha - 1) * alpha ^ 2 * (1 - omega) * lambda * 
              (alpha * lambda + 1) ^ (-1 / alpha - 3))) / 
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) -
              (z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) + (-1 / alpha - 1) * 
              alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) ^ 2) /
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) ^ 2
               
              G22 <- G22 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 + 
                      G2 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
              
              
              predNumbers <- cumsum(predNumbers)
              res[(predNumbers[2]+1):predNumbers[3], (predNumbers[2]+1):predNumbers[3]] <- #omega
                t(as.data.frame(X[(nrow(X) * 2 / 3 + 1):nrow(X), (predNumbers[2]+1):predNumbers[3]]) *
                  G00 * weight) %*% X[(nrow(X) * 2 / 3 + 1):nrow(X), (predNumbers[2]+1):predNumbers[3]]
              
              res[(predNumbers[1]+1):predNumbers[2], (predNumbers[2]+1):predNumbers[3]] <- #omega alpha
                t(as.data.frame(X[(nrow(X) / 3 + 1):(nrow(X) * 2 / 3), (predNumbers[1] + 1):predNumbers[2]]) * 
                  G01 * weight) %*% as.matrix(X[(nrow(X) * 2/ 3 + 1):nrow(X), (predNumbers[2]+1):predNumbers[3]])
              
              res[(predNumbers[2]+1):predNumbers[3], (predNumbers[1]+1):predNumbers[2]] <- #omega alpha
                t(t(as.data.frame(X[(nrow(X) / 3 + 1):(nrow(X) * 2 / 3), (predNumbers[1] + 1):predNumbers[2]]) * 
                    G01 * weight) %*% as.matrix(X[(nrow(X) * 2/ 3 + 1):nrow(X), (predNumbers[2]+1):predNumbers[3]]))
              
              res[1:predNumbers[1], (predNumbers[2]+1):predNumbers[3]] <- #omega lambda
                t(as.data.frame(X[1:(nrow(X) / 3), 1:predNumbers[1]]) * 
                  G02 * weight) %*% as.matrix(X[(nrow(X) * 2 / 3 + 1):nrow(X), (predNumbers[2]+1):predNumbers[3]])
              
              res[(predNumbers[2]+1):predNumbers[3], 1:predNumbers[1]] <- #omega lambda
                t(t(as.data.frame(X[1:(nrow(X) / 3), 1:predNumbers[1]]) * 
                    G02 * weight) %*% as.matrix(X[(nrow(X) * 2 / 3 + 1):nrow(X), (predNumbers[2]+1):predNumbers[3]]))
              
              res[(predNumbers[1]+1):predNumbers[2], (predNumbers[1]+1):predNumbers[2]] <- #alpha
                t(as.data.frame(X[(nrow(X) / 3 + 1):(nrow(X) * 2 / 3), (predNumbers[1]+1):predNumbers[2]]) *
                  G11 * weight) %*% X[(nrow(X) / 3 + 1):(nrow(X) * 2 / 3), (predNumbers[1]+1):predNumbers[2]]
              
              res[1:predNumbers[1], 1:predNumbers[1]] <-  # lambda
                t(as.data.frame(X[1:(nrow(X) / 3), 1:predNumbers[1]]) * 
                  G22 * weight) %*% X[1:(nrow(X) / 3), 1:predNumbers[1]]
              
              res[1:predNumbers[1], (predNumbers[1]+1):predNumbers[2]] <- #alpha lambda
                t(as.data.frame(X[1:(nrow(X) / 3), 1:predNumbers[1]]) * 
                  G12 * weight) %*% as.matrix(X[(nrow(X) / 3 + 1):(nrow(X) * 2 / 3), (predNumbers[1]+1):predNumbers[2]])
              
              res[(predNumbers[1]+1):predNumbers[2], 1:predNumbers[1]] <- #alpha lambda
                t(t(as.data.frame(X[1:(nrow(X) / 3), 1:predNumbers[1]]) * 
                    G12 * weight) %*% as.matrix(X[(nrow(X) / 3 + 1):(nrow(X) * 2 / 3), (predNumbers[1]+1):predNumbers[2]]))
              
              res
            }
    )
  }
  
  validmu <- function(mu) {
    all(is.finite(mu)) && all(0 < mu)
  }
  
  devResids <- function (y, eta, wt, ...) {
    #TODO
    # AAAAAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaa
    0
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    N <- pw / (1 - (1 + alpha * lambda) ^ (- 1 / alpha))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    pr <- 1 - (1 + alpha * lambda) ^ (- 1 / alpha)
    
    bigTheta0 <- pw * 0 # w.r to omega
    bigTheta1 <- pw  *  alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
    ((lambda * alpha + 1) ^ (1 / alpha - 1) * ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha)) /
    (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)# w.r to alpha
    bigTheta2 <- -pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
    ((alpha * lambda + 1) ^ (1 / alpha - 1) / ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)# w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1, bigTheta0) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * (1 - pr) / (pr ^ 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    
    ifelse(x == 1, 
           omega + (1 - omega) * lambda / (1 + alpha * lambda) ^ (1 / alpha + 1), 
           (1 - omega) * stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / (1 - P0))
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    CDF <- function(x) {
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, (1 + alpha * lambda) ^ (- 1 / alpha), 
      (1 + alpha * lambda) ^ (- 1 / alpha) + 
      omega * (1 - (1 + alpha * lambda) ^ (- 1 / alpha)) + 
      (1 - omega) * (stats::pnbinom(q = x, mu = lambda, size = 1 / alpha) -
                     (1 + alpha * lambda) ^ (- 1 / alpha)))))
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
    if (attr(family$links, "linkNames")[1] == "neglog") start <- -start,
    if (!is.null(controlMethod$alphaStart)) {
      start <- c(start, controlMethod$alphaStart)
    } else {
      if (controlModel$alphaFormula == ~ 1) {
        start <- c(start, family$links[[2]](abs(mean(observed[wch$reg] ^ 2) - mean(observed[wch$reg])) / (mean(observed[wch$reg]) ^ 2 + .25)))
      } else {
        cc <- colnames(Xvlm)
        cc <- cc[grepl(x = cc, pattern = "alpha$")]
        cc <- unlist(strsplit(x = cc, ":alpha"))
        cc <- sapply(cc, FUN = function(x) {
          ifelse(x %in% names(start), start[x], 0) # TODO: gosh this is terrible pick a better method
        })
        if (attr(family$links, "linkNames")[1] == attr(family$links, "linkNames")[2])
          start <- c(start, cc)
        else
          start <- c(start, -cc)
      }
    },
    if (controlModel$omegaFormula == ~ 1) {
      omg <- (length(observed[wch$reg]) - sum(observed == 1)) / (sum(observed[wch$reg]) - length(observed[wch$reg]))
      start <- c(start, family$links[[3]](omg / (1 - omg)))
    } else {
      cc <- colnames(Xvlm)
      cc <- cc[grepl(x = cc, pattern = "omega$")]
      cc <- unlist(strsplit(x = cc, ":omega"))
      cc <- sapply(cc, FUN = function(x) {
        ifelse(x %in% names(start), start[x], 0) # TODO: gosh this is terrible pick a better method
      })
      start <- c(start, cc)
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

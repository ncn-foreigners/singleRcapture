#' @rdname singleRmodels
#' @importFrom stats uniroot
#' @importFrom stats dnbinom
#' @importFrom rootSolve multiroot
#' @export
oiztnegbin <- function(nSim = 1000, epsSim = 1e-8, eimStep = 6,
                       lambdaLink = c("log", "neglog"), 
                       alphaLink = c("log", "neglog"),
                       omegaLink = c("logit", "cloglog", "probit"), ...) {
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
    "logit" = singleRinternallogitLink,
    "cloglog" = singleRinternalcloglogLink,
    "probit" = singleRinternalprobitLink
  )
  
  links[1:3] <- c(lambdaLink, alphaLink, omegaLink)
  
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    
    if (!deriv) {
      switch (type,
        "nontrunc" =  omega + (1 - omega) * lambda,
        "trunc"    = (omega + (1 - omega) * lambda) / 
        (1 - (1 - omega) * (1 + alpha * lambda) ^ (-1 / alpha))
      )
    } else {
      switch (type,
        "nontrunc" = {
          matrix(c(
            1 - omega,
            0,
            1 - lambda
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
             alphaLink(eta[, 2], inverse = TRUE, deriv = 1),
             omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
          ), ncol = 3)
        },
        "trunc" = {
          matrix(c(
            (1 - omega) * (alpha * lambda + 1) ^ (1 / alpha - 1) *
            ((alpha * lambda + 1) ^ (1 / alpha + 1) + 
            ((alpha + 1) * omega - alpha - 1) * lambda - 1) /
            ((alpha * lambda + 1) ^ (1 / alpha) + omega - 1) ^ 2,
            (1 - omega) * ((1 - omega) * lambda + omega) * 
            (lambda * alpha + 1) ^ (1 / alpha - 1) * 
            ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha) /
            (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) + omega - 1) ^ 2),
            -(alpha * lambda + 1) ^ (1 / alpha) *
            ((lambda - 1) * (alpha * lambda + 1) ^ (1 / alpha) + 1) /
            (omega + (alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
             alphaLink(eta[, 2], inverse = TRUE, deriv = 1),
             omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
          ), ncol = 3)
        }
      )
    }
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    switch (type,
    nontrunc =  omega + (1 - omega) * lambda * (1 + lambda * (1 + alpha)),
    trunc    = (omega + (1 - omega) * lambda * (1 + lambda * (1 + alpha))) /
    (1 - (1 - omega) * (1 + alpha * lambda) ^ (-1 / alpha))      
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
    P0 <- (1 - omega) * (1 + alpha * lambda) ^ (-1 / alpha)
    res <- c(0, 0)
    k <- 2 # 1 is the first possible y value for 0 truncated distribution 1 inflated
    # but here we compute the (1 - z) * psi function which takes 0 at y = 1
    finished <- c(FALSE, FALSE)
    while ((k < nSim) & !all(finished)) {
      prob <- (1 - omega) * stats::dnbinom(x = k:(k + eimStep), 
                                           size = 1 / alpha, 
                                           mu = lambda) / (1 - P0)
      if (any(!is.finite(prob))) {prob <- 0}
      toAdd <- cbind(compdigamma(y = k:(k + eimStep), alpha = alpha),
                     comptrigamma(y = k:(k + eimStep), alpha = alpha)) * prob
      toAdd <- colSums(toAdd)
      k <- k + eimStep + 1
      res <- res + toAdd
      finished <- abs(toAdd) < epsSim
    }
    res
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    
    z  <- (omega + (1 - omega) * lambda * (1 + alpha * lambda) ^ (-1 / alpha - 1)) / 
    (1 - (1 - omega) * (1 + alpha * lambda) ^ (-1 / alpha))
    
    XXX <- mu.eta(eta, type = "trunc") - z
    
    Edig  <- apply(X = eta, MARGIN = 1, FUN = function(x) {compExpectG1(x)})
    Etrig <- Edig[2, ]
    Edig  <- Edig[1, ]
    
    # omega
    G0 <- z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1)) /
    (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) - 
    1 / ((alpha * lambda + 1) ^ (1 / alpha) * 
    (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))) - 
    (1 - z) / (1 - omega)
    
    G00 <- -z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha- 1 )) ^ 2 /
    (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) ^ 2 + 
    1 / ((alpha * lambda + 1) ^ (2 / alpha) * 
    (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2) -
    (1 - z) / (1 - omega) ^ 2
    
    G00 <- G00 * omegaLink(eta[, 3], inverse = TRUE, deriv = 1) ^ 2 +
            G0 * omegaLink(eta[, 3], inverse = TRUE, deriv = 2)
    
    # omega alpha
    G01 <- -(log(lambda * alpha + 1) / alpha ^ 2 - 
    lambda / (alpha * (lambda * alpha + 1))) /
    ((lambda * alpha + 1) ^ (1 / alpha) * 
    (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) -
    ((1 - omega) * (log(lambda * alpha + 1) / alpha ^ 2 -
    lambda / (alpha * (lambda * alpha + 1)))) /
    ((lambda * alpha + 1) ^ (2 / alpha) * 
    (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2) -
    z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
    ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * 
    (1 - lambda * (lambda * alpha + 1) ^ (-1 / alpha -1 )) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1))) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) ^ 2
    
    G01 <- G01 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
                 omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
    
    # omega lambda
    G02 <- (alpha * lambda + 1) ^ (-1 / alpha - 1) / 
    (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha)) +
    ((1 - omega) * (alpha * lambda + 1) ^ (-2 / alpha - 1)) / 
    (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2 +
    z * (-(alpha * lambda + 1) ^ (-1 / alpha - 1) - 
    (-1 / alpha - 1) * alpha * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) /
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) -
    z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 
    (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) *
    (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1)) / 
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) ^ 2
    
    G02 <- G02 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                  omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
    # alpha
    G1 <- (1 - omega) * (log(lambda * alpha + 1) / alpha ^ 2 - 
    lambda / (alpha * (lambda * alpha + 1))) / 
    ((lambda * alpha + 1) ^ (1 / alpha) * (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) +
    ((1 - z) * log(lambda * alpha + 1) / alpha ^ 2 + 
    (lambda * (-(1 - z) / alpha - XXX)) / (lambda * alpha + 1) + XXX / alpha + Edig) + 
    (1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega)
    
    G11 <- (1 - omega) * (log(lambda * alpha + 1) / alpha ^ 2 - 
    lambda / (alpha * (lambda * alpha + 1))) ^ 2 / 
    ((lambda * alpha + 1) ^ (1 / alpha) * (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) +
    (1 - omega) ^ 2 * (log(lambda * alpha + 1) / alpha ^ 2 - 
    lambda / (alpha * (lambda * alpha+ 1 ))) ^ 2 / ((lambda * alpha + 1) ^ (2 / alpha) * 
    (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2) + 
    ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) ^ 2) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
    (1 - omega) ^ 2 * z * lambda ^ 2 * (lambda * alpha + 1) ^ (-2 / alpha - 2) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) ^ 2 / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha-1) + omega) ^ 2 +
    ((1 - omega) * (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 +
    (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) + 
    lambda ^ 2 / (alpha * (lambda * alpha + 1) ^ 2))) /
    ((lambda * alpha + 1) ^ (1 / alpha) * (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) +
    (-(1 - z) * (2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
    (1 - z) * (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) - 
    (lambda ^ 2 * (-(1 - z) / alpha - XXX)) / (lambda * alpha + 1) ^ 2 - XXX / alpha ^ 2 + Etrig) + 
    ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) -
    (lambda ^ 2 * (-1 / alpha - 1)) / (lambda * alpha + 1) ^ 2)) /
    ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega)
    
    G11 <- G11 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
            G1 * alphaLink(eta[, 2], inverse = TRUE, deriv = 2)
    
    # alpha lambda
    G12 <- z * ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) + (1 - omega) * lambda * (-1 / alpha - 1) *
    alpha * (lambda * alpha + 1) ^ (-1 / alpha - 2) * 
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 2)) /
    (lambda * alpha + 1)) + ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 2)) / alpha +
    (1 - omega) * lambda * (-1 / alpha - 1) * (lambda * alpha + 1) ^ (-1 / alpha - 2)) /
    ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
    ((1 - omega) ^ 2 * (lambda * alpha + 1) ^ (-2 / alpha- 1 ) * 
    (log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1)))) /
    (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2 -
    ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1))) / (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha)) -
    ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * 
    ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 
    (1 - omega) * lambda * (-1 / alpha - 1) * alpha * (lambda * alpha + 1) ^ (-1 / alpha - 2)) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
    ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) ^ 2 +
    ((1 - z) / (alpha * (lambda * alpha + 1)) + (-(1 - z) / alpha - XXX) / (lambda * alpha + 1) -
    (lambda * (-(1 - z) / alpha - XXX) * alpha) / (lambda * alpha + 1) ^ 2)
    
    G12 <- G12 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                  alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    #lambda
    G2 <- z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) +
    (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) /
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) +
    ((alpha * (-XXX - (1 - z) / alpha)) / (alpha * lambda + 1) + XXX / lambda) - 
    ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1)) / 
    (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))
    
    G22 <- (1 - omega) ^ 2 * (alpha * lambda + 1) ^ (-2 / alpha - 2) /
    (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2 +
    z * (2 * (-1 / alpha - 1) * alpha * (1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 2) +
    (-1 / alpha - 2) * (-1 / alpha - 1) * alpha ^ 2 * (1 - omega) * 
    lambda * (alpha * lambda + 1) ^ (-1 / alpha - 3)) / 
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) -
    z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha- 1 ) +
    (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) ^ 2 /
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) ^ 2 +
    (-(alpha ^ 2 * (-XXX - (1 - z) / alpha)) / (alpha * lambda + 1) ^ 2 - XXX / lambda ^ 2) -
    ((-1 / alpha - 1) * alpha * (1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 2)) /
    (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))
    
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
    
    G0 <- z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1)) /
    (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) - 
    1 / ((alpha * lambda + 1) ^ (1 / alpha) * 
    (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))) - 
    (1 - z) / (1 - omega)
    
    G1 <- (1 - omega) * (log(lambda * alpha + 1) / alpha ^ 2 - 
    lambda / (alpha * (lambda * alpha + 1))) / ((lambda * alpha + 1) ^ (1 / alpha) * 
    (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) +
    (1 - z) * (log(lambda * alpha + 1) / alpha ^ 2 + 
    (lambda * (-1 / alpha - y)) / (lambda * alpha + 1) + y / alpha + dig) + 
    (1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega)
    
    G2 <- z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) +
    (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) /
    ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) +
    (1 - z) * ((alpha * (-y - 1 / alpha)) / (alpha * lambda + 1) + y / lambda) - 
    ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1)) / 
    (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))
    
    G0 <- G0 *  omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
    G1 <- G1 *  alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    G2 <- G2 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    uMatrix <- cbind(G2, G1, G0)
    
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
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
            function(beta) {
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3)
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              omega  <-  omegaLink(eta[, 3], inverse = TRUE)
              
              -sum(weight * (z * log(omega + (1 - omega) *
              lambda * (1 + alpha * lambda) ^ (-1 / alpha - 1)) + 
              (1 - z) * (log(1 - omega) + lgamma(y + 1 / alpha) - 
              lgamma(1 / alpha) - lgamma(y + 1) - 
              (y + 1 / alpha) * log(1 + lambda * alpha) + 
              y * log(lambda * alpha)) -
              log(1 - (1 - omega) * (1 + alpha * lambda) ^ (-1 / alpha))))
            },
            function(beta) {
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3)
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              omega  <-  omegaLink(eta[, 3], inverse = TRUE)
              
              dig <-  compdigamma(y = y, alpha = alpha)
              
              G0 <- z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1)) /
              (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) - 
              1 / ((alpha * lambda + 1) ^ (1 / alpha) * 
              (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))) - 
              (1 - z) / (1 - omega)
              
              G1 <- (1 - omega) * (log(lambda * alpha + 1) / alpha ^ 2 - 
              lambda / (alpha * (lambda * alpha + 1))) / 
              ((lambda * alpha + 1) ^ (1 / alpha) * (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) +
              (1 - z) * (log(lambda * alpha + 1) / alpha ^ 2 + 
              (lambda * (-1 / alpha - y)) / (lambda * alpha + 1) + y / alpha + dig) + 
              (1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega)
              
              G2 <- z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) +
              (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) /
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) +
              (1 - z) * ((alpha * (-y - 1 / alpha)) / (alpha * lambda + 1) + y / lambda) - 
              ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1)) / 
              (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))
              
              
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
              G0 <- z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1)) /
              (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) - 
              1 / ((alpha * lambda + 1) ^ (1 / alpha) * 
              (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))) - 
              (1 - z) / (1 - omega)
              
              G00 <- -z * (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha- 1 )) ^ 2 /
              (omega + lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) * (1 - omega)) ^ 2 + 
              1 / ((alpha * lambda + 1) ^ (2 / alpha) * 
              (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2) -
              (1 - z) / (1 - omega) ^ 2
              
              G00 <- G00 * omegaLink(eta[, 3], inverse = TRUE, deriv = 1) ^ 2 +
                      G0 * omegaLink(eta[, 3], inverse = TRUE, deriv = 2)
              
              # omega alpha
              G01 <- -(log(lambda * alpha + 1) / alpha ^ 2 - 
              lambda / (alpha * (lambda * alpha + 1))) /
              ((lambda * alpha + 1) ^ (1 / alpha) * 
              (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) -
              ((1 - omega) * (log(lambda * alpha + 1) / alpha ^ 2 -
              lambda / (alpha * (lambda * alpha + 1)))) /
              ((lambda * alpha + 1) ^ (2 / alpha) * 
              (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2) -
              z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
              ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * 
              (1 - lambda * (lambda * alpha + 1) ^ (-1 / alpha -1 )) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1))) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) ^ 2
              
              G01 <- G01 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
                           omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
              
              # omega lambda
              G02 <- (alpha * lambda + 1) ^ (-1 / alpha - 1) / 
              (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha)) +
              ((1 - omega) * (alpha * lambda + 1) ^ (-2 / alpha - 1)) / 
              (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2 +
              z * (-(alpha * lambda + 1) ^ (-1 / alpha - 1) - 
              (-1 / alpha - 1) * alpha * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) /
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) -
              z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 
              (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) *
              (1 - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1)) / 
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) ^ 2
              
              G02 <- G02 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                            omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
              # alpha
              G1 <- (1 - omega) * (log(lambda * alpha + 1) / alpha ^ 2 - 
              lambda / (alpha * (lambda * alpha + 1))) / 
              ((lambda * alpha + 1) ^ (1 / alpha) * (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) +
              (1 - z) * (log(lambda * alpha + 1) / alpha ^ 2 + 
              (lambda * (-1 / alpha - y)) / (lambda * alpha + 1) + y / alpha + dig) + 
              (1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega)
              
              G11 <- (1 - omega) * (log(lambda * alpha + 1) / alpha ^ 2 - 
              lambda / (alpha * (lambda * alpha + 1))) ^ 2 / 
              ((lambda * alpha + 1) ^ (1 / alpha) * (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) +
              (1 - omega) ^ 2 * (log(lambda * alpha + 1) / alpha ^ 2 - 
              lambda / (alpha * (lambda * alpha+ 1 ))) ^ 2 / ((lambda * alpha + 1) ^ (2 / alpha) * 
              (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2) + 
              ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) ^ 2) / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
              (1 - omega) ^ 2 * z * lambda ^ 2 * (lambda * alpha + 1) ^ (-2 / alpha - 2) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) ^ 2 / ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha-1) + omega) ^ 2 +
              ((1 - omega) * (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 +
              (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) + 
              lambda ^ 2 / (alpha * (lambda * alpha + 1) ^ 2))) /
              ((lambda * alpha + 1) ^ (1 / alpha) * (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha))) +
              (1 - z) * (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
              (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) - 
              (lambda ^ 2 * (-1 / alpha - y)) / (lambda * alpha + 1) ^ 2 - y / alpha ^ 2 +
              trig) + ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) -
              (lambda ^ 2 * (-1 / alpha - 1)) / (lambda * alpha + 1) ^ 2)) /
              ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega)
              
              G11 <- G11 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
                      G1 * alphaLink(eta[, 2], inverse = TRUE, deriv = 2)
              
              # alpha lambda
              G12 <- z * ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) + (1 - omega) * lambda * (-1 / alpha - 1) *
              alpha * (lambda * alpha + 1) ^ (-1 / alpha - 2) * 
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 2)) /
              (lambda * alpha + 1)) + ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 2)) / alpha +
              (1 - omega) * lambda * (-1 / alpha - 1) * (lambda * alpha + 1) ^ (-1 / alpha - 2)) /
              ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) -
              ((1 - omega) ^ 2 * (lambda * alpha + 1) ^ (-2 / alpha- 1 ) * 
              (log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1)))) /
              (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2 -
              ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1))) / (1 - (1 - omega) / (lambda * alpha + 1) ^ (1 / alpha)) -
              ((1 - omega) * z * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * 
              ((1 - omega) * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 
              (1 - omega) * lambda * (-1 / alpha - 1) * alpha * (lambda * alpha + 1) ^ (-1 / alpha - 2)) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
              ((1 - omega) * lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + omega) ^ 2 +
              (1 - z) * (1 / (alpha * (lambda * alpha + 1)) + (-1 / alpha - y) / (lambda * alpha + 1) -
              (lambda * (-1 / alpha- y) * alpha) / (lambda * alpha + 1) ^ 2)
              
              G12 <- G12 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                            alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
              
              #lambda
              G2 <- z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1) +
              (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) /
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) +
              (1 - z) * ((alpha * (-y - 1 / alpha)) / (alpha * lambda + 1) + y / lambda) - 
              ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 1)) / 
              (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))
              
              G22 <- (1 - omega) ^ 2 * (alpha * lambda + 1) ^ (-2 / alpha - 2) /
              (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2 +
              z * (2 * (-1 / alpha - 1) * alpha * (1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 2) +
              (-1 / alpha - 2) * (-1 / alpha - 1) * alpha ^ 2 * (1 - omega) * 
              lambda * (alpha * lambda + 1) ^ (-1 / alpha - 3)) / 
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) -
              z * ((1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha- 1 ) +
              (-1 / alpha - 1) * alpha * (1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) ^ 2 /
              ((1 - omega) * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + omega) ^ 2 +
              (1 - z) * (-(alpha ^ 2 * (-y - 1 / alpha)) / (alpha * lambda + 1) ^ 2 - y / lambda ^ 2) -
              ((-1 / alpha - 1) * alpha * (1 - omega) * (alpha * lambda + 1) ^ (-1 / alpha - 2)) /
              (1 - (1 - omega) / (alpha * lambda + 1) ^ (1 / alpha))
              
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
    0
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    N <- pw / (1 - (1 - omega) * (1 + alpha * lambda) ^ (- 1 / alpha))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    pr <- 1 - (1 - omega) * (1 + alpha * lambda) ^ (- 1 / alpha)
    
    # w.r to omega
    bigTheta0 <- -pw * omegaLink(eta[, 3], inverse = TRUE, deriv = 1) * 
      (alpha * lambda + 1) ^ (1 / alpha) / 
      (omega + (alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2
    # w.r to alpha
    bigTheta1 <- pw  *  alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
      ((1 - omega) * (lambda * alpha + 1) ^ (1 / alpha - 1) *
      ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha)) /
      (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) + omega - 1) ^ 2)
    # w.r to lambda
    bigTheta2 <- pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
      (((omega - 1) * (alpha * lambda + 1) ^ (1 / alpha - 1)) /
      ((alpha * lambda + 1) ^ (1 / alpha) + omega - 1) ^ 2)
    
    bigTheta <- t(c(bigTheta2, bigTheta1, bigTheta0) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * (1 - pr) / (pr ^ 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    P0 <- (1 - omega) * (1 + alpha * lambda) ^ (-1 / alpha)
    
    switch (type,
      "trunc" = ifelse(x == 1, 
        (omega + (1 - omega) * 
        stats::dnbinom(x = 1, mu = lambda, size = 1 / alpha)) / (1 - P0), 
        (1 - omega) * stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / (1 - P0)
      ),
      "nontrunc" = ifelse(x == 0, 
        (1 - omega) * P0, 
        ifelse(x == 1, 
        omega + (1 - omega) * stats::dnbinom(x = 1, mu = lambda, size = 1 / alpha), 
        (1 - omega) * stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / (1 - P0)
        )
      )
    )
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    P0 <- (1 - omega) * (1 + alpha * lambda) ^ (-1 / alpha)
    
    CDF <- function(x) {
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, (1 - omega) * P0, 
      omega + (1 - omega) * 
      stats::pnbinom(q = x, mu = lambda, size = 1 / alpha))))
    }
    lb <- CDF(lower)
    ub <- CDF(upper)
    p_u <- stats::runif(n, lb, ub)
    sims <- rep(0, n)
    cond <- CDF(sims) < p_u
    while (any(cond)) {
      sims[cond] <- sims[cond] + 1
      cond <- CDF(sims) < p_u
    }
    sims
  }
  
  # new starting points
  
  getStart <- expression(
    if (!is.null(controlMethod$start)) {
      start <- controlMethod$start
    } else {
      init <- c(
        family$links[[1]](mean(observed)),
        family$links[[2]](abs((var(observed) / mean(observed) - 1) / mean(observed)) + .1),
        family$links[[3]](mean(observed == 1) + .01)
      )
      if (attr(terms, "intercept")) {
        start <- c(init[1], rep(0, attr(Xvlm, "hwm")[1] - 1))
      } else {
        start <- rep(init[1] / attr(Xvlm, "hwm")[1], attr(Xvlm, "hwm")[1])
      }
      if ("(Intercept):alpha" %in% colnames(Xvlm)) {
        start <- c(start, init[2], rep(0, attr(Xvlm, "hwm")[2] - 1))
      } else {
        start <- c(start, rep(init[2] / attr(Xvlm, "hwm")[2], attr(Xvlm, "hwm")[2]))
      }
      if ("(Intercept):omega" %in% colnames(Xvlm)) {
        start <- c(start, init[3], rep(0, attr(Xvlm, "hwm")[3] - 1))
      } else {
        start <- c(start, rep(init[3] / attr(Xvlm, "hwm")[3], attr(Xvlm, "hwm")[3]))
      }
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
      family    = "oiztnegbin",
      etaNames  = c("lambda", "alpha", "omega"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

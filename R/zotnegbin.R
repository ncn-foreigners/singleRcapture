#' @rdname singleRmodels
#' @export
zotnegbin <- function(nSim = 1000, epsSim = 1e-8, eimStep = 6, 
                      lambdaLink = c("log", "neglog"), 
                      alphaLink = c("log", "neglog"),
                      ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  if (missing(alphaLink))  alphaLink <- "log"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink, alphaLink)
  
  lambdaLink <- switch(lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  alphaLink <- switch(alphaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  links[1:2] <- c(lambdaLink, alphaLink)
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    if (!deriv) {
      switch (
        type,
        "nontrunc" = lambda,
        "trunc" = (lambda - lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha))) / 
        (1 - (1 + alpha * lambda) ^ (-1 / alpha) - lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha)))
      )
    } else {
      switch (
        type,
        "nontrunc" = {
          matrix(c(1, 0) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
             alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        },
        "trunc" = {
          matrix(c(
            ((alpha * lambda + 1) ^ (2 / alpha) * (alpha ^ 2 * lambda ^ 2 + 2 * alpha * lambda + 1) +
            (alpha * lambda + 1) ^ (1 / alpha) * ((-alpha ^ 2 - 2 * alpha - 1) * lambda ^ 2 - 2 * alpha * lambda - 2) + 1) /
            ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2,
            lambda ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) * 
            (lambda * alpha ^ 2 + (lambda + 1) * alpha + 1) * log(lambda * alpha + 1) +
            (lambda * alpha + 1) ^ (1 / alpha) * ((1 - 2 * lambda) * alpha ^ 2 - lambda * alpha) - alpha ^ 2) /
            (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha + 1) - lambda * alpha - lambda - 1) ^ 2)
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
             alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        }
      )
    }
  }

  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    switch (type,
    "nontrunc" = lambda * (1 + alpha * lambda),
    "trunc" = (lambda * (1 + alpha * lambda) + lambda ^ 2 - 
    lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha))) / 
    (1 - (1 + alpha * lambda) ^ (-1 / alpha) - 
    lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha))) - 
    mu.eta(eta, type = "trunc") ^ 2
    )
  }
  
  compExpect <- function(eta) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    P1 <- lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1)
    #P0 <- stats::dnbinom(x = 0, size = 1 / alpha, mu = lambda)
    res <- rep(0, NROW(eta))
    k <- 2 
    finished <- rep(FALSE, NROW(eta))
    while ((k < nSim) & !all(finished)) {
      prob <- apply(cbind(k:(k + eimStep)), MARGIN = 1, FUN = function(x) {
        stats::dnbinom(
          x = x, 
          size = 1 / alpha, 
          mu = lambda
        ) / (1 - P0 - P1)
      })
      trg <- apply(cbind(k:(k + eimStep)), MARGIN = 1, FUN = function(x) {
        (2 * (digamma(x + 1 / alpha) - digamma(1 / alpha)) * alpha +
           trigamma(x + 1 / alpha) - trigamma(1 / alpha)) / alpha ^ 4
      })
      prob[!(is.finite(prob))] <- 0
      trg[!(is.finite(trg))] <- 0
      toAdd <- trg * prob
      toAdd <- rowSums(toAdd)
      k <- k + eimStep + 1
      res <- res + toAdd
      finished <- abs(toAdd) < epsSim
    }
    res
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    Ey <- mu.eta(eta = eta)
    prob <- 1 - (1 + lambda * alpha) ^ (-1 / alpha) - 
    lambda * (1 + lambda * alpha) ^ (-1 - 1 / alpha)
    
    Etrig <- compExpect(eta)
    
    # 2nd log(alpha) derivative
    
    G00 <- Etrig + (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
    (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) ^ 2 / 
    (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) ^ 2 -
    (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) ^ 2 / 
    (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) ^ 2 - (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
    (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) + lambda ^ 2 / 
    (alpha * (lambda * alpha + 1) ^ 2)) / (lambda * alpha + 1) ^ (1 / alpha) -
    lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * 
    (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) -
    (lambda ^ 2 * (-1 / alpha - 1)) / (lambda * alpha + 1) ^ 2)) / 
    (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) -
    (2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) -
    (lambda ^ 2 * (-1 / alpha - Ey)) / (lambda * alpha + 1) ^ 2 - Ey / alpha ^ 2
    
    G00 <- G00 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2
    
    # mixed derivative
    G01 <- -((alpha * lambda + 1) ^ (1 / alpha) * ((alpha ^ 3 + alpha ^ 2) * lambda ^ 3 +
    (2 * alpha ^ 2 + 2 * alpha) * lambda ^ 2 + (alpha + 1) * lambda) * 
    log(alpha * lambda + 1) + (alpha * lambda + 1) ^ (1 / alpha) * 
    ((alpha ^ 4 - alpha ^ 3 - alpha ^ 2) * lambda ^ 3 + 
    ((-2 * alpha ^ 4 - 2 * alpha ^ 3) * Ey + 4 * alpha ^ 3 - alpha ^ 2 - alpha) * lambda ^ 2 + 
    ((-4 * alpha ^ 3 - 2 * alpha ^ 2) * Ey + 3 * alpha ^ 2) * lambda - 2 * alpha ^ 2 * Ey) +
    (alpha * lambda + 1) ^ (2 / alpha) * (-alpha ^ 4 * lambda ^ 3 + 
    (alpha ^ 4 * Ey - 2 * alpha ^ 3) * lambda ^ 2 + 
    (2 * alpha ^ 3 * Ey - alpha ^ 2) * lambda + alpha ^ 2 * Ey) +
    ((alpha ^ 4 + 2 * alpha ^ 3 + alpha ^ 2) * Ey - 2 * alpha ^ 3 - alpha ^ 2) * lambda ^ 2 +
    ((2 * alpha ^ 3 + 2 * alpha ^ 2) * Ey - 2 * alpha ^ 2) * lambda + alpha ^ 2 * Ey) /
    (alpha ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2)
    
    G01 <- G01 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                  alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    # second beta derivative
    G11 <- ((alpha * lambda + 1) ^ (2 / alpha) * 
    (alpha ^ 3 * lambda ^ 4 + (2 * alpha ^ 2 - 2 * alpha ^ 3 * Ey) * lambda ^ 3 +
    (alpha - 5 * alpha ^ 2 * Ey) * lambda ^ 2 - 4 * alpha * Ey * lambda - Ey) + 
    (alpha * lambda + 1) ^ (1 / alpha) * ((alpha - alpha ^ 3) * lambda ^ 4 +
    ((4 * alpha ^ 3 + 4 * alpha ^ 2) * Ey - 4 * alpha ^ 2 - alpha + 1) * lambda ^ 3 +
    ((10 * alpha ^ 2 + 6 * alpha) * Ey - 3 * alpha - 1) * lambda ^ 2 + 
    (8 * alpha + 2) * Ey * lambda + 2 * Ey) + 
    ((-2 * alpha ^ 3 - 4 * alpha ^ 2 - 2 * alpha) * Ey + 2 * alpha ^ 2 + 2 * alpha) * lambda ^ 3 +
    ((-5 * alpha ^ 2 - 6 * alpha - 1) * Ey + 2 * alpha + 1) * lambda ^ 2 + 
    (-4 * alpha - 2) * Ey * lambda - Ey) / 
    (lambda ^ 2 * (alpha * lambda + 1) ^ 2 * 
    ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2)
    
    G11 <- G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
    matrix(
      -c(G11 * prior, # lambda predictor derivative without X matrix,
        G01 * prior,  # mixed derivative without X matrix
        G01 * prior,  # mixed derivative without X matrix
        G00 * prior  # alpha predictor derivative without X matrix
      ),
      dimnames = list(rownames(eta), 
                      c("lambda", "mixed", "mixed", "alpha")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    G0 <- y / alpha + (digamma(1 / alpha) - digamma(1 / alpha + y)) / alpha ^ 2 +
    ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
    (lambda * alpha + 1) ^ (1 / alpha) + lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1)) +
    log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) / (lambda * alpha + 1)
    
    G0 <- G0 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    G1 <- y / lambda - ((1 / alpha + 1) * alpha * lambda * 
    (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
    (-1 / (alpha * lambda + 1) ^ (1 / alpha) - 
    lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 1) -
    (alpha * (y + 1 / alpha)) / (alpha * lambda + 1)
    
    G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- solve(weight[[x]]) # more stable
      xx %*% uMatrix[x, ]
    })
    pseudoResid <- t(pseudoResid)
    dimnames(pseudoResid) <- dimnames(eta)
    pseudoResid
  }

  minusLogLike <- function(y, X, 
                           weight    = 1, 
                           NbyK      = FALSE, 
                           vectorDer = FALSE, 
                           deriv     = 0,
                           offset, 
                           ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X) / 2), rep(0, NROW(X) / 2))
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)

        -sum(weight * (lgamma(y + 1 / alpha) - lgamma(1 / alpha) -
        lgamma(y + 1) - (y + 1 / alpha) * log(1 + lambda * alpha) +
        y * log(lambda * alpha) - log(1 - (1 + lambda * alpha) ^ (-1 / alpha) - 
        lambda * (1 + lambda * alpha) ^ (-1 - 1 / alpha)))
        )
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
        
        G0 <- y / alpha + (digamma(1 / alpha) - digamma(1 / alpha + y)) / alpha ^ 2  +
        log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) / (lambda * alpha + 1) -
        (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
        (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
        (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
        (-(lambda * alpha + 1) ^ (-1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1)
        
        G0 <- G0 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) * weight
        
        G1 <- y / lambda  + (alpha * (-y - 1 / alpha)) / (alpha * lambda + 1) -
        ((1 / alpha + 1) * alpha * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
        (-(alpha * lambda + 1) ^ (-1 / alpha) - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 1)
        
        G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * weight
        
        if (NbyK) {
          XX <- 1:(attr(X, "hwm")[1])
          return(cbind(as.data.frame(X[1:nrow(eta), XX]) * G1, 
                       as.data.frame(X[-(1:nrow(eta)), -XX]) * G0))
        }
        if (vectorDer) {
          return(cbind(G1, G0))
        }
        
        as.numeric(c(G1, G0) %*% X)
      },
      function(beta) {
        lambdaPredNumber <- attr(X, "hwm")[1]
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
        res <- matrix(nrow = length(beta), ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        G0 <- y / alpha + (digamma(1 / alpha) - digamma(1 / alpha + y)) / alpha ^ 2 -
        (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
        (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
        (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
        (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) +
        log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) / (lambda * alpha + 1)
        
        G1 <- y / lambda + ((-1 / alpha-1) * alpha * lambda * 
        (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
        (-1 / (alpha * lambda + 1) ^ (1 / alpha) - 
        lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 1) + 
        (alpha * (-y - 1 / alpha)) / (alpha * lambda + 1)
        
        # 2nd log(alpha) derivative
        
        G00 <- (2 * digamma(1 / alpha + y)) / alpha ^ 3 - (2 * digamma(1 / alpha)) / alpha ^ 3 +
        trigamma(1 / alpha + y) / alpha ^ 4 - trigamma(1 / alpha) / alpha ^ 4 + 
        (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
        (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
        (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) ^ 2 / 
        (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) ^ 2 -
        (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) ^ 2 / 
        (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
        (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
        (lambda * alpha + 1)) ^ 2 - (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
        (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) + lambda ^ 2 / 
        (alpha * (lambda * alpha + 1) ^ 2)) / (lambda * alpha + 1) ^ (1 / alpha) -
        lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * 
        (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) -
        (lambda ^ 2 * (-1 / alpha - 1)) / (lambda * alpha + 1) ^ 2)) / 
        (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) -
        (2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) -
        (lambda ^ 2 * (-1 / alpha - y)) / (lambda * alpha + 1) ^ 2 - y / alpha ^ 2
          
        G00 <- G00 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
               G0  * alphaLink(eta[, 2], inverse = TRUE, deriv = 2)
        
        # mixed derivative
        G01 <- -((alpha * lambda + 1) ^ (1 / alpha) * ((alpha ^ 3 + alpha ^ 2) * lambda ^ 3 +
        (2 * alpha ^ 2 + 2 * alpha) * lambda ^ 2 + (alpha + 1) * lambda) * 
        log(alpha * lambda + 1) + (alpha * lambda + 1) ^ (1 / alpha) * 
        ((alpha ^ 4 - alpha ^ 3 - alpha ^ 2) * lambda ^ 3 + 
        ((-2 * alpha ^ 4 - 2 * alpha ^ 3) * y + 4 * alpha ^ 3 - alpha ^ 2 - alpha) * lambda ^ 2 + 
        ((-4 * alpha ^ 3 - 2 * alpha ^ 2) * y + 3 * alpha ^ 2) * lambda - 2 * alpha ^ 2 * y) +
        (alpha * lambda + 1) ^ (2 / alpha) * (-alpha ^ 4 * lambda ^ 3 + 
        (alpha ^ 4 * y - 2 * alpha ^ 3) * lambda ^ 2 + (2 * alpha ^ 3 * y - alpha ^ 2) * lambda + alpha ^ 2 * y) +
        ((alpha ^ 4 + 2 * alpha ^ 3 + alpha ^ 2) * y - 2 * alpha ^ 3 - alpha ^ 2) * lambda ^ 2 +
        ((2 * alpha ^ 3 + 2 * alpha ^ 2) * y - 2 * alpha ^ 2) * lambda + alpha ^ 2 * y) /
        (alpha ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2)
        
        
        G01 <- G01 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
               alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
        
        # second beta derivative
        G11 <- ((alpha * lambda + 1) ^ (2 / alpha) * 
        (alpha ^ 3 * lambda ^ 4 + (2 * alpha ^ 2 - 2 * alpha ^ 3 * y) * lambda ^ 3 +
        (alpha - 5 * alpha ^ 2 * y) * lambda ^ 2 - 4 * alpha * y * lambda - y) + 
        (alpha * lambda + 1) ^ (1 / alpha) * ((alpha - alpha ^ 3) * lambda ^ 4 +
        ((4 * alpha ^ 3 + 4 * alpha ^ 2) * y - 4 * alpha ^ 2 - alpha + 1) * lambda ^ 3 +
        ((10 * alpha ^ 2 + 6 * alpha) * y - 3 * alpha - 1) * lambda ^ 2 + 
        (8 * alpha + 2) * y * lambda + 2 * y) + 
        ((-2 * alpha ^ 3 - 4 * alpha ^ 2 - 2 * alpha) * y + 2 * alpha ^ 2 + 2 * alpha) * lambda ^ 3 +
        ((-5 * alpha ^ 2 - 6 * alpha - 1) * y + 2 * alpha + 1) * lambda ^ 2 + 
        (-4 * alpha - 2) * y * lambda - y) / 
        (lambda ^ 2 * (alpha * lambda + 1) ^ 2 * 
        ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2)
        
        G11 <- G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 +
               G1  * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
        
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]) * weight * G00) %*% X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * weight * G11) %*% X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * weight * G01) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]))
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * weight * G01) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        res
      }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  devResids <- function (y, eta, wt, ...) {
    # lambda <- exp(eta[, 1])
    # alpha <- exp(eta[, 2])
    
    # 
    # logLikFit <- (lgamma(y + 1/alpha) - lgamma(1/alpha) -
    # lgamma(y + 1) - (y + 1/alpha) * log(1+alpha * lambda) +
    # y * log(lambda * alpha) - log(1 - (1+alpha * lambda) ^ (-1/alpha) -
    # lambda * ((1+alpha * lambda) ^ (-1-1/alpha))))
    # 
    # 
    # yUnq <- unique(y) # see comments in zotpoisson
    # findL <- function(yNow) {
    #   root <- rootSolve::multiroot(
    #     start = c(1, mean(lambda), mean(alpha)),# maybe pick better starting points
    #     f = function(x) {# TODO:: provide analytic jacobian matrix will make it faster and more reliable
    #       s <- log(x[1]) # this is the lagrange multiplier and has no constraints of positivity
    #       l <- x[2]
    #       a <- x[3] # including constraints
    #       
    #       # c(log(l-l*((1+a*l)^(-1-1/a)))-log(1-(1+a*l)^(-1/a)-l*((1+a*l)^(-1-1/a)))-log(yNow),#sder
    #       #   yNow/l-(1+a*yNow)/(1+a*l)+s/l+s*(1+a)*((1+a*l)^(-2-1/a))/(1-(1+a*l)^(-1-1/a))-(1+s)*l*(1+a)*((1+a*l)^(-2-1/a))/(1-(1+a*l)^(-1/a)-l*((1+a*l)^(-1-1/a))),#lambda der
    #       #   (digamma(1/a)-digamma(yNow+1/a))/(a^2)-l*(yNow+1/a)/(1+a*l)+log(1+a*l)/(a^2)+
    #       #   yNow/a-s*((1+a*l)^(-1-1/a))*(log(1+a*l)/(a^2)-l*(1+1/a)/(1+a*l))/(1-(1+a*l)^(-1-1/a))+
    #       #   (1+s)*((log(1+a*l)/(a^2)-l/(a*(1+a*l)))/((1+a*l)^(1/a))+l*((1+a*l)^(-1-1/a))*
    #       #   (log(1+a*l)/(a^2)-l*(1+1/a)/(1+a*l)))/(1-(1+a*l)^(-1/a)-l*((1+a*l)^(-1-1/a))))#alpha der
    #       c(log(l-l*(a*l+1)^(-1/a-1))-log(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)-log(yNow),#sder
    #         -((-1/a-1)*a*(-s-1)*l*(a*l+1)^(-1/a-2))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)+(s*(-(a*l+1)^(-1/a-1)-(-1/a-1)*a*l*(a*l+1)^(-1/a-2)+1))/(l-l*(a*l+1)^(-1/a-1))+(a*(-yNow-1/a))/(a*l+1)+yNow/l,#lambda der
    #         ((-s-1)*(-(log(l*a+1)/a^2-l/(a*(l*a+1)))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1))))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)-(l*s*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1)))/(l-l*(l*a+1)^(-1/a-1))+(log(l*a+1)-digamma(1/a+yNow)+digamma(1/a))/a^2+(l*(-1/a-yNow))/(l*a+1)+yNow/a)
    #     }, maxiter = 10000, positive = TRUE
    #   )$root
    #   print(root)
    #   root <- log(root)
    #   
    #   (lgamma(yNow + exp(-root[3])) - lgamma(exp(-root[3])) -
    #   log(factorial(yNow)) - (yNow + exp(-root[3])) * log(1+exp(root[2] + root[3])) +
    #   yNow * (root[2]+root[3]) - log(1 - (1+exp(root[2] + root[3])) ^ (-exp(-root[3])) -
    #   exp(root[2]) * ((1+exp(root[2] + root[3])) ^ (-1-exp(-root[3])))))
    # }
    # logLikIdeal <- sapply(yUnq, FUN = function(x) {
    #   ifelse(x == 2, 0, findL(x))
    # })
    # 
    # logLikIdeal <- sapply(y, FUN = function(x) logLikIdeal[yUnq == x])
    # 
    # sign(y - mu.eta(eta = eta)) * sqrt(-2 * wt * (logLikFit - logLikIdeal))
    NULL
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    N <- pw * (1 - lambda * ((1 + lambda * alpha) ^ (-1 - 1 / alpha))) / 
    (1 - (1 + lambda * alpha) ^ (- 1 / alpha) - lambda * 
    ((1 + lambda * alpha) ^ (-1 - 1 / alpha)))
    
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    prob <- 1 - (1 + lambda * alpha) ^ (-1 / alpha) - 
      lambda * (1 + lambda * alpha) ^ (-1 - 1 / alpha)
    
    bigTheta1 <- pw * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) * 
    as.numeric(((lambda*alpha+1)^(1/alpha)*(lambda^2*alpha^2+2*lambda*alpha+1)*log(lambda*alpha+1)+(lambda*alpha+1)^(1/alpha)*(-lambda^2*alpha^2-lambda*alpha)-lambda^2*alpha^2)/(alpha^2*((lambda*alpha+1)^(1/alpha+1)-lambda*alpha-lambda-1)^2))
    
    bigTheta2 <- pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
    as.numeric(-((alpha*lambda+1)^(1/alpha+1)-1)/((alpha*lambda+1)^(1/alpha+1)+(-alpha-1)*lambda-1)^2)
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * ((1 - prob) / (prob ^ 2)) * 
    (1 - lambda * ((1 + lambda * alpha) ^ (-1 - 1 / alpha))) ^ 2)
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / 
        (1 - stats::dnbinom(x = 0, mu = lambda, size = 1 / alpha) - 
        stats::dnbinom(x = 1, mu = lambda, size = 1 / alpha))
      },
      "nontrunc" = stats::dnbinom(x = x, mu = lambda, size = 1 / alpha)
    )
  }

  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    lb <- stats::pnbinom(lower, mu = lambda, size = 1 / alpha)
    ub <- stats::pnbinom(upper, mu = lambda, size = 1 / alpha)
    p_u <- stats::runif(n, lb, ub)
    
    sims <- stats::qnbinom(p_u, mu = lambda, size = 1 / alpha)
    sims
  }
  
  getStart <- expression(
    if (method == "IRLS") {
      # init <- log(abs((observed / mean(observed) - 1) / mean(observed)) + .1)
      init <- log(abs((observed / mean(observed) - 1) / observed) + .1)
      etaStart <- cbind(
        pmin(family$links[[1]](observed), family$links[[1]](12)),
        family$links[[2]](ifelse(init < -.5, .1, init + .55))
      ) + offset
      etaStart <- etaStart[(observed > 1), , drop = FALSE]
    } else if (method == "optim") {
      init <- c(
        family$links[[1]](mean(observed)),
        family$links[[2]](abs((var(observed) / mean(observed) - 1) / mean(observed)) + .1)
      )
      if (attr(terms, "intercept")) {
        coefStart <- c(init[1], rep(0, attr(Xvlm, "hwm")[1] - 1))
      } else {
        coefStart <- rep(init[1] / attr(Xvlm, "hwm")[1], attr(Xvlm, "hwm")[1])
      }
      if ("(Intercept):alpha" %in% colnames(Xvlm)) {
        coefStart <- c(coefStart, init[2], rep(0, attr(Xvlm, "hwm")[2] - 1))
      } else {
        coefStart <- c(coefStart, rep(init[2] / attr(Xvlm, "hwm")[2], attr(Xvlm, "hwm")[2]))
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
      family    = "zotnegbin",
      etaNames  = c("lambda", "alpha"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

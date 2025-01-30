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
  
  Wfun <- function(prior, eta, y, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)[iddx]
    Ey <- mu.eta(eta = eta)[iddx]
    
    prob <- (1 - (1 + lambda * alpha) ^ (-1 / alpha) - 
    lambda * (1 + lambda * alpha) ^ (-1 - 1 / alpha))[iddx]
    
    Etrig <- compExpect(eta[iddx, , drop = FALSE])
    
    # 2nd alpha derivative
    G00 <- iddx
    G00[iddx] <- Etrig + (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
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
    
    G00[iddx] <- G00[iddx] * alphaLink(eta[, 2], inverse = TRUE, deriv = 1)[iddx] ^ 2
    
    # mixed derivative
    G01 <- iddx
    G01[iddx] <- -((alpha * lambda + 1) ^ (1 / alpha) * ((alpha ^ 3 + alpha ^ 2) * lambda ^ 3 +
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
    
    G01[iddx] <- G01[iddx] * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] *
      alphaLink(eta[, 2], inverse = TRUE, deriv = 1)[iddx]
    
    # second lambda derivative
    G11 <- iddx
    G11[iddx] <- ((alpha * lambda + 1) ^ (2 / alpha) * 
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
    
    G11[iddx] <- G11[iddx] * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] ^ 2
    
    matrix(
      -c(G11 * prior, # lambda
        G01 * prior,  # mixed 
        G01 * prior,  # mixed 
        G00 * prior  # alpha 
      ),
      dimnames = list(rownames(eta), 
                      c("lambda", "mixed", "mixed", "alpha")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)[iddx]
    weight[iddx, ] <- weight[iddx, ] / prior[iddx]
    
    G0 <- iddx
    G0[iddx] <- y[iddx] / alpha + (digamma(1 / alpha) - digamma(1 / alpha + y[iddx])) / alpha ^ 2 +
    ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
    (lambda * alpha + 1) ^ (1 / alpha) + lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1)) +
    log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y[iddx])) / (lambda * alpha + 1)
    
    G0[iddx] <- G0[iddx] * alphaLink(eta[, 2], inverse = TRUE, deriv = 1)[iddx]
    
    G1 <- iddx
    G1[iddx] <- y[iddx] / lambda - ((1 / alpha + 1) * alpha * lambda * 
    (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
    (-1 / (alpha * lambda + 1) ^ (1 / alpha) - 
    lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 1) -
    (alpha * (y[iddx] + 1 / alpha)) / (alpha * lambda + 1)
    
    G1[iddx] <- G1[iddx] * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx]
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      if (iddx[x]) {
        xx <- solve(weight[[x]])
        xx %*% uMatrix[x, ]
      } else {
        uMatrix[x, ]
      }
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
    iddx <- y > 1
    X <- as.matrix(X)

    if (!(deriv %in% c(0, 1, 2))) 
      stop("Only score function and derivatives up to 2 are supported.")
    
    # to make it conform to how switch in R works, i.e. indexing begins with 1
    deriv <- deriv + 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)

        -sum(weight * iddx * (lgamma(y + 1 / alpha) - lgamma(1 / alpha) -
        lgamma(y + 1) - (y + 1 / alpha) * log(1 + lambda * alpha) +
        y * log(lambda * alpha) - log(1 - (1 + lambda * alpha) ^ (-1 / alpha) - 
        lambda * (1 + lambda * alpha) ^ (-1 - 1 / alpha)))
        )
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)[iddx]
        
        G0 <- iddx
        G0[iddx] <- y[iddx] / alpha + (digamma(1 / alpha) - digamma(1 / alpha + y[iddx])) / alpha ^ 2  +
        log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y[iddx])) / (lambda * alpha + 1) -
        (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
        (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
        (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
        (-(lambda * alpha + 1) ^ (-1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1)
        
        G0[iddx] <- G0[iddx] * alphaLink(eta[, 2], inverse = TRUE, deriv = 1)[iddx] * weight[iddx]
        
        G1 <- iddx
        G1[iddx] <- y[iddx] / lambda  + (alpha * (-y[iddx] - 1 / alpha)) / (alpha * lambda + 1) -
        ((1 / alpha + 1) * alpha * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
        (-(alpha * lambda + 1) ^ (-1 / alpha) - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 1)
        
        G1[iddx] <- G1[iddx] * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] * weight[iddx]
        
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
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)[iddx]
        
        res <- matrix(nrow = length(beta), ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        G0 <- iddx
        G0[iddx] <- y[iddx] / alpha + (digamma(1 / alpha) - digamma(1 / alpha + y[iddx])) / alpha ^ 2 -
        (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
        (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
        (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
        (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) +
        log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y[iddx])) / (lambda * alpha + 1)
        
        G1 <- iddx
        G1[iddx] <- y[iddx] / lambda + ((-1 / alpha-1) * alpha * lambda * 
        (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
        (-1 / (alpha * lambda + 1) ^ (1 / alpha) - 
        lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 1) + 
        (alpha * (-y[iddx] - 1 / alpha)) / (alpha * lambda + 1)
        
        # 2nd alpha derivative
        
        G00 <- iddx
        G00[iddx] <- (2 * digamma(1 / alpha + y[iddx])) / alpha ^ 3 - (2 * digamma(1 / alpha)) / alpha ^ 3 +
        trigamma(1 / alpha + y[iddx]) / alpha ^ 4 - trigamma(1 / alpha) / alpha ^ 4 + 
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
        (lambda ^ 2 * (-1 / alpha - y[iddx])) / (lambda * alpha + 1) ^ 2 - y[iddx] / alpha ^ 2
          
        G00[iddx] <- G00[iddx] * alphaLink(eta[, 2], inverse = TRUE, deriv = 1)[iddx] ^ 2 +
                     G0[iddx]  * alphaLink(eta[, 2], inverse = TRUE, deriv = 2)[iddx]
        
        # mixed derivative
        G01 <- iddx
        G01[iddx] <- -((alpha * lambda + 1) ^ (1 / alpha) * ((alpha ^ 3 + alpha ^ 2) * lambda ^ 3 +
        (2 * alpha ^ 2 + 2 * alpha) * lambda ^ 2 + (alpha + 1) * lambda) * 
        log(alpha * lambda + 1) + (alpha * lambda + 1) ^ (1 / alpha) * 
        ((alpha ^ 4 - alpha ^ 3 - alpha ^ 2) * lambda ^ 3 + 
        ((-2 * alpha ^ 4 - 2 * alpha ^ 3) * y[iddx] + 4 * alpha ^ 3 - alpha ^ 2 - alpha) * lambda ^ 2 + 
        ((-4 * alpha ^ 3 - 2 * alpha ^ 2) * y[iddx] + 3 * alpha ^ 2) * lambda - 2 * alpha ^ 2 * y[iddx]) +
        (alpha * lambda + 1) ^ (2 / alpha) * (-alpha ^ 4 * lambda ^ 3 + 
        (alpha ^ 4 * y[iddx] - 2 * alpha ^ 3) * lambda ^ 2 + (2 * alpha ^ 3 * y[iddx] - alpha ^ 2) * lambda + alpha ^ 2 * y[iddx]) +
        ((alpha ^ 4 + 2 * alpha ^ 3 + alpha ^ 2) * y[iddx] - 2 * alpha ^ 3 - alpha ^ 2) * lambda ^ 2 +
        ((2 * alpha ^ 3 + 2 * alpha ^ 2) * y[iddx] - 2 * alpha ^ 2) * lambda + alpha ^ 2 * y[iddx]) /
        (alpha ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2)
        
        
        G01[iddx] <- G01[iddx] * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] *
                                  alphaLink(eta[, 2], inverse = TRUE, deriv = 1)[iddx]
        
        # second lambda derivative
        G11 <- iddx
        G11[iddx] <- ((alpha * lambda + 1) ^ (2 / alpha) * 
        (alpha ^ 3 * lambda ^ 4 + (2 * alpha ^ 2 - 2 * alpha ^ 3 * y[iddx]) * lambda ^ 3 +
        (alpha - 5 * alpha ^ 2 * y[iddx]) * lambda ^ 2 - 4 * alpha * y[iddx] * lambda - y[iddx]) + 
        (alpha * lambda + 1) ^ (1 / alpha) * ((alpha - alpha ^ 3) * lambda ^ 4 +
        ((4 * alpha ^ 3 + 4 * alpha ^ 2) * y[iddx] - 4 * alpha ^ 2 - alpha + 1) * lambda ^ 3 +
        ((10 * alpha ^ 2 + 6 * alpha) * y[iddx] - 3 * alpha - 1) * lambda ^ 2 + 
        (8 * alpha + 2) * y[iddx] * lambda + 2 * y[iddx]) + 
        ((-2 * alpha ^ 3 - 4 * alpha ^ 2 - 2 * alpha) * y[iddx] + 2 * alpha ^ 2 + 2 * alpha) * lambda ^ 3 +
        ((-5 * alpha ^ 2 - 6 * alpha - 1) * y[iddx] + 2 * alpha + 1) * lambda ^ 2 + 
        (-4 * alpha - 2) * y[iddx] * lambda - y[iddx]) / 
        (lambda ^ 2 * (alpha * lambda + 1) ^ 2 * 
        ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2)
        
        G11[iddx] <- G11[iddx] * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] ^ 2 +
                     G1[iddx]  * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)[iddx]
        
        
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
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    mu <- mu.eta(eta = eta)

    logLikFit <- (
      lgamma(y + 1 / alpha) - lgamma(1 / alpha) - lgamma(y + 1) - (y + 1 / alpha) * log(1 + lambda * alpha) +
      y * log(lambda * alpha) - log(1 - (1 + lambda * alpha) ^ (-1 / alpha) - lambda * (1 + lambda * alpha) ^ (-1 - 1 / alpha))
    )[iddx]

    yUnq <- unique(y[iddx])

    if (any(yUnq > 77)) {
      warning("Curently numerical deviance is unreliable for counts greater than 78.")
    }
    
    findL <- function(yNow) {
      stats::optim(
        par = if(TRUE) c(0, log(yNow), -15 + 2 * (yNow %in% 3:5) - 2 * (yNow %in% 6:9)) else c(0, log(yNow), -6),
        fn = function(x) {
          s <- x[1]
          l <- exp(x[2])
          a <- exp(x[3])
          
          sum(c(((l-l*(a*l+1)^(-1/a-1))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)-yNow) ,# s der
                s*((-(a*l+1)^(-1/a-1)-(-1/a-1)*a*l*(a*l+1)^(-1/a-2)+1)/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)+((-1/a-1)*a*l*(a*l+1)^(-1/a-2)*(l-l*(a*l+1)^(-1/a-1)))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)^2)+((-1/a-1)*a*l*(a*l+1)^(-1/a-2))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)+(a*(-yNow-1/a))/(a*l+1)+yNow/l,# lambda der
                s*(-((l-l*(l*a+1)^(-1/a-1))*(-(log(l*a+1)/a^2-l/(a*(l*a+1)))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1))))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)^2-(l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1)))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1))-(-(log(l*a+1)/a^2-l/(a*(l*a+1)))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1)))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)+log(l*a+1)/a^2+(l*(-1/a-yNow))/(l*a+1)+yNow/a-digamma(1/a+yNow)/a^2+digamma(1/a)/a^2) ^ 2)#alpha der
        },
        gr = function(x) {
          s <- x[1]
          l <- exp(x[2])
          a <- exp(x[3])
          d12.21 <- ((a*l+1)^(2/a)*(a^2*l^2+2*a*l+1)+(a*l+1)^(1/a)*((-a^2-2*a-1)*l^2-2*a*l-2)+1)/((a*l+1)^(1/a+1)+(-a-1)*l-1)^2
          d13.31 <- (l^2*((l*a+1)^(1/a)*(l*a^2+(l+1)*a+1)*log(l*a+1)+(l*a+1)^(1/a)*((1-2*l)*a^2-l*a)-a^2))/(a^2*((l*a+1)^(1/a+1)-l*a-l-1)^2)
          d22.22 <- s*((-2*(-1/a-1)*a*(a*l+1)^(-1/a-2)-(-1/a-2)*(-1/a-1)*a^2*l*(a*l+1)^(-1/a-3))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)+((-1/a-1)*a*(a*l+1)^(-1/a-2)*(l-l*(a*l+1)^(-1/a-1)))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)^2+((-1/a-2)*(-1/a-1)*a^2*l*(a*l+1)^(-1/a-3)*(l-l*(a*l+1)^(-1/a-1)))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)^2+(2*(-1/a-1)*a*l*(a*l+1)^(-1/a-2)*(-(a*l+1)^(-1/a-1)-(-1/a-1)*a*l*(a*l+1)^(-1/a-2)+1))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)^2+(2*(-1/a-1)^2*a^2*l^2*(a*l+1)^(-2/a-4)*(l-l*(a*l+1)^(-1/a-1)))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)^3)+((-1/a-1)*a*(a*l+1)^(-1/a-2))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)+((-1/a-2)*(-1/a-1)*a^2*l*(a*l+1)^(-1/a-3))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)+((-1/a-1)^2*a^2*l^2*(a*l+1)^(-2/a-4))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)^2-(a^2*(-yNow-1/a))/(a*l+1)^2-yNow/l^2
          d23.32 <- (((l*a+1)^(2/a)*(l^5*s*a^5+(5*l^4*s-l^4)*a^4+((-l^5+2*l^4+9*l^3)*s-l^4-3*l^3)*a^3+((-3*l^4+6*l^3+7*l^2)*s-3*l^3-3*l^2)*a^2+((-3*l^3+6*l^2+2*l)*s-3*l^2-l)*a+(2*l-l^2)*s-l)+(l*a+1)^(1/a)*(-l^5*s*a^5+((-3*l^5-5*l^4)*s+l^4)*a^4+((-3*l^5-10*l^4-9*l^3)*s+2*l^4+3*l^3)*a^3+((-l^5-7*l^4-13*l^3-7*l^2)*s+l^4+5*l^3+3*l^2)*a^2+((-2*l^4-5*l^3-8*l^2-2*l)*s+2*l^3+4*l^2+l)*a+(-l^3-l^2-2*l)*s+l^2+l))*log(l*a+1)+(l*a+1)^(2/a)*((3*l^3*yNow-l^5*s-2*l^4)*a^5+((3*l^3+9*l^2)*yNow+(2*l^5-8*l^4+2*l^3)*s-8*l^3)*a^4+((6*l^2+9*l)*yNow+(l^5+2*l^4-13*l^3+4*l^2)*s+l^4-10*l^2)*a^3+((3*l+3)*yNow+(2*l^4-2*l^3-6*l^2+2*l)*s+2*l^3-4*l)*a^2+((l^3-2*l^2)*s+l^2)*a)+(l*a+1)^(3/a)*((l^4-l^3*yNow)*a^5+(3*l^3-3*l^2*yNow)*a^4+(3*l^2-3*l*yNow)*a^3+(l-yNow)*a^2)+(l*a+1)^(1/a)*((-3*l^3*yNow+l^5*s+l^4)*a^5+((-6*l^3-9*l^2)*yNow+(3*l^5+8*l^4-4*l^3)*s+7*l^3)*a^4+((-3*l^3-12*l^2-9*l)*yNow+(3*l^5+7*l^4+13*l^3-8*l^2)*s-2*l^4+3*l^3+11*l^2)*a^3+((-3*l^2-6*l-3)*yNow+(l^5+4*l^4+6*l^3+6*l^2-4*l)*s-l^4-3*l^3+3*l^2+5*l)*a^2+((l^4+l^3+2*l^2)*s-l^3-l^2)*a)+l^3*yNow*a^5+((3*l^3+3*l^2)*yNow+2*l^3*s-2*l^3)*a^4+((3*l^3+6*l^2+3*l)*yNow+4*l^2*s-3*l^3-4*l^2)*a^3+((l^3+3*l^2+3*l+1)*yNow+2*l*s-l^3-3*l^2-2*l)*a^2)/(a^2*(l*a+1)^2*((l*a+1)^(1/a+1)-l*a-l-1)^3)
          d33.33 <- s*((2*(l-l*(l*a+1)^(-1/a-1))*(-(log(l*a+1)/a^2-l/(a*(l*a+1)))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1)))^2)/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)^3-((l-l*(l*a+1)^(-1/a-1))*(-(log(l*a+1)/a^2-l/(a*(l*a+1)))^2/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1))^2-(-(2*log(l*a+1))/a^3+(2*l)/(a^2*(l*a+1))+l^2/(a*(l*a+1)^2))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(-(2*log(l*a+1))/a^3+(2*l)/(a^2*(l*a+1))-(l^2*(-1/a-1))/(l*a+1)^2)))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)^2-(l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1))^2)/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)+(2*l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1))*(-(log(l*a+1)/a^2-l/(a*(l*a+1)))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1))))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)^2-(l*(l*a+1)^(-1/a-1)*(-(2*log(l*a+1))/a^3+(2*l)/(a^2*(l*a+1))-(l^2*(-1/a-1))/(l*a+1)^2))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1))+(-(log(l*a+1)/a^2-l/(a*(l*a+1)))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1)))^2/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)^2-(-(log(l*a+1)/a^2-l/(a*(l*a+1)))^2/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1))^2-(-(2*log(l*a+1))/a^3+(2*l)/(a^2*(l*a+1))+l^2/(a*(l*a+1)^2))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(-(2*log(l*a+1))/a^3+(2*l)/(a^2*(l*a+1))-(l^2*(-1/a-1))/(l*a+1)^2))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)-(2*log(l*a+1))/a^3+(2*l)/(a^2*(l*a+1))-(l^2*(-1/a-yNow))/(l*a+1)^2-yNow/a^2+(2*digamma(1/a+yNow))/a^3-(2*digamma(1/a))/a^3+trigamma(1/a+yNow)/a^4-trigamma(1/a)/a^4
          
          f2 <- 2*c(((l-l*(a*l+1)^(-1/a-1))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)-yNow) ,# s der
                    s*((-(a*l+1)^(-1/a-1)-(-1/a-1)*a*l*(a*l+1)^(-1/a-2)+1)/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)+((-1/a-1)*a*l*(a*l+1)^(-1/a-2)*(l-l*(a*l+1)^(-1/a-1)))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)^2)+((-1/a-1)*a*l*(a*l+1)^(-1/a-2))/(-1/(a*l+1)^(1/a)-l*(a*l+1)^(-1/a-1)+1)+(a*(-yNow-1/a))/(a*l+1)+yNow/l,# lambda der
                    s*(-((l-l*(l*a+1)^(-1/a-1))*(-(log(l*a+1)/a^2-l/(a*(l*a+1)))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1))))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)^2-(l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1)))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1))-(-(log(l*a+1)/a^2-l/(a*(l*a+1)))/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)*(log(l*a+1)/a^2+(l*(-1/a-1))/(l*a+1)))/(-1/(l*a+1)^(1/a)-l*(l*a+1)^(-1/a-1)+1)+log(l*a+1)/a^2+(l*(-1/a-yNow))/(l*a+1)+yNow/a-digamma(1/a+yNow)/a^2+digamma(1/a)/a^2)
          
          c(sum(f2 * c(0, d12.21, d13.31)),
            sum(f2 * c(d12.21, d22.22, d23.32)),
            sum(f2 * c(d13.31, d23.32, d33.33)))
        },
        method = "CG",
        control = list(maxit = 50000, abstol = .Machine$double.eps)
      )
    }
    fff <- function(yNow) {
      if (yNow == 2) return(0)
      
      xx <- findL(yNow)
      a <- exp(xx$par[3])
      l <- exp(xx$par[2])
      
      lgamma(yNow+1/a)-lgamma(1/a)-lgamma(yNow+1)-(yNow+1/a)*log(1+l*a)+yNow*log(l*a)-log(1-(1+l*a)^(-1/a)-l*(1+l*a)^(-1-1/a))
    }
    
    suppressWarnings(logLikIdeal <- sapply(yUnq, FUN = fff))
    names(logLikIdeal) <- yUnq
    
    logLikIdeal <- sapply(y[iddx], FUN = function(x) {
      logLikIdeal[names(logLikIdeal) == x]
    })
    
    diff <- iddx
    diff[iddx] <- logLikIdeal - logLikFit

    if (any(logLikFit > 0)) {
      warning("Dispertion parameter values are on the boundary of parameter space. Deviance residuals will be asigned 0 on these observations.")
      diff[iddx][logLikFit > 0]   <- 0
    } else if (any(diff < 0)) {
      warning("Numerical deviance finder found worse saturated likelihood than fitted model. Expect NA's in deviance/deviance residuals.")
    }

    diff[diff < 0] <- 0

    sign(y - mu) * sqrt(2 * wt * diff)
  }

  pointEst <- function (pw, eta, contr = FALSE, y, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    N <- pw * (iddx * (1 - lambda * ((1 + lambda * alpha) ^ (-1 - 1 / alpha))) / 
    (1 - (1 + lambda * alpha) ^ (- 1 / alpha) - lambda * 
    ((1 + lambda * alpha) ^ (-1 - 1 / alpha))) + (1 - iddx))
    
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, y, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    prob <- 1 - (1 + lambda * alpha) ^ (-1 / alpha) - 
      lambda * (1 + lambda * alpha) ^ (-1 - 1 / alpha)
    
    bigTheta1 <- pw * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) * 
    as.numeric(((lambda*alpha+1)^(1/alpha)*(lambda^2*alpha^2+2*lambda*alpha+1)*log(lambda*alpha+1)+(lambda*alpha+1)^(1/alpha)*(-lambda^2*alpha^2-lambda*alpha)-lambda^2*alpha^2)/(alpha^2*((lambda*alpha+1)^(1/alpha+1)-lambda*alpha-lambda-1)^2))
    
    bigTheta2 <- pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
    as.numeric(-((alpha*lambda+1)^(1/alpha+1)-1)/((alpha*lambda+1)^(1/alpha+1)+(-alpha-1)*lambda-1)^2)
    
    bigTheta <- t(c(bigTheta2, bigTheta1)[rep(iddx, 2)] %*% Xvlm[rep(iddx, 2), , drop = FALSE])
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum((pw * ((1 - prob) / (prob ^ 2)) * 
    (1 - lambda * ((1 + lambda * alpha) ^ (-1 - 1 / alpha))) ^ 2)[iddx])
    
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
      init <- log(abs((observed / weighted.mean(observed, priorWeights) - 1) / observed) + .1)
      etaStart <- cbind(
        pmin(family$links[[1]](observed), family$links[[1]](12)),
        family$links[[2]](ifelse(init < -.5, .1, init + .55))
      ) + offset
    } else if (method == "optim") {
      init <- c(
        family$links[[1]](weighted.mean(observed, priorWeights)),
        family$links[[2]](abs((cov.wt(cbind(observed, observed), wt = priorWeights, method = "ML")$cov[1,1] / weighted.mean(observed, priorWeights) - 1) / weighted.mean(observed, priorWeights)) + .1)
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

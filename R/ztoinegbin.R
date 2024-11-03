#' @rdname singleRmodels
#' @importFrom stats uniroot
#' @importFrom stats dnbinom
#' @export
ztoinegbin <- function(nSim = 1000, epsSim = 1e-8, eimStep = 6,
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
        "nontrunc" = omega * (1 - (1 + alpha * lambda) ^ (-1 / alpha)) + (1 - omega) * lambda,
        "trunc" = omega + (1 - omega) * lambda / (1 - (1 + alpha * lambda) ^ (-1 / alpha))
      )
    } else {
      switch (
        type,
        "nontrunc" = {
          matrix(c(
            1 + omega * (alpha * lambda + 1) ^ (-1 / alpha - 1) - omega,
            -omega * (log(lambda * alpha + 1) / alpha ^ 2 -
            lambda / (alpha * (lambda * alpha + 1))) / 
            (lambda * alpha + 1) ^ (1 / alpha),
            1 - lambda - 1 / (alpha * lambda + 1) ^ (1 / alpha)
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
            (-alpha - 1) * lambda - 1) / 
            ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2,
            (1 - omega) * lambda * (lambda * alpha + 1) ^ (1 / alpha - 1) *
            ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha) /
            (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2),
            -lambda / (1 - (1 + alpha * lambda) ^ (-1 / alpha))
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
      nontrunc = omega * (1 - (1 + alpha * lambda) ^ (-1 / alpha)) + 
      (1 - omega) * (lambda + (1 + alpha) * lambda ^ 2),
      trunc = omega + (1 - omega) * (lambda + (1 + alpha) * lambda ^ 2) / 
      (1 - (1 + alpha * lambda) ^ (-1 / alpha))
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
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    #P0 <- stats::dnbinom(x = 0, size = 1 / alpha, mu = lambda)
    res <- rep(0, NROW(eta))
    k <- 2
    finished <- rep(FALSE, NROW(eta))
    while ((k < nSim) & !all(finished)) {
      prob <- apply(cbind(k:(k + eimStep)), MARGIN = 1, FUN = function(x) {
        (1 - omega) * stats::dnbinom(
          x = x, 
          size = 1 / alpha, 
          mu = lambda
        ) / (1 - P0)
      })
      trg <- apply(cbind(k:(k + eimStep)), MARGIN = 1, FUN = function(x) {
        comptrigamma(y = x, alpha = alpha)
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
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    
    z  <- omega + (1 - omega) * lambda * 
    ((1 + alpha * lambda) ^ (-1 / alpha - 1)) / 
    (1 - (1 + alpha * lambda) ^ (-1 / alpha))
    
    XXX <- mu.eta(eta, type = "trunc") - z
    
    Etrig <- compExpectG1(eta)
    
    # omega
    G00 <- (-z * ((alpha * lambda + 1) ^ (1 / alpha + 1) +
    (-alpha - 1) * lambda - 1) ^ 2) / (((alpha * lambda + 1) ^ (1 / alpha + 1) +
    (-alpha - 1) * lambda - 1) * omega + lambda) ^ 2 -
    (1 - z) / (1 - omega) ^ 2
    
    G00 <- G00 * omegaLink(eta[, 3], inverse = TRUE, deriv = 1) ^ 2
    
    # omega alpha
    G01 <- (-z) * lambda * ((alpha * lambda + 1) ^ (1 / alpha + 1) *
    log(alpha * lambda + 1) + (-alpha ^ 2 - alpha) * lambda * 
    (alpha * lambda + 1) ^ (1 / alpha) + alpha ^ 2 * lambda) /
    (alpha ^ 2 * (((alpha * lambda + 1) ^ (1 / alpha + 1) +
    (-alpha - 1) * lambda - 1) * omega + lambda) ^ 2)
    
    G01 <- G01 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
                 omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
    
    # omega lambda
    G02 <- z * ((lambda - 1) * (alpha * lambda + 1) ^ (1 / alpha) + 1) /
    (((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) * omega + lambda) ^ 2
    
    G02 <- G02 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                  omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
    
    # alpha
    G11 <- (1 - z) * ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / 
    (alpha * (lambda * alpha + 1))) ^ 2 / ((lambda * alpha + 1) ^ (1 / alpha) *
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) + (log(lambda * alpha + 1) / alpha ^ 2 -
    lambda / (alpha * (lambda * alpha + 1))) ^ 2 / ((lambda * alpha + 1) ^ (2 / alpha) *
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha)) ^ 2) + (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 +
    (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) + lambda ^ 2 / 
    (alpha * (lambda * alpha + 1) ^ 2)) / ((lambda * alpha + 1) ^ (1 / alpha) * 
    (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) - (2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
    2 * lambda / (alpha ^ 2 * (lambda * alpha + 1))) - 
    (lambda ^ 2 * (-(1 - z) / alpha - XXX)) / 
    (lambda * alpha + 1) ^ 2 - XXX / alpha ^ 2 + Etrig +
    (1 - omega) * z * lambda * ((lambda * alpha + 1) ^ (1 / alpha + 1) *
    log(lambda * alpha + 1) * ((lambda * (1 / alpha + 1)) / (lambda * alpha + 1) -
    log(lambda * alpha + 1) / alpha ^ 2) + (lambda * alpha + 1) ^ (1 / alpha) * 
    (-lambda * alpha ^ 2 - lambda * alpha) * (lambda / (alpha * (lambda * alpha + 1)) - 
    log(lambda * alpha + 1) / alpha ^ 2) + (-2 * lambda * alpha - lambda) * 
    (lambda * alpha + 1) ^ (1 / alpha) + lambda * (lambda * alpha + 1) ^ (1 / alpha) +
    2 * lambda * alpha) / (alpha ^ 2 * (lambda * alpha + 1) * 
    ((lambda * alpha + 1) ^ (1 / alpha) * (omega * lambda * alpha + omega) - 
    omega * lambda * alpha + (1 - omega) * lambda - omega) * ((lambda * alpha + 1) ^ (1 / alpha) - 1)) -
    (1 - omega) * z * lambda * ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) +
    (lambda * alpha + 1) ^ (1 / alpha) * (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2) *
    ((lambda * alpha + 1) ^ (1 / alpha) * (omega * lambda * alpha + omega) *
    (lambda / (alpha * (lambda * alpha + 1)) - log(lambda * alpha + 1) / alpha ^ 2) + 
    omega * lambda * (lambda * alpha + 1) ^ (1 / alpha) - omega * lambda) / 
    (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) *
    (omega * lambda * alpha + omega) - omega * lambda * alpha + (1 - omega) * lambda - omega) ^ 2 *
    ((lambda * alpha + 1) ^ (1 / alpha) - 1)) + ((omega - 1) * z * lambda * 
    (lambda * alpha + 1) ^ (1 / alpha - 1) * (lambda / (alpha * (lambda * alpha + 1)) -
    log(lambda * alpha + 1) / alpha ^ 2) * ((lambda * alpha + 1) ^ (1 / alpha + 1) * 
    log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (1 / alpha) *
    (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2)) /
    (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) * (omega * lambda * alpha + omega) -
    omega * lambda * alpha + (1 - omega) * lambda - omega) * ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2) +
    (2 * (omega - 1) * z * lambda * ((lambda * alpha + 1) ^ (1 / alpha + 1) * 
    log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (1 / alpha) * 
    (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2)) /
    (alpha ^ 3 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) *
    (omega * lambda * alpha + omega) - omega * lambda * alpha + (1 - omega) * lambda - omega) *
    ((lambda * alpha + 1) ^ (1 / alpha) - 1)) + ((omega - 1) * z * lambda ^ 2 * 
    ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) + 
    (lambda * alpha + 1) ^ (1 / alpha) * (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2)) /
    (alpha ^ 2 * (lambda * alpha + 1) ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) *
    (omega * lambda * alpha + omega) - omega * lambda * alpha + 
    (1 - omega) * lambda - omega) * ((lambda * alpha + 1) ^ (1 / alpha) - 1))

    G11 <- G11 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2
    
    # alpha lambda
    G12 <- ((omega * (lambda * alpha + 1) ^ (1 + 1 / alpha) -
    lambda * (omega * alpha + omega - 1) - omega) * 
    (((lambda * alpha + 1) ^ (1 / alpha) - 1) * 
    ((omega - 1) * z * (lambda - 1) * (lambda * alpha + 1) ^ (1 / alpha) *
    (lambda * alpha - (lambda * alpha + 1) * log(lambda * alpha + 1)) -
    (omega - 1) * z * lambda * alpha ^ 2 * ((lambda - 1) * 
    (lambda * alpha + 1) ^ (1 / alpha) + 1)) - (omega - 1) * z * 
    (lambda * alpha + 1) ^ (1 / alpha) * ((lambda - 1) * (lambda * alpha + 1) ^ (1 / alpha) + 1) *
    (lambda * alpha - (lambda * alpha + 1) * log(lambda * alpha + 1))) -
    (omega - 1) * z * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) - 1) *
    ((lambda - 1) * (lambda * alpha + 1) ^ (1 / alpha) + 1) * 
    (omega * (lambda * alpha + 1) ^ (1 / alpha) * 
    (lambda * alpha * (alpha + 1) - (lambda * alpha + 1) * log(lambda * alpha + 1)) -
    omega * lambda * alpha ^ 2)) / (alpha ^ 2 * (lambda * alpha + 1) ^ 2 * 
    ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2 *
    (omega * (lambda * alpha + 1) ^ (1 + 1 / alpha) -
    lambda * (omega * alpha + omega - 1) - omega) ^ 2) +
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
    G22 <- (1 - z) * (-((-1 / alpha - 1) * alpha * (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
    (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) + (alpha * lambda + 1) ^ (-2 / alpha - 2) /
    (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2) - XXX / lambda ^ 2 - 
    (alpha ^ 2 * (-XXX - (1 - z) / alpha)) / (alpha * lambda + 1) ^ 2 - 
    ((omega - 1) * z * ((alpha * lambda + 1) ^ (1 / alpha) *
    (((alpha ^ 2 - 1) * omega - alpha + 1) * lambda ^ 2 + 
    ((-6 * alpha ^ 2 - 5 * alpha - 1) * omega + 4 * alpha) * lambda +
    (-6 * alpha - 4) * omega + 2) + (alpha * lambda + 1) ^ (2 / alpha) *
    (((-2 * alpha ^ 2 - alpha) * omega + alpha) * lambda ^ 2 +
    ((6 * alpha ^ 2 + 4 * alpha) * omega - 2 * alpha) * lambda + (6 * alpha + 5) * omega - 1) +
    (alpha * lambda + 1) ^ (3 / alpha) * ((alpha ^ 2 + alpha) * omega * lambda ^ 2 +
    (-2 * alpha ^ 2 - alpha + 1) * omega * lambda + (-2 * alpha - 2) * omega) +
    ((2 * alpha ^ 2 + 2 * alpha) * omega - 2 * alpha) * lambda + (2 * alpha + 1) * omega - 1)) /
    ((alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) * (alpha * omega * lambda + omega) +
    ((-alpha - 1) * omega + 1) * lambda - omega) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)
    
    
    G22 <- G22 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
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
  
  funcZ <- function(eta, weight, y, prior, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    z <- as.numeric(y == 1)
    weight <- weight / prior
    
    dig <-  compdigamma(y = y, alpha = alpha)
    
    G0 <- z * ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) *
    lambda - 1) / (((alpha * lambda + 1) ^ (1 / alpha + 1) +
    (-alpha - 1) * lambda - 1) * omega + lambda) -
    (1 - z) / (1 - omega)
    
    G1 <- ((1 - omega) * z * lambda * ((lambda * alpha + 1) ^ (1 / alpha + 1) *
    log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (1 / alpha) *
    (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2)) /
    (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) *
    (omega * lambda * alpha + omega) - omega * lambda * alpha + 
    (1 - omega) * lambda - omega) * ((lambda * alpha + 1) ^ (1 / alpha) - 1)) +
    
      (1 - z) * ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
    ((lambda * alpha + 1) ^ (1 / alpha) * (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) +
    log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) /
    (lambda * alpha + 1) + y / alpha + dig)
    
    G2 <- (1 - z) * (-(alpha * lambda + 1) ^ (-1 / alpha - 1) / 
    (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) + (alpha * (-y - 1 / alpha)) /
    (alpha * lambda + 1) + y / lambda) + 
    ((omega - 1) * z * ((lambda - 1) * (alpha * lambda + 1) ^ (1 / alpha) + 1)) /
    ((alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) * 
    (alpha * omega * lambda + omega) + 
    ((-alpha - 1) * omega + 1) * lambda - omega) * 
    ((alpha * lambda + 1) ^ (1 / alpha) - 1))
    
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
    y <- as.numeric(y)
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X) / 3), rep(0, NROW(X) / 3), rep(0, NROW(X) / 3))
    }
    
    z <- as.numeric(y == 1)
    X <- as.matrix(X)
    
    if (!(deriv %in% c(0, 1, 2))) 
      stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
            function(beta) {
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3) + offset
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              omega  <-  omegaLink(eta[, 3], inverse = TRUE)
              
              -sum(weight * (z * log(omega + (1 - omega) *
              lambda * (1 + alpha * lambda) ^ (-1 / alpha - 1) / 
              (1 - (1 + lambda * alpha) ^ (-1 / alpha))) + 
              (1 - z) * (log(1 - omega) + lgamma(y + 1 / alpha) - 
              lgamma(1 / alpha) - lgamma(y + 1) - 
              (y + 1 / alpha) * log(1 + lambda * alpha) + y * log(lambda * alpha) - 
              log(1 - (1 + lambda * alpha) ^ (-1 / alpha)))))
            },
            function(beta) {
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3) + offset
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              omega  <-  omegaLink(eta[, 3], inverse = TRUE)
              
              dig <-  compdigamma(y = y, alpha = alpha)
              
              G0 <- z * ((alpha * lambda+1) ^ (1 / alpha + 1) + (-alpha - 1) *
              lambda - 1) / (((alpha * lambda + 1) ^ (1 / alpha + 1) +
              (-alpha - 1) * lambda - 1) * omega + lambda) -
              (1 - z) / (1 - omega)
              
              G1 <- ((1 - omega) * z * lambda * ((lambda * alpha + 1) ^ (1 / alpha + 1) *
              log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (1 / alpha) *
              (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2)) /
              (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) *
              (omega * lambda * alpha + omega) - omega * lambda * alpha + 
              (1 - omega) * lambda - omega) * ((lambda * alpha + 1) ^ (1 / alpha) - 1)) +
              (1 - z) * ((log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
              ((lambda * alpha + 1) ^ (1 / alpha) * (1 - 1 / (lambda * alpha + 1) ^ (1 / alpha))) +
              log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) /
              (lambda * alpha + 1) + y / alpha + dig)
              
              G2 <- (1 - z) * (-(alpha * lambda + 1) ^ (-1 / alpha - 1) / 
              (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) + (alpha * (-y - 1 / alpha)) /
              (alpha * lambda + 1) + y / lambda) + 
              ((omega - 1) * z * ((lambda - 1) * (alpha * lambda + 1) ^ (1 / alpha) + 1)) /
              ((alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) * 
              (alpha * omega * lambda + omega) + 
              ((-alpha - 1) * omega + 1) * lambda - omega) * 
              ((alpha * lambda + 1) ^ (1 / alpha) - 1))
              
              
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
              eta <- matrix(as.matrix(X) %*% beta, ncol = 3) + offset
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              omega  <-  omegaLink(eta[, 3], inverse = TRUE)
              res <- matrix(nrow = length(beta), 
                            ncol = length(beta), 
                            dimnames = list(names(beta), names(beta)))
              
              trig <- comptrigamma(y = y, alpha = alpha)
              dig <-  compdigamma(y = y, alpha = alpha)
              
              # omega
              G0 <- z * ((alpha * lambda+1) ^ (1 / alpha + 1) + (-alpha - 1) *
              lambda - 1) / (((alpha * lambda + 1) ^ (1 / alpha + 1) +
              (-alpha - 1) * lambda - 1) * omega + lambda) -
              (1 - z) / (1 - omega)
              
              G00 <- (-z * ((alpha * lambda + 1) ^ (1 / alpha + 1) +
              (-alpha - 1) * lambda - 1) ^ 2) / (((alpha * lambda + 1) ^ (1 / alpha + 1) +
              (-alpha - 1) * lambda - 1) * omega + lambda) ^ 2 -
              (1 - z) / (1 - omega) ^ 2
              
              G00 <- G00 * omegaLink(eta[, 3], inverse = TRUE, deriv = 1) ^ 2 +
                      G0 * omegaLink(eta[, 3], inverse = TRUE, deriv = 2)
              
              # omega alpha
              G01 <- (-z) * lambda * ((alpha * lambda + 1) ^ (1 / alpha + 1) *
              log(alpha * lambda + 1) + (-alpha ^ 2 - alpha) * lambda * 
              (alpha * lambda + 1) ^ (1 / alpha) + alpha ^ 2 * lambda) /
              (alpha ^ 2 * (((alpha * lambda + 1) ^ (1 / alpha + 1) +
              (-alpha - 1) * lambda - 1) * omega + lambda) ^ 2)
              
              G01 <- G01 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
                           omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
              
              # omega lambda
              G02 <- z * ((lambda - 1) * (alpha * lambda + 1) ^ (1 / alpha) + 1) /
              (((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) * omega + lambda) ^ 2
              
              G02 <- G02 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                            omegaLink(eta[, 3], inverse = TRUE, deriv = 1)
              # alpha
              G1 <- ((1 - omega) * z * lambda * ((lambda * alpha + 1) ^ (1 / alpha + 1) *
              log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (1 / alpha) *
              (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2)) /
              (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) *
              (omega * lambda * alpha + omega) - omega * lambda * alpha + 
              (1 - omega) * lambda - omega) * ((lambda * alpha + 1) ^ (1 / alpha) - 1)) +
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
              (lambda * alpha + 1) ^ 2 - y / alpha ^ 2 + trig) - 
              ((omega - 1) * z * lambda * ((lambda * alpha + 1) ^ (1 / alpha + 1) *
              log(lambda * alpha + 1) * ((lambda * (1 / alpha + 1)) / (lambda * alpha + 1) -
              log(lambda * alpha + 1) / alpha ^ 2) + (lambda * alpha + 1) ^ (1 / alpha) * 
              (-lambda * alpha ^ 2 - lambda * alpha) * (lambda / (alpha * (lambda * alpha + 1)) - 
              log(lambda * alpha + 1) / alpha ^ 2) + (-2 * lambda * alpha - lambda) * 
              (lambda * alpha + 1) ^ (1 / alpha) + lambda * (lambda * alpha + 1) ^ (1 / alpha) +
              2 * lambda * alpha)) / (alpha ^ 2 * (lambda * alpha + 1) * 
              ((lambda * alpha + 1) ^ (1 / alpha) * (omega * lambda * alpha + omega) - 
              omega * lambda * alpha + (1 - omega) * lambda - omega) * ((lambda * alpha + 1) ^ (1 / alpha) - 1)) +
              ((omega - 1) * z * lambda * ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) +
              (lambda * alpha + 1) ^ (1 / alpha) * (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2) *
              ((lambda * alpha + 1) ^ (1 / alpha) * (omega * lambda * alpha + omega) *
              (lambda / (alpha * (lambda * alpha + 1)) - log(lambda * alpha + 1) / alpha ^ 2) + 
              omega * lambda * (lambda * alpha + 1) ^ (1 / alpha) - omega * lambda)) / 
              (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) *
              (omega * lambda * alpha + omega) - omega * lambda * alpha + (1 - omega) * lambda - omega) ^ 2 *
              ((lambda * alpha + 1) ^ (1 / alpha) - 1)) + ((omega - 1) * z * lambda * 
              (lambda * alpha + 1) ^ (1 / alpha - 1) * (lambda / (alpha * (lambda * alpha + 1)) -
              log(lambda * alpha + 1) / alpha ^ 2) * ((lambda * alpha + 1) ^ (1 / alpha + 1) * 
              log(lambda * alpha + 1) + (lambda * alpha+ 1 ) ^ (1 / alpha) *
              (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2)) /
              (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) * (omega * lambda * alpha + omega) -
              omega * lambda * alpha + (1 - omega) * lambda - omega) * ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2) +
              (2 * (omega - 1) * z * lambda * ((lambda * alpha + 1) ^ (1 / alpha + 1) * 
              log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (1 / alpha) * 
              (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2)) /
              (alpha ^ 3 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) *
              (omega * lambda * alpha + omega) - omega * lambda * alpha + (1 - omega) * lambda - omega) *
              ((lambda * alpha + 1) ^ (1 / alpha) - 1)) + ((omega - 1) * z * lambda ^ 2 * 
              ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) + 
              (lambda * alpha + 1) ^ (1 / alpha) * (-lambda * alpha ^ 2 - lambda * alpha) + lambda * alpha ^ 2)) /
              (alpha ^ 2 * (lambda * alpha + 1) ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) *
              (omega * lambda * alpha + omega) - omega * lambda * alpha + 
              (1 - omega) * lambda - omega) * ((lambda * alpha + 1) ^ (1 / alpha) - 1))
              
              G11 <- G11 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
                      G1 * alphaLink(eta[, 2], inverse = TRUE, deriv = 2)
              
              # alpha lambda
              G12 <- ((omega * (lambda * alpha + 1) ^ (1 + 1 / alpha) -
              lambda * (omega * alpha + omega - 1) - omega) * 
              (((lambda * alpha + 1) ^ (1 / alpha) - 1) * 
              ((omega - 1) * z * (lambda - 1) * (lambda * alpha + 1) ^ (1 / alpha) *
              (lambda * alpha - (lambda * alpha + 1) * log(lambda * alpha + 1)) -
              (omega - 1) * z * lambda * alpha ^ 2 * ((lambda - 1) * 
              (lambda * alpha + 1) ^ (1 / alpha) + 1)) - (omega - 1) * z * 
              (lambda * alpha + 1) ^ (1 / alpha) * ((lambda - 1) * (lambda * alpha + 1) ^ (1 / alpha) + 1) *
              (lambda * alpha - (lambda * alpha + 1) * log(lambda * alpha + 1))) -
              (omega - 1) * z * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) - 1) *
              ((lambda - 1) * (lambda * alpha + 1) ^ (1 / alpha) + 1) * 
              (omega * (lambda * alpha + 1) ^ (1 / alpha) * 
              (lambda * alpha * (alpha + 1) - (lambda * alpha + 1) * log(lambda * alpha + 1)) -
              omega * lambda * alpha ^ 2)) / (alpha ^ 2 * (lambda * alpha + 1) ^ 2 * 
              ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2 *
              (omega * (lambda * alpha + 1) ^ (1 + 1 / alpha) -
              lambda * (omega * alpha + omega - 1) - omega) ^ 2) +
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
              ((omega - 1) * z * ((lambda - 1) * (alpha * lambda + 1) ^ (1 / alpha) + 1)) /
              ((alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) * 
              (alpha * omega * lambda + omega) + 
              ((-alpha - 1) * omega + 1) * lambda - omega) * 
              ((alpha * lambda + 1) ^ (1 / alpha) - 1))
              
              G22 <- (1 - z) * (-((-1 / alpha - 1) * alpha * (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
              (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) + (alpha * lambda + 1) ^ (-2 / alpha - 2) /
              (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) ^ 2 - (alpha ^ 2 * (-y - 1 / alpha)) /
              (alpha * lambda + 1) ^ 2 - y / lambda ^ 2) - ((omega - 1) * z * 
              ((alpha * lambda + 1) ^ (1 / alpha) * (((alpha ^ 2 - 1) * omega - alpha + 1) * lambda ^ 2 +
              ((-6 * alpha ^ 2 - 5 * alpha - 1) * omega + 4 * alpha) * lambda +
              (-6 * alpha - 4) * omega + 2) + (alpha * lambda + 1) ^ (2 / alpha) *
              (((-2 * alpha ^ 2 - alpha) * omega + alpha) * lambda ^ 2 + 
              ((6 * alpha ^ 2 + 4 * alpha) * omega - 2 * alpha) * lambda +
              (6 * alpha + 5) * omega - 1) + (alpha * lambda + 1) ^ (3 / alpha) *
              ((alpha ^ 2 + alpha) * omega * lambda ^ 2 + 
              (-2 * alpha ^ 2 - alpha + 1) * omega * lambda + (-2 * alpha - 2) * omega) +
              ((2 * alpha ^ 2 + 2 * alpha) * omega - 2 * alpha) * lambda +
              (2 * alpha + 1) * omega - 1)) / ((alpha * lambda + 1) ^ 2 *
              ((alpha * lambda + 1) ^ (1 / alpha) * (alpha * omega * lambda + omega) +
              ((-alpha - 1) * omega + 1) * lambda - omega) ^ 2 * 
              ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)
               
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
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    omega  <-  omegaLink(eta[, 3], inverse = TRUE)
    z <- (y == 1)
    mu <- mu.eta(eta = eta)
    
    logLikFit <- (
      z * log(omega + (1 - omega) * lambda * (1 + alpha * lambda) ^ (-1 / alpha - 1) / 
      (1 - (1 + lambda * alpha) ^ (-1 / alpha))) + (1 - z) * (log(1 - omega) + lgamma(y + 1 / alpha) - 
      lgamma(1 / alpha) - lgamma(y + 1) - (y + 1 / alpha) * log(1 + lambda * alpha) + y * log(lambda * alpha) - 
      log(1 - (1 + lambda * alpha) ^ (-1 / alpha)))
    )
    
    yUnq <- unique(y)
    
    if (any(yUnq > 77)) {
      warning("Curently numerical deviance is unreliable for counts greater than 78.")
    }
    
    findL <- function(t) {
      yNow <- yUnq[t]
      stats::optim(
        #par = if(yNow < 26) c(0, .6, 0) else c(-.5, log(yNow), -20),
        par = c(0, log(yNow), -10),
        fn = function(x) {
          s <- x[1]
          l <- exp(x[2])
          a <- exp(x[3])
          
          prob <- 1 - (1+a*l)^(-1/a)
          prob <- 1 / prob
          sum(c((l*prob - yNow) * 4.5,# s der
                yNow/l+(-yNow*a-1)/(1+a*l)-(1+a*l)^(-1-1/a)*prob+s*(prob-prob^2*(l*(1+a*l)^(-1-1/a))),# lambda der
                (log(l*a+1)/a^2-l/(a*(l*a+1)))/((l*a+1)^(1/a)*(1-1/(l*a+1)^(1/a)))+(s*l*(log(l*a+1)/a^2-l/(a*(l*a+1))))/((l*a+1)^(1/a)*(1-1/(l*a+1)^(1/a))^2)+log(l*a+1)/a^2+(l*(-1/a-yNow))/(l*a+1)+yNow/a-digamma(yNow+1/a)/a^2+digamma(1/a)/a^2,#alpha der
                #this is experimental
                lgamma(yNow+1/a)-lgamma(1/a) - lgamma(yNow+1)-(yNow+1/a)*log(1+a*l)+yNow*log(l*a)-log(1-(1+a*l)^(-1/a))) ^ 2) ^ .5
        },
        method = "BFGS",
        control = list(maxit = 10000, abstol = .Machine$double.eps, reltol = .Machine$double.eps)
      )$par
    }
    
    suppressWarnings({
      logLikIdeal <- sapply(1:length(yUnq), FUN = function(x) {
        ifelse(yUnq[x] == 1, 0, {
          xx <- findL(x)
          lagrange <- xx[1]
          l <- exp(xx[2])
          a <- exp(xx[3])
          (lgamma(yUnq[x] + 1 / a) - lgamma(1 / a) -
              lgamma(yUnq[x] + 1) - (yUnq[x] + 1 / a) * log(1 + a * l) +
              yUnq[x] * log(l * a) - log(1 - (1 + a * l) ^ (-1 / a)))
        })
      })
    })
    
    logLikIdeal <- sapply(1:length(y), FUN = function(x) {
      logLikIdeal[yUnq == y[x]]
    })
    
    diff <- logLikIdeal - logLikFit
    
    if (any(logLikFit > 0)) {
      warning("Dispertion parameter values are on the boundary of parameter space. Deviance residuals will be asigned 0 on these observations.")
      diff[logLikFit > 0]   <- 0
    } else if (any(diff < 0)) {
      warning("Numerical deviance finder found worse saturated likelihood than fitted model. Expect NA's in deviance/deviance residuals.")
      diff[diff < 0] <- 0
    }
    
    #diff <- ifelse(abs(diff) < 1e-1 & diff > 0, 0, diff)
    
    sign(y - mu) * sqrt(2 * wt * diff)
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
    ((lambda * alpha + 1) ^ (1 / alpha - 1) * 
    ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha)) /
    (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)# w.r to alpha
    bigTheta2 <- -pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
    ((alpha * lambda + 1) ^ (1 / alpha - 1) / ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)# w.r to lambda
    
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
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    
    switch (type,
      "trunc" = {
        (1 - omega) * stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / 
        (1 - P0) + omega * as.numeric(x == 1)
      },
      "nontrunc" = {
        stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) * 
        (as.numeric(x == 0) + as.numeric(x > 0) * (1 - omega)) +
        omega * (1 - P0) * as.numeric(x == 1)
      }
    )
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
    cond <- CDF(sims) < p_u
    while (any(cond)) {
      sims[cond] <- sims[cond] + 1
      cond <- CDF(sims) < p_u
    }
    sims
  }
  
  # new starting points
  
  getStart <- expression(
    if (method == "IRLS") {
      init <- log(abs((observed / weighted.mean(observed, priorWeights) - 1) / observed) + .1)
      etaStart <- cbind(
        pmin(family$links[[1]](observed), family$links[[1]](12)),
        family$links[[2]](ifelse(init < -.5, .1, init + .55)),
        family$links[[3]](weighted.mean(observed == 1, priorWeights) * (.5 + .5 * (observed == 1)) + .01)
      ) + offset
      # print(summary(etaStart))
      # print(summary(cbind(family$links[[1]](etaStart[,1], inverse = TRUE),family$links[[2]](etaStart[,2], inverse = TRUE),family$links[[3]](etaStart[,3], inverse = TRUE))))
      # stop("abc")
    } else if (method == "optim") {
      init <- c(
        family$links[[1]](weighted.mean(observed, priorWeights)),
        family$links[[2]](abs((cov.wt(cbind(observed, observed), wt = priorWeights, method = "ML")$cov[1,1] / weighted.mean(observed, priorWeights) - 1) / weighted.mean(observed, priorWeights)) + .1),
        family$links[[3]](weighted.mean(observed == 1, priorWeights) + .01)
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
      if ("(Intercept):omega" %in% colnames(Xvlm)) {
        coefStart <- c(coefStart, init[3], rep(0, attr(Xvlm, "hwm")[3] - 1))
      } else {
        coefStart <- c(coefStart, rep(init[3] / attr(Xvlm, "hwm")[3], attr(Xvlm, "hwm")[3]))
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
      family    = "ztoinegbin",
      etaNames  = c("lambda", "alpha", "omega"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

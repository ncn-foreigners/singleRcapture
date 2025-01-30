#' @rdname singleRmodels
#' @importFrom stats uniroot
#' @importFrom stats dnbinom
#' @export
ztHurdlenegbin <- function(nSim = 1000, epsSim = 1e-8, eimStep = 6,
                           lambdaLink = c("log", "neglog"), 
                           alphaLink = c("log", "neglog"),
                           piLink = c("logit", "cloglog", "probit"), 
                           ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  if (missing(alphaLink))  alphaLink  <- "log"
  if (missing(piLink))     piLink     <- "logit"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink, alphaLink, piLink)
  
  lambdaLink <- switch(lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  alphaLink <- switch(alphaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  piLink <- switch(piLink,
    "logit" = singleRinternallogitLink,
    "cloglog" = singleRinternalcloglogLink,
    "probit" = singleRinternalprobitLink
  )
  
  links[1:3] <- c(lambdaLink, alphaLink, piLink)
  
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    PI     <-     piLink(eta[, 3], inverse = TRUE)
    
    if (!deriv) {
      switch (type,
        "nontrunc" = (1 - (1 + alpha * lambda) ^ (-1 / alpha)) * (PI + (1 - PI) * lambda),
        "trunc" = PI + (1 - PI) * 
        (lambda - lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha))) / 
        (1 - (1 + alpha * lambda) ^ (-1 / alpha) - 
        lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha)))
      )
    } else {
      switch (
        type,
        "nontrunc" = {
          matrix(c(
            (1 - PI) * (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha)) +
            ((1 - PI) * lambda + PI) * (alpha * lambda + 1) ^ (-1 / alpha - 1),
            -((1 - PI) * lambda + PI) * 
            (log(lambda * alpha + 1) / alpha ^ 2 - 
            lambda / (alpha * (lambda * alpha + 1))) /
            (lambda * alpha + 1) ^ (1 / alpha),
            (1 - lambda) * (1 - 1 / (alpha * lambda + 1) ^ (1 / alpha))
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
              alphaLink(eta[,2], inverse = TRUE, deriv = 1),
                piLink(eta[, 3], inverse = TRUE, deriv = 1)
          ), ncol = 3)
        },
        "trunc" = {
          matrix(c(
            (1 - PI) * ((alpha * lambda + 1) ^ (2 / alpha) * (alpha ^ 2 * lambda ^ 2 + 2 * alpha * lambda + 1) +
            (alpha * lambda + 1) ^ (1 / alpha) * ((-alpha ^ 2 - 2 * alpha - 1) * lambda ^ 2 - 2 * alpha * lambda - 2) + 1) /
            ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2,
            (1 - PI) * lambda ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) * 
            (lambda * alpha ^ 2 + (lambda + 1) * alpha + 1) * log(lambda * alpha + 1) +
            (lambda * alpha + 1) ^ (1 / alpha) * ((1 - 2 * lambda) * alpha ^ 2 - lambda * alpha) - alpha ^ 2) /
            (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha + 1) - lambda * alpha - lambda - 1) ^ 2),
            1 - (lambda - lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha))) / 
            (1 - (1 + alpha * lambda) ^ (-1 / alpha) - 
            lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha)))
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
             alphaLink(eta[, 2], inverse = TRUE, deriv = 1),
                piLink(eta[, 3], inverse = TRUE, deriv = 1)
          ), ncol = 3)
        }
      )
    }
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    PI     <-     piLink(eta[, 3], inverse = TRUE)
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    
    switch (type,
      nontrunc = (1 - P0) * (PI + (1 - PI) * lambda + (1 + alpha) * lambda ^ 2),
      trunc = PI + (1 - PI) * (lambda * (1 + alpha * lambda) + lambda ^ 2 - 
      lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha))) / 
      (1 - (1 + alpha * lambda) ^ (-1 / alpha) - 
      lambda * ((1 + alpha * lambda) ^ (-1 - 1 / alpha)))
    ) - mu.eta(eta = eta, type = type) ^ 2
  }
  
  compdigamma <- function(y, alpha) {
    (-digamma(y + 1 / alpha) + digamma(1 / alpha)) / (alpha ^ 2)
  }
  
  
  comptrigamma <- function(y, alpha) {
    (2 * (digamma(y + 1 / alpha) - digamma(1 / alpha)) * alpha +
       trigamma(y + 1 / alpha) - trigamma(1 / alpha)) / (alpha ^ 4)
  }
  
  # Computing the expected value of di/trigamma functions on (y + 1/alpha)
  
  compExpectG1 <- function(eta) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    PI     <-  piLink(eta[, 3], inverse = TRUE)
    
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    P1 <- lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1)
    res <- rep(0, NROW(eta))
    
    # 1 is the first possible y value for 0 truncated hurdle distribution
    # but here we compute the (1 - z) * psi function which takes 0 at y = 1
    k <- 2 
    finished <- rep(FALSE, NROW(eta))
    while ((k < nSim) & !all(finished)) {
      prob <- apply(cbind(k:(k + eimStep)), MARGIN = 1, FUN = function(x) {
        (1 - PI) * stats::dnbinom(
          x = x, 
          size = 1 / alpha, 
          mu = lambda
        ) / (1 - P0 - P1)
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
    PI     <-     piLink(eta[, 3], inverse = TRUE)
    
    z  <- PI
    
    ## expected for (1-I)Y
    XXX <- mu.eta(eta, type = "trunc") - z
    
    Etrig <- compExpectG1(eta)
    
    # PI
    G00 <- (-z / PI ^ 2 - (1 - z) / (1 - PI) ^ 2)
    
    G00 <- G00 * piLink(eta[, 3], inverse = TRUE, deriv = 1) ^ 2
    
    G01 <- rep(0, NROW(eta))
    G02 <- rep(0, NROW(eta))
    
    # alpha
    G11 <- Etrig + (1 - z) * (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
    (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) ^ 2 / 
    (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) ^ 2 -
    (1 - z) * (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) ^ 2 / 
    (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
    (lambda * alpha + 1)) ^ 2 - (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
    (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) + lambda ^ 2 / 
    (alpha * (lambda * alpha + 1) ^ 2)) / (lambda * alpha + 1) ^ (1 / alpha) -
    lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * 
    (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) -
    (lambda ^ 2 * (-1 / alpha - 1)) / (lambda * alpha + 1) ^ 2)) / 
    (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) -
    (1 - z) * (2 * log(lambda * alpha + 1)) / alpha ^ 3 + (1 - z) * (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) -
    (lambda ^ 2 * (-(1 - z) / alpha - XXX)) / (lambda * alpha + 1) ^ 2 - XXX / alpha ^ 2
    
    G11 <- G11 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2
    
    # alpha lambda
    G12 <- (-1) * ((1 - z) * (alpha * lambda + 1) ^ (1 / alpha) * 
    ((alpha ^ 3 + alpha ^ 2) * lambda ^ 3 + (2 * alpha ^ 2 + 2 * alpha) * lambda ^ 2 + (alpha + 1) * lambda) * 
    log(alpha * lambda + 1) + (alpha * lambda + 1) ^ (1 / alpha) * 
    ((1 - z) * (alpha ^ 4 - alpha ^ 3 - alpha ^ 2) * lambda ^ 3 + 
    ((-2 * alpha ^ 4 - 2 * alpha ^ 3) * XXX + (1 - z) * 4 * alpha ^ 3 - 
    (1 - z) *  alpha ^ 2 - (1 - z) *  alpha) * lambda ^ 2 + 
    ((-4 * alpha ^ 3 - 2 * alpha ^ 2) * XXX + (1 - z) *  3 * alpha ^ 2) * lambda - 
    2 * alpha ^ 2 * XXX) + (alpha * lambda + 1) ^ (2 / alpha) * (-(1 - z) * alpha ^ 4 * lambda ^ 3 + 
    (alpha ^ 4 * XXX - (1 - z) * 2 * alpha ^ 3) * lambda ^ 2 + (2 * alpha ^ 3 * XXX - 
    (1 - z) * alpha ^ 2) * lambda + alpha ^ 2 * XXX) + ((alpha ^ 4 + 2 * alpha ^ 3 + alpha ^ 2) * XXX - 
    (1 - z) * 2 * alpha ^ 3 - (1 - z) * alpha ^ 2) * lambda ^ 2 +
    ((2 * alpha ^ 3 + 2 * alpha ^ 2) * XXX - (1 - z) * 2 * alpha ^ 2) * lambda + alpha ^ 2 * XXX) /
    (alpha ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2)
    
    G12 <- G12 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                  alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    #lambda
    G22 <- ((alpha * lambda + 1) ^ (2 / alpha) * ((1 - z) * alpha ^ 3 * lambda ^ 4 + 
    ((1 - z) * 2 * alpha ^ 2 - 2 * alpha ^ 3 * XXX) * lambda ^ 3 +
    ((1 - z) * alpha - 5 * alpha ^ 2 * XXX) * lambda ^ 2 - 4 * alpha * XXX * lambda - XXX) + 
    (alpha * lambda + 1) ^ (1 / alpha) * ((1 - z) * (alpha - alpha ^ 3) * lambda ^ 4 +
    ((4 * alpha ^ 3 + 4 * alpha ^ 2) * XXX - (1 - z) * 4 * alpha ^ 2 - (1 - z) * alpha + (1 - z)) * lambda ^ 3 +
    ((10 * alpha ^ 2 + 6 * alpha) * XXX - (1 - z) * 3 * alpha - (1 - z)) * lambda ^ 2 + 
    (8 * alpha + 2) * XXX * lambda + 2 * XXX) + 
    ((-2 * alpha ^ 3 - 4 * alpha ^ 2 - 2 * alpha) * XXX + (1 - z) * 2 * alpha ^ 2 + 
    (1 - z) *  2 * alpha) * lambda ^ 3 + ((-5 * alpha ^ 2 - 6 * alpha - 1) * XXX + 
    (1 - z) * 2 * alpha + (1 - z)) * lambda ^ 2 + (-4 * alpha - 2) * XXX * lambda - XXX) / 
    (lambda ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2)
    
    G22 <- G22 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
    matrix(
      -c(G22 * prior, G12 * prior, G02 * prior,
         G12 * prior, G11 * prior, G01 * prior,
         G02 * prior, G01 * prior, G00 * prior),
      dimnames = list(
        rownames(eta), 
        c("lambda", "lambda:alpha", "lambda:PI",
          "lambda:alpha", "alpha", "alpha:PI",
          "lambda:PI", "alpha:PI", "PI")
      ),
      ncol = 9
    )
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    PI     <-     piLink(eta[, 3], inverse = TRUE)
    z <- as.numeric(y == 1)
    weight <- weight / prior
    
    dig <- compdigamma(y = y, alpha = alpha)
    
    G0 <- z / PI - (1 - z) / (1 - PI)
    
    G1 <- (1 - z) * (y / alpha + (digamma(1 / alpha) - digamma(1 / alpha + y)) / alpha ^ 2  +
    log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) / (lambda * alpha + 1) -
    (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
    (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
    (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
    (-(lambda * alpha + 1) ^ (-1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1))
    
    G2 <- (1 - z) * (y / lambda  + (alpha * (-y - 1 / alpha)) / (alpha * lambda + 1) -
    ((1 / alpha + 1) * alpha * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
    (-(alpha * lambda + 1) ^ (-1 / alpha) - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 1))
    
    G0 <- G0 *  piLink(eta[, 3], inverse = TRUE, deriv = 1)
    G1 <- G1 *  alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    G2 <- G2 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    uMatrix <- cbind(G2, G1, G0)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 3)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      xx <- solve(weight[[x]])
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
      offset <- cbind(rep(0, NROW(X) / 3), rep(0, NROW(X) / 3), rep(0, NROW(X) / 3))
    }
    
    y <- as.numeric(y)
    z <- as.numeric(y == 1)
    X <- as.matrix(X)
    
    if (!(deriv %in% c(0, 1, 2))) 
      stop("Only score function and derivatives up to 2 are supported.")
    
    # to make it conform to how switch in R works, i.e. indexing begins with 1
    deriv <- deriv + 1
    
    switch (deriv,
            function(beta) {
              eta    <- matrix(as.matrix(X) %*% beta, ncol = 3) + offset
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              PI     <-     piLink(eta[, 3], inverse = TRUE)
              
              -sum(weight * (z * log(PI) + (1 - z) * log(1 - PI) + (1 - z) *
              (lgamma(y + 1 / alpha) - lgamma(1 / alpha) -
              lgamma(y + 1) - (y + 1 / alpha) * log(1 + lambda * alpha) +
              y * log(lambda * alpha) - log(1 - (1 + lambda * alpha) ^ (-1 / alpha) - 
              lambda * (1 + lambda * alpha) ^ (-1 - 1 / alpha)))))
            },
            function(beta) {
              eta    <- matrix(as.matrix(X) %*% beta, ncol = 3) + offset
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              PI     <-     piLink(eta[, 3], inverse = TRUE)
              
              dig <- compdigamma(y = y, alpha = alpha)
              
              G0 <- z / PI - (1 - z) / (1 - PI)
              
              G1 <- (1 - z) * (y / alpha + (digamma(1 / alpha) - digamma(1 / alpha + y)) / alpha ^ 2  +
              log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) / (lambda * alpha + 1) -
              (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
              (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
              (-(lambda * alpha + 1) ^ (-1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1))
              
              G2 <- (1 - z) * (y / lambda  + (alpha * (-y - 1 / alpha)) / (alpha * lambda + 1) -
              ((1 / alpha + 1) * alpha * lambda * (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
              (-(alpha * lambda + 1) ^ (-1 / alpha) - lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 1))
              
              G0 <- G0 * weight *     piLink(eta[, 3], inverse = TRUE, deriv = 1)
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
              eta    <- matrix(as.matrix(X) %*% beta, ncol = 3) + offset
              lambda <- lambdaLink(eta[, 1], inverse = TRUE)
              alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
              PI     <-     piLink(eta[, 3], inverse = TRUE)
              
              res <- matrix(
                nrow = length(beta), 
                ncol = length(beta), 
                dimnames = list(names(beta), names(beta))
              )
              
              trig <- comptrigamma(y = y, alpha = alpha)
              dig  <-  compdigamma(y = y, alpha = alpha)
              
              # PI
              G0 <- z / PI - (1 - z) / (1 - PI)
              
              G00 <- (-z / PI ^ 2 - (1 - z) / (1 - PI) ^ 2)
              
              G00 <- G00 * piLink(eta[, 3], inverse = TRUE, deriv = 1) ^ 2 +
                      G0 * piLink(eta[, 3], inverse = TRUE, deriv = 2)
              
              # alpha
              G1 <- (1 - z) * (y / alpha + (digamma(1 / alpha) - digamma(1 / alpha + y)) / alpha ^ 2 -
              (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) /
              (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) /
              (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) +
              log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - y)) / (lambda * alpha + 1))
              
              G11 <- (1 - z) * ((2 * digamma(1 / alpha + y)) / alpha ^ 3 - 
              (2 * digamma(1 / alpha)) / alpha ^ 3 + trigamma(1 / alpha + y) / alpha ^ 4 - 
              trigamma(1 / alpha) / alpha ^ 4 + (-(log(lambda * alpha + 1) / alpha ^ 2 - 
              lambda / (alpha * (lambda * alpha + 1))) /(lambda * alpha + 1) ^ (1 / alpha) - 
              lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * (log(lambda * alpha + 1) / alpha ^ 2 + 
              (lambda * (-1 / alpha - 1)) / (lambda * alpha + 1))) ^ 2 / 
              (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) ^ 2 -
              (-(log(lambda * alpha + 1) / alpha ^ 2 - lambda / (alpha * (lambda * alpha + 1))) ^ 2 / 
              (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) *
              (log(lambda * alpha + 1) / alpha ^ 2 + (lambda * (-1 / alpha - 1)) /
              (lambda * alpha + 1)) ^ 2 - (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
              (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) + lambda ^ 2 / 
              (alpha * (lambda * alpha + 1) ^ 2)) / (lambda * alpha + 1) ^ (1 / alpha) -
              lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) * (-(2 * log(lambda * alpha + 1)) / alpha ^ 3 + 
              (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) - (lambda ^ 2 * (-1 / alpha - 1)) / (lambda * alpha + 1) ^ 2)) / 
              (-1 / (lambda * alpha + 1) ^ (1 / alpha) - lambda * (lambda * alpha + 1) ^ (-1 / alpha - 1) + 1) -
              (2 * log(lambda * alpha + 1)) / alpha ^ 3 + (2 * lambda) / (alpha ^ 2 * (lambda * alpha + 1)) -
              (lambda ^ 2 * (-1 / alpha - y)) / (lambda * alpha + 1) ^ 2 - y / alpha ^ 2)
              
              G11 <- G11 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
                      G1 * alphaLink(eta[, 2], inverse = TRUE, deriv = 2)
              
              # alpha lambda
              G12 <- (z - 1) * ((alpha * lambda + 1) ^ (1 / alpha) * ((alpha ^ 3 + alpha ^ 2) * lambda ^ 3 +
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
              
              G12 <- G12 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
                            alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
              
              #lambda
              G2 <- (1 - z) * (y / lambda + ((-1 / alpha-1) * alpha * lambda * 
              (alpha * lambda + 1) ^ (-1 / alpha - 2)) / 
              (-1 / (alpha * lambda + 1) ^ (1 / alpha) - 
              lambda * (alpha * lambda + 1) ^ (-1 / alpha - 1) + 1) + 
              (alpha * (-y - 1 / alpha)) / (alpha * lambda + 1))
              
              G22 <- (1 - z) * ((alpha * lambda + 1) ^ (2 / alpha) * 
              (alpha ^ 3 * lambda ^ 4 + (2 * alpha ^ 2 - 2 * alpha ^ 3 * y) * lambda ^ 3 +
              (alpha - 5 * alpha ^ 2 * y) * lambda ^ 2 - 4 * alpha * y * lambda - y) + 
              (alpha * lambda + 1) ^ (1 / alpha) * ((alpha - alpha ^ 3) * lambda ^ 4 +
              ((4 * alpha ^ 3 + 4 * alpha ^ 2) * y - 4 * alpha ^ 2 - alpha + 1) * lambda ^ 3 +
              ((10 * alpha ^ 2 + 6 * alpha) * y - 3 * alpha - 1) * lambda ^ 2 + 
              (8 * alpha + 2) * y * lambda + 2 * y) + 
              ((-2 * alpha ^ 3 - 4 * alpha ^ 2 - 2 * alpha) * y + 2 * alpha ^ 2 + 2 * alpha) * lambda ^ 3 +
              ((-5 * alpha ^ 2 - 6 * alpha - 1) * y + 2 * alpha + 1) * lambda ^ 2 + 
              (-4 * alpha - 2) * y * lambda - y) / (lambda ^ 2 * (alpha * lambda + 1) ^ 2 * 
              ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2)
              
              G22 <- G22 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 + 
                      G2 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
              
              
              predNumbers <- cumsum(predNumbers)
              res[(predNumbers[2]+1):predNumbers[3], (predNumbers[2]+1):predNumbers[3]] <- #PI
                t(as.data.frame(X[(nrow(X) * 2 / 3 + 1):nrow(X), (predNumbers[2]+1):predNumbers[3]]) *
                    G00 * weight) %*% X[(nrow(X) * 2 / 3 + 1):nrow(X), (predNumbers[2]+1):predNumbers[3]]
              
              res[(predNumbers[1]+1):predNumbers[2], (predNumbers[2]+1):predNumbers[3]] <- 0
              
              res[(predNumbers[2]+1):predNumbers[3], (predNumbers[1]+1):predNumbers[2]] <- 0
              
              res[1:predNumbers[1], (predNumbers[2]+1):predNumbers[3]] <- 0
              
              res[(predNumbers[2]+1):predNumbers[3], 1:predNumbers[1]] <- 0
              
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
    PI     <-     piLink(eta[, 3], inverse = TRUE)
    mu <- mu.eta(eta = eta)
    
    logLikFit <- (
      (y == 1) * log(PI) + (y>1) * log(1 - PI) + (y>1) *
      (lgamma(y + 1 / alpha) - lgamma(1 / alpha) - lgamma(y + 1) - 
      (y + 1 / alpha) * log(1 + lambda * alpha) + y * log(lambda * alpha) - 
      log(1 - (1 + lambda * alpha) ^ (-1 / alpha) - lambda * (1 + lambda * alpha) ^ (-1 - 1 / alpha)))
    )
    
    yUnq <- unique(y)
    
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
      if (yNow < 3) return(0)
      
      xx <- findL(yNow)
      a <- exp(xx$par[3])
      l <- exp(xx$par[2])
      
      lgamma(yNow+1/a)-lgamma(1/a)-lgamma(yNow+1)-(yNow+1/a)*log(1+l*a)+yNow*log(l*a)-log(1-(1+l*a)^(-1/a)-l*(1+l*a)^(-1-1/a))
    }
    
    suppressWarnings(logLikIdeal <- sapply(yUnq, FUN = fff))
    names(logLikIdeal) <- yUnq
    
    logLikIdeal <- sapply(y, FUN = function(x) {
      logLikIdeal[names(logLikIdeal) == x]
    })
    diff <- logLikIdeal - logLikFit
    
    if (any(logLikFit > 0)) {
      warning("Dispertion parameter values are on the boundary of parameter space. Deviance residuals will be asigned 0 on these observations.")
      diff[logLikFit > 0]   <- 0
    } else if (any(diff < 0)) {
      warning("Numerical deviance finder found worse saturated likelihood than fitted model. Expect NA's in deviance/deviance residuals.")
    }
    
    diff[diff < 0] <- 0
    
    sign(y - mu) * sqrt(2 * wt * diff)
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    N <- pw * (1 - lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1)) / 
    (1 - (1 + alpha * lambda) ^ (- 1 / alpha) - 
    lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    bigTheta0 <- pw * 0 # w.r to PI
    bigTheta1 <- pw  *  alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
      ((lambda * alpha + 1) ^ (1 / alpha) * (lambda ^ 2 * alpha ^ 2 + 2 * lambda * alpha + 1) * 
      log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (1 / alpha) * 
      (-lambda ^ 2 * alpha ^ 2 - lambda * alpha) - lambda ^ 2 * alpha ^ 2) / 
      (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha + 1) - lambda * alpha - lambda - 1) ^ 2)# w.r to alpha
    bigTheta2 <- -pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
      ((alpha * lambda + 1) ^ (1 / alpha + 1) - 1) /
      ((alpha * lambda + 1) ^ (1 / alpha + 1) + (-alpha - 1) * lambda - 1) ^ 2# w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1, bigTheta0) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * (1 - lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1)) * 
    (1 + alpha * lambda) ^ (- 1 / alpha) / (1 - (1 + alpha * lambda) ^ (- 1 / alpha) - 
    lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1)) ^ 2)
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    PI     <-  piLink(eta[, 3], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        as.numeric(x == 1) * PI + as.numeric(x > 0) * 
        (1 - PI) * stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / 
        (1 - (1 + alpha * lambda) ^ (-1 / alpha) - 
        lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1))
      },
      "nontrunc" = {
        stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / 
        (1 - lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1)) *
        (as.numeric(x == 0) + as.numeric(x > 1) * (1 - PI)) + 
        as.numeric(x == 1) * PI * (1 - (1 + alpha * lambda) ^ (-1 / alpha) / 
        (1 - lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1)))
      }
    )
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <- alphaLink(eta[, 2], inverse = TRUE)
    PI     <- piLink(eta[, 3], inverse = TRUE)
    P0 <- (1 + alpha * lambda) ^ (- 1 / alpha)
    P1 <- lambda * (1 + alpha * lambda) ^ (- 1 / alpha - 1)
    CDF <- function(x) {
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, P0 / (1 - P1), 
      P0 / (1 - P1) + PI * (1 - P0 / (1 - P1)) + 
      (1 - PI) * (stats::pnbinom(x, mu = lambda, size = 1 / alpha) - P1 - P0) / (1 - P1))))
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
  
  getStart <- expression(
    if (method == "IRLS") {
      init <- log(abs((observed / weighted.mean(observed, priorWeights) - 1) / observed) + .1)
      etaStart <- cbind(
        pmin(family$links[[1]](observed), family$links[[1]](12)),
        family$links[[2]](ifelse(init < -.5, .1, init + .55)),
        family$links[[3]](weighted.mean(observed == 1, priorWeights) * (.5 + .5 * (observed == 1)) + .01)
      ) + offset
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
      if ("(Intercept):pi" %in% colnames(Xvlm)) {
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
      family    = "ztHurdlenegbin",
      etaNames  = c("lambda", "alpha", "pi"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

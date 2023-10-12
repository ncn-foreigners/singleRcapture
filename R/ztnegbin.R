#' @rdname singleRmodels
#' @importFrom stats uniroot
#' @importFrom stats dnbinom
#' @importFrom stats optim
#' @export
ztnegbin <- function(nSim = 1000, epsSim = 1e-8, eimStep = 6,
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
      switch (type,
          "nontrunc" = lambda,
          "trunc" = lambda / (1 - (1 + alpha * lambda) ^ (-1 / alpha))
      )
    } else {
      switch (type,
          "nontrunc" = {
            matrix(c(1, 0) * c(
              lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
              alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
            ), ncol = 2)
          },
          "trunc" = {
            matrix(c(
              (alpha * lambda + 1) ^ (1 / alpha - 1) *
              ((alpha * lambda + 1) ^ (1 / alpha + 1) +
              (-alpha - 1) * lambda - 1) / 
              ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2, 
              lambda * (lambda * alpha + 1) ^ (1 / alpha - 1) *
              ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha) /
              (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)
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
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    switch (type,
    nontrunc = lambda * (1 + alpha * lambda),
    trunc = (lambda + alpha * (lambda ^ 2) - alpha * (lambda ^ 2) * P0) / ((1 - P0) ^ 2)
    )
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
  
  compExpect <- function(eta) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    
    P0 <- (1 + alpha * lambda) ^ (-1 / alpha)
    #P0 <- stats::dnbinom(x = 0, size = 1 / alpha, mu = lambda)
    res <- rep(0, NROW(eta))
    k <- 1
    finished <- rep(FALSE, NROW(eta))
    while ((k < nSim) & !all(finished)) {
      prob <- apply(cbind(k:(k + eimStep)), MARGIN = 1, FUN = function(x) {
        stats::dnbinom(
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
    Ey <- mu.eta(eta = eta)
    Etrig <- compExpect(eta)
    
    G00 <- prior * (Etrig + ((lambda * alpha + 1) ^ (1 / alpha) * 
    (lambda ^ 2 * alpha ^ 2 + 2 * lambda * alpha + 1) * log(lambda * alpha + 1) ^ 2 +
    ((lambda * alpha + 1) ^ (1 / alpha) * (2 * lambda ^ 2 * alpha ^ 3 + 
    (4 * lambda - 2 * lambda ^ 2) * alpha ^ 2 + (2 - 2 * lambda) * alpha) + 
    (lambda * alpha + 1) ^ (2 / alpha) * (-2 * lambda ^ 2 * alpha ^ 3 - 4 * lambda * alpha ^ 2 - 2 * alpha)) * 
    log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (2 / alpha) * 
    ((3 * lambda ^ 2 - 2 * Ey * lambda) * alpha ^ 3 +
    (2 * lambda - Ey) * alpha ^ 2) + (lambda * alpha + 1) ^ (1 / alpha) *
    ((4 * Ey * lambda - 3 * lambda ^ 2) * alpha ^ 3 +
    (lambda ^ 2 - 2 * lambda + 2 * Ey) * alpha ^ 2) -
    2 * Ey * lambda * alpha ^ 3 - Ey * alpha ^ 2) /
    (alpha ^ 4 * (lambda * alpha + 1) ^ 2 * 
    ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)) * 
    alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2
    
    # mixed derivative
    G01 <- -((alpha * lambda + 1) ^ (1 / alpha + 1) *
    log(alpha * lambda + 1) + (alpha * lambda + 1) ^ (1 / alpha) *
    ((alpha ^ 2 - alpha) * lambda - 2 * alpha ^ 2 * Ey) +
    (alpha * lambda + 1) ^ (2 / alpha) * (alpha ^ 2 * Ey - alpha ^ 2 * lambda) + alpha ^ 2 * Ey) /
    (alpha ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)
    
    G01 <- G01 * prior * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
           lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    # second beta derivative
    
    G11 <- prior * ((alpha * lambda + 1) ^ (2 / alpha) * (alpha * lambda ^ 2 - 2 * alpha * Ey * lambda - Ey) +
    (alpha * lambda + 1) ^ (1 / alpha) * ((1 - alpha) * lambda ^ 2 + 4 * alpha * Ey * lambda + 2 * Ey) - 2 * alpha * Ey * lambda - Ey) /
    (lambda ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
    matrix(
      -c(G11, # lambda
         G01, # mixed
         G01, # mixed
         G00  # alpha
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "alpha")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
    dig <- compdigamma(y = y, alpha = alpha)
    weight <- weight / prior

    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    # log(alpha) derivative
    dig <- compdigamma(y = y, alpha = alpha)
    
    G0 <- (dig + ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) +
    (y - lambda) * alpha * (lambda * alpha + 1) ^ (1 / alpha) - y * alpha) /
    (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) - 1))) *
    alphaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    # Beta derivative
    G1 <- -((lambda - y) * (alpha * lambda + 1) ^ (1 / alpha) + y) /
    (lambda * (alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) - 1)) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
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
                           offset, ...) {
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
        
        -sum(weight * (lgamma(y + 1 / alpha) - lgamma(1 / alpha) - lgamma(y + 1) - 
        (y + 1 / alpha) * log(1 + lambda * alpha) + y * log(lambda * alpha) - 
        log(1 - (1 + lambda * alpha) ^ (-1 / alpha))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        alpha  <-  alphaLink(eta[, 2], inverse = TRUE)
        
        # log(alpha) derivative
        dig <- compdigamma(y = y, alpha = alpha)
        
        G0 <- (dig + ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) +
        (y - lambda) * alpha * (lambda * alpha + 1) ^ (1 / alpha) - y * alpha) /
        (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) - 1))) *
        alphaLink(eta[, 2], inverse = TRUE, deriv = 1) * weight
        
        # Beta derivative
        G1 <- -((lambda - y) * (alpha * lambda + 1) ^ (1 / alpha) + y) /
        (lambda * (alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) - 1)) *
        lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * weight
        
        if (NbyK) {
          XX <- 1:(attr(X, "hwm")[1])
          return(cbind(as.data.frame(X[1:nrow(eta), XX]) * G1, as.data.frame(X[-(1:nrow(eta)), -XX]) * G0))
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
        res <- matrix(nrow = length(beta), 
                      ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        trig <- comptrigamma(y = y, alpha = alpha)
        dig <-  compdigamma(y = y, alpha = alpha)
        
        # 2nd log(alpha) derivative
        G00 <- ((trig + ((lambda * alpha + 1) ^ (1 / alpha) * 
        (lambda ^ 2 * alpha ^ 2 + 2 * lambda * alpha + 1) * log(lambda * alpha + 1) ^ 2 +
        ((lambda * alpha + 1) ^ (1 / alpha) * (2 * lambda ^ 2 * alpha ^ 3 + 
        (4 * lambda - 2 * lambda ^ 2) * alpha ^ 2 + (2 - 2 * lambda) * alpha) + 
        (lambda * alpha + 1) ^ (2 / alpha) * (-2 * lambda ^ 2 * alpha ^ 3 - 4 * lambda * alpha ^ 2 - 2 * alpha)) * 
        log(lambda * alpha + 1) + (lambda * alpha + 1) ^ (2 / alpha) * 
        ((3 * lambda ^ 2 - 2 * y * lambda) * alpha ^ 3 +
        (2 * lambda - y) * alpha ^ 2) + (lambda * alpha + 1) ^ (1 / alpha) *
        ((4 * y * lambda - 3 * lambda ^ 2) * alpha ^ 3 +
        (lambda ^ 2 - 2 * lambda + 2 * y) * alpha ^ 2) -
        2 * y * lambda * alpha ^ 3 - y * alpha ^ 2) /
        (alpha ^ 4 * (lambda * alpha + 1) ^ 2 * 
        ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)) * 
        alphaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
        (dig + ((lambda * alpha + 1) ^ (1 / alpha + 1) * log(lambda * alpha + 1) +
        (y - lambda) * alpha * (lambda * alpha + 1) ^ (1 / alpha) - y * alpha) /
        (alpha ^ 2 * (lambda * alpha + 1) * ((lambda * alpha + 1) ^ (1 / alpha) - 1))) *
        alphaLink(eta[, 2], inverse = TRUE, deriv = 2))
        
        # mixed derivative
        G01 <- -((alpha * lambda + 1) ^ (1 / alpha + 1) *
        log(alpha * lambda + 1) + (alpha * lambda + 1) ^ (1 / alpha) *
        ((alpha ^ 2 - alpha) * lambda - 2 * alpha ^ 2 * y) +
        (alpha * lambda + 1) ^ (2 / alpha) * (alpha ^ 2 * y - alpha ^ 2 * lambda) + alpha ^ 2 * y) /
        (alpha ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)
        
        G01 <- G01 * alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
               lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
        
        # second beta derivative
        
        G11 <- (((alpha * lambda + 1) ^ (2 / alpha) * (alpha * lambda ^ 2 - 2 * alpha * y * lambda - y) +
        (alpha * lambda + 1) ^ (1 / alpha) * ((1 - alpha) * lambda ^ 2 + 4 * alpha * y * lambda + 2 * y) - 2 * alpha * y * lambda - y) /
        (lambda ^ 2 * (alpha * lambda + 1) ^ 2 * ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2) *
        lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 -
        ((lambda - y) * (alpha * lambda + 1) ^ (1 / alpha) + y) /
        (lambda * (alpha * lambda + 1) * ((alpha * lambda + 1) ^ (1 / alpha) - 1)) *
        lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)) * weight
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- 
          t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]) *
            G00 * weight) %*% X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
        
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- 
          t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * 
            G11 * weight) %*% X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- 
          t(t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * 
            G01 * weight) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]))
        
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- 
          t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * 
            G01 * weight) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
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
    mu <- mu.eta(eta = eta)
    
    logLikFit <- (
      lgamma(y + 1 / alpha) - lgamma(1 / alpha) - lgamma(y + 1) - 
      (y + 1 / alpha) * log(1 + alpha * lambda) +
      y * log(lambda * alpha) - log(1 - (1 + alpha * lambda) ^ (-1 / alpha))
    )
    
    yUnq <- unique(y)
    
    if (any(yUnq > 77)) {
      warning("Curently numerical deviance is unreliable for counts greater than 78.")
    }
    
    # findL <- function(t) {
    #   yNow <- yUnq[t]
    #   nleqslv::nleqslv(
    #     x = -c(.5, log(yNow), ifelse(yNow < 10, 1, 2)),
    #     fn = function(x) {
    #       s <- x[1]
    #       l <- exp(x[2])
    #       a <- exp(x[3])
    #        
    #       prob <- 1 - (1+a*l)^(-1/a)
    #       prob <- 1 / prob
    #       c(l*prob - yNow,# s der
    #         yNow/l+(1+yNow*a)/(1+l*a)+s*yNow*l-(yNow-l)/(l*(1+l*a))-s*yNow*(yNow-l)/(l*(1+l*a)),# lambda der
    #         (1+yNow*a)*l+(yNow-l)*(1+s*yNow)*((1+a*l)*log(1+a*l)-a*l)+l*(1+a*l)*(digamma(1/a)+log(a)+log(yNow+1/a))-(digamma(yNow+1/a)+1)*(1+a*l)*l)#alpha der
    #     },
    #     control = list(
    #       xtol = .Machine$double.eps
    #     )
    #   )$x
    # }
    
    ## This could be more stable but this will do for now
    
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
    
    ### for testing 
    # findL <- function(yNow) {
    #   stats::optim(
    #     par = c(0, .6, -14),
    #     fn = function(x) {
    #       s <- x[1]
    #       l <- exp(x[2])
    #       a <- exp(x[3])
    # 
    #       prob <- 1 - (1+a*l)^(-1/a)
    #       prob <- 1 / prob
    #       sum(c((l*prob - yNow) * 10,# s der
    #             yNow/l+(-yNow*a-1)/(1+a*l)-(1+a*l)^(-1-1/a)*prob+s*(prob-prob^2*(l*(1+a*l)^(-1-1/a))),# lambda der
    #             (log(l*a+1)/a^2-l/(a*(l*a+1)))/((l*a+1)^(1/a)*(1-1/(l*a+1)^(1/a)))+(s*l*(log(l*a+1)/a^2-l/(a*(l*a+1))))/((l*a+1)^(1/a)*(1-1/(l*a+1)^(1/a))^2)+log(l*a+1)/a^2+(l*(-1/a-yNow))/(l*a+1)+yNow/a-digamma(yNow+1/a)/a^2+digamma(1/a)/a^2) ^ 2) ^ .5#alpha der
    #     },
    #     method = "BFGS",
    #     control = list(maxit = 100000, abstol = .Machine$double.eps, reltol = .Machine$double.eps)
    #   )
    # }
    
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

    bigTheta1 <- pw  *  alphaLink(eta[, 2], inverse = TRUE, deriv = 1) *
    ((lambda * alpha + 1) ^ (1 / alpha - 1) * ((lambda * alpha + 1) * log(lambda * alpha + 1) - lambda * alpha)) /
    (alpha ^ 2 * ((lambda * alpha + 1) ^ (1 / alpha) - 1) ^ 2)# w.r to alpha
    bigTheta2 <- -pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
    ((alpha * lambda + 1) ^ (1 / alpha - 1) / ((alpha * lambda + 1) ^ (1 / alpha) - 1) ^ 2)# w.r to lambda

    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)

    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * (1 - pr) / (pr ^ 2))

    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    alpha  <-  alphaLink(eta[, 2], inverse = TRUE)

    switch (type,
      "trunc" = {
        stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / 
        (1 - (1 + alpha * lambda) ^ (-1 / alpha))
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
      family    = "ztnegbin",
      etaNames  = c("lambda", "alpha"),
      simulate  = simulate,
      getStart  = getStart,
      extraInfo = c(
        mean       = "lambda",
        variance   = "lambda * (1 + alpha * lambda)",
        popSizeEst = "1 / (1 - (1 + alpha * lambda) ^ (- 1 / alpha))",
        meanTr     = "lambda / (1 - (1 + alpha * lambda) ^ (-1 / alpha))",
        varianceTr =
          "(lambda + alpha * (lambda ^ 2) - alpha * (lambda ^ 2) * (1 + alpha * lambda) ^ (-1 / alpha)) / ((1 - (1 + alpha * lambda) ^ (-1 / alpha)) ^ 2)"
      )
    ),
    class = c("singleRfamily", "family")
  )
}

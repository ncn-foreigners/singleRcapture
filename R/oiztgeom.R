#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @export
oiztgeom <- function(lambdaLink = c("log", "neglog"), 
                     omegaLink = c("logit", "cloglog", "probit"), 
                     ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  if (missing(omegaLink))  omegaLink <- "logit"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink, omegaLink)
  
  lambdaLink <- switch (lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  omegaLink <- switch (omegaLink,
    "logit"   = singleRinternallogitLink,
    "cloglog" = singleRinternalcloglogLink,
    "probit" = singleRinternalprobitLink
  )
  
  links[1:2] <- c(lambdaLink, omegaLink)
  
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    if (!deriv) {
      switch (type,
        "nontrunc" = omega + lambda * (1 - omega),
        "trunc" = (lambda + omega - lambda * omega) * (1 + lambda) / (lambda + omega)
      )
    } else {
      switch (type,
        "nontrunc" = {
          matrix(c(1 - omega, 1 - lambda) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
             omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        },
        "trunc" = {
          matrix(c(
            (1 - omega) * (lambda + 2 * omega) * lambda / (lambda + omega) ^ 2, 
            -(1 + lambda) * lambda ^ 2 / (lambda + omega) ^ 2
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
             omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        }
      )
    }
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    p1 <- (lambda ^ 2) * omega + lambda * omega + omega + lambda
    p1 <- p1 / (lambda ^ 2 + lambda * omega + omega + lambda)
    switch (type,
    "nontrunc" = omega + (1 - omega) * lambda * (2 * lambda + 1),
    "trunc" = p1 + (1 - omega ) * (lambda ^ 2) * (2 * (lambda ^ 2) + 5 * lambda + 4) / ((1 + lambda) * (lambda + omega))
    ) - (mu.eta(eta = eta, type = type) ^ 2)
  }
  
  Wfun <- function(prior, eta, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- (lambda ^ 2) * omega + lambda * omega + lambda + omega
    z <- z / (lambda ^ 2 + lambda * omega + lambda + omega) # expected for I's
    Ey <- mu.eta(eta)
    XXX <- Ey - z ## z here is the prob of 1 and XXX is expected for (1-z)*y
    
    # omega^2 derivative
    G00 <- (1 / ((lambda + omega) ^ 2) - z * ((lambda ^ 2 + lambda + 1) ^ 2) / 
    ((omega * (lambda ^ 2 + lambda + 1) + lambda) ^ 2) - 
    (1 - z) / ((1 - omega) ^ 2))
    
    G00 <- prior * G00 * omegaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2
    
    # mixed derivative
    G01 <- ((1 - z) / ((omega + lambda) ^ 2) + 
    z * lambda * (omega * omega * (lambda ^ 3 + 2 * (lambda ^ 2) + 4 * lambda + 2) + 
                  omega * (4 * (lambda ^ 2) + 2 * lambda) + lambda ^ 3) / 
    (((omega + lambda) * (omega * (lambda ^ 2 + lambda + 1) + lambda)) ^ 2))
    
    G01 <- G01 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
    omegaLink(eta[, 2], inverse = TRUE, deriv = 1) * prior
    
    # Beta^2 derivative
    G11 <- (z * (((omega + 2 * lambda + 1) ^ 2) / 
    ((omega * lambda + omega + lambda ^ 2 + lambda) ^ 2) - 
    2 / (omega * lambda + omega + lambda ^ 2 + lambda) + 
    (2 * omega) / (omega * (lambda ^ 2) + omega * lambda + omega + lambda) - 
    ((2 * omega * lambda + omega + 1) ^ 2) / 
    ((omega * (lambda ^ 2) + omega * lambda + omega + lambda) ^ 2)) + 
    (XXX / ((1 + lambda) ^ 2) - 
    XXX / (lambda ^ 2) + (1 - z) / ((omega + lambda) ^ 2)))
    
    G11 <- prior * G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
    matrix(
      -c(G11, # lambda
         G01, # mixed
         G01, # mixed
         G00  # omega
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "omega")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- ifelse(y == 1, y, 0)
    weight <- weight / prior
    
    G1 <- z * (omega - 1) * lambda * (omega * (2 + lambda) + lambda) / ((1 + lambda) * (lambda + omega) * (omega * (lambda ^ 2 + lambda + 1) + lambda))
    G1 <- G1 + (1 - z) * (y / lambda - y / (1 + lambda) - 1 / (omega + lambda))
    G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    
    G0 <- z * (lambda + 1) * (lambda ^ 2) / ((omega + lambda) * (omega * (lambda ^ 2 + lambda + 1) + lambda))
    G0 <- G0 - (1 - z) * (lambda + 1) / ((1 - omega) * (omega + lambda)) # omega derivative
    G0 <- G0 * omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      #xx <- solve(weight[[x]])
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
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X) / 2), rep(0, NROW(X) / 2))
    }
    z <- as.numeric(y == 1)
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        -sum(weight * (z * (log((lambda ^ 2) * omega + lambda * omega + lambda + omega) - 
        log(lambda ^ 2 + lambda * omega + lambda + omega)) + 
        (1 - z) * (log(1 - omega) + y * log(lambda) - y * log(1 + lambda) - log(lambda + omega))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- (z * (omega - 1) * lambda * (omega * (2 + lambda) + lambda) / 
        ((1 + lambda) * (lambda + omega) * (omega * (lambda ^ 2 + lambda + 1) + lambda)) + 
        (1 - z) * (y / lambda - y / (1 + lambda) - 1 / (omega + lambda)))
        
        G0 <- (z * (lambda + 1) * (lambda ^ 2) / 
        ((omega + lambda) * (omega * (lambda ^ 2 + lambda + 1) + lambda)) - 
        (1 - z) * (lambda + 1) / ((1 - omega) * (omega + lambda)))
        
        G1 <- G1 * weight * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
        G0 <- G0 * weight *  omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
        
        if (NbyK) {
          XX <- 1:(attr(X, "hwm")[1])
          return(cbind(as.data.frame(X[1:nrow(eta), XX]) * G1, as.data.frame(X[-(1:nrow(eta)), -XX]) * G0))
        } 
        
        if (vectorDer) {
          return(cbind(G1, G0))
        }
        as.numeric(c(G1, G0) %*% X)
      },
      function (beta) {
          lambdaPredNumber <- attr(X, "hwm")[1]
          eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
          omega  <-  omegaLink(eta[, 2], inverse = TRUE)
          lambda <- lambdaLink(eta[, 1], inverse = TRUE)
          
          res <- matrix(nrow = length(beta), ncol = length(beta), 
                        dimnames = list(names(beta), names(beta)))
          
          G1 <- (z * (omega - 1) * lambda * (omega * (2 + lambda) + lambda) / 
          ((1 + lambda) * (lambda + omega) * (omega * (lambda ^ 2 + lambda + 1) + lambda)) + 
          (1 - z) * (y / lambda - y / (1 + lambda) - 1 / (omega + lambda)))
          
          G0 <- (z * (lambda + 1) * (lambda ^ 2) / 
          ((omega + lambda) * (omega * (lambda ^ 2 + lambda + 1) + lambda)) - 
          (1 - z) * (lambda + 1) / ((1 - omega) * (omega + lambda)))
          
          # omega^2 derivative
          G00 <- (1 / ((lambda + omega) ^ 2) - z * ((lambda ^ 2 + lambda + 1) ^ 2) / 
          ((omega * (lambda ^ 2 + lambda + 1) + lambda) ^ 2) - 
          (1 - z) / ((1 - omega) ^ 2))
          
          G00 <- weight * (G0 *  omegaLink(eta[, 2], inverse = TRUE, deriv = 2) + 
                           G00 * omegaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2)
          
          G00 <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)] * G00)) %*% 
          as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
          
          # mixed derivative
          G01 <- ((1 - z) / ((omega + lambda) ^ 2) + 
          z * lambda * (omega * omega * (lambda ^ 3 + 2 * (lambda ^ 2) + 4 * lambda + 2) + 
                        omega * (4 * (lambda ^ 2) + 2 * lambda) + lambda ^ 3) / 
          (((omega + lambda) * (omega * (lambda ^ 2 + lambda + 1) + lambda)) ^ 2))
          
          G01 <- G01 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
          omegaLink(eta[, 2], inverse = TRUE, deriv = 1) * weight
          
          G01 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * G01) %*% 
          as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
          
          # Beta^2 derivative
          G11 <- (z * (((omega + 2 * lambda + 1) ^ 2) / 
          ((omega * lambda + omega + lambda ^ 2 + lambda) ^ 2) - 
          2 / (omega * lambda + omega + lambda ^ 2 + lambda) + 
          (2 * omega) / (omega * (lambda ^ 2) + omega * lambda + omega + lambda) - 
          ((2 * omega * lambda + omega + 1) ^ 2) / 
          ((omega * (lambda ^ 2) + omega * lambda + omega + lambda) ^ 2)) + 
          (1 - z) * (y / ((1 + lambda) ^ 2) - 
          y / (lambda ^ 2) + 1 / ((omega + lambda) ^ 2)))
          
          G11 <- (G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 + 
                  G1 *  lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)) * weight # second derivative of log link
          
          G11 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * G11)) %*% 
          X[1:(nrow(X) / 2), 1:lambdaPredNumber]
          
          res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
          res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
          res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(G01)
          res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- G01
          
          res
        }
    )
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  devResids <- function(y, eta, wt, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    mu1 <- mu.eta(eta = eta)
    idealOmega <- ifelse(y == 1, 1, 0)
    idealLambda <- ifelse(y > 1, y - 1, 0)
    diff <- ifelse(
      y == 1,
      -(log((lambda ^ 2) * omega + lambda * omega + lambda + omega) - log(lambda ^ 2 + lambda * omega + lambda + omega)),
      (y - 1) * log(idealLambda) - y * log(1 + idealLambda) - log(1 - omega) + log(omega + lambda) - y * log(lambda / (1 + lambda))
    )
    sign(y - mu1) * sqrt(2 * wt * diff)
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    N <- pw * (1 + lambda) / (lambda + omega)
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    bigTheta1 <- -pw * omegaLink(eta[, 2], inverse = TRUE, deriv = 1) * (1 + lambda) / ((lambda + omega) ^ 2) # w.r to omega
    bigTheta2 <- pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * (omega - 1) / ((lambda + omega) ^ 2) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <- sum(pw * (1 - omega) / (lambda + omega))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        ifelse(x == 1, ((lambda ^ 2) * omega + lambda * omega + lambda + omega) / 
        (lambda ^ 2 + lambda * omega + lambda + omega), 
        ((lambda / (1 + lambda)) ^ x) * (1 - omega) / (lambda + omega))
      },
      "nontrunc" = {
        (1 - omega) * stats::dnbinom(x = x, size = 1, mu = lambda) +
        omega * as.numeric(x == 1)
      }
    )
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    CDF <- function(x) {
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, 
      (1 - omega) * stats::pnbinom(x, mu = lambda, size = 1), 
      omega +  (1 - omega) * stats::pnbinom(x, mu = lambda, size = 1))))
    }
    lb <- CDF(rep(lower, n))
    ub <- CDF(rep(upper, n))
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
    if (method == "IRLS") {
      etaStart <- cbind(
        pmin(family$links[[1]](observed), family$links[[1]](12)),
        (sizeObserved * priorWeights * (observed == 1) + .5) / (sizeObserved * priorWeights * sum(observed == 1) + 1)
      ) + offset
    } else if (method == "optim") {
      init <- c(
        family$links[[1]](weighted.mean(observed, priorWeights)),
        family$links[[2]](weighted.mean(observed == 1, priorWeights) + .01)
      )
      if (attr(terms, "intercept")) {
        coefStart <- c(init[1], rep(0, attr(Xvlm, "hwm")[1] - 1))
      } else {
        coefStart <- rep(init[1] / attr(Xvlm, "hwm")[1], attr(Xvlm, "hwm")[1])
      }
      if ("(Intercept):omega" %in% colnames(Xvlm)) {
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
      family    = "oiztgeom",
      etaNames  = c("lambda", "omega"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

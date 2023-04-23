#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @export
ztoipoisson <- function(lambdaLink = c("log", "neglog"), 
                        omegaLink = c("logit", "cloglog", "probit"), 
                        ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  if (missing(omegaLink))  omegaLink <- "logit"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink, omegaLink)

  lambdaLink <- switch(lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  omegaLink <- switch(omegaLink,
    "logit" = singleRinternallogitLink,
    "cloglog" = singleRinternalcloglogLink,
    "probit" = singleRinternalprobitLink
  )
  
  links[1:2] <- c(lambdaLink, omegaLink)
  
  mu.eta <- function(eta, type = "trunc", ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    switch (type,
    "nontrunc" = omega * (1 - exp(-lambda)) + lambda * (1 - omega),
    "trunc" = omega + (1 - omega) * lambda / (1 - exp(-lambda))
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    switch (type,
    "nontrunc" = omega * (1 - exp(-lambda)) + (1 - omega) * (lambda + lambda ^ 2),
    "trunc" = omega + (1 - omega) * (lambda ^ 2 + lambda) / (1 - exp(-lambda))
    ) - mu.eta(type = type, eta = eta) ^ 2
  }
  
  Wfun <- function(prior, y, eta, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- omega + (1 - omega) * lambda / (exp(lambda) - 1) # expected for I's
    #z <- ifelse(y == 1, y, 0)
    
    G00 <- prior * ((-(z * (1 - lambda / (exp(lambda) - 1)) ^ 2) / 
    (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) ^ 2 - 
    (1 - z) / (1 - omega) ^ 2) * (omegaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2) + 
    ((z * (1 - lambda / (exp(lambda) - 1))) / 
    (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) - 
    (1 - z) / (1 - omega)) * omegaLink(eta[, 2], inverse = TRUE, deriv = 2))
    
    # mixed derivative
    G01 <- prior * (((z * ((lambda - 1) * exp(lambda) + 1)) / 
    (omega * exp(lambda) + (1 - omega) * lambda - omega) ^ 2) * 
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
    omegaLink(eta[, 2], inverse = TRUE, deriv = 1))
    
    #expected value of (1-z)*y
    XXXX <- (1 - omega) * lambda * (exp(lambda) - 1) / (exp(lambda) - 1)
    
    G11 <- prior * (lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)  * 
    ((1 - z) * (-exp(lambda) / (exp(lambda) - 1)) + XXXX / lambda +
    (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) /
    (exp(lambda) - 1) ^ 2)) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega)) +
    (lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2) * 
    ((1 - z) * (exp(2 * lambda) / (exp(lambda) - 1) ^ 2 - exp(lambda) / (exp(lambda) - 1))-
    XXXX / lambda ^ 2 + (z * ((2 * (1 - omega) * lambda * exp(2 * lambda)) / (exp(lambda) - 1) ^ 3-
    ((1 - omega) * lambda * exp(lambda)) / (exp(lambda) - 1) ^ 2-
    (2 * (1 - omega) * exp(lambda)) / (exp(lambda) - 1) ^ 2)) /
    (((1 - omega) * lambda) / (exp(lambda) - 1) + omega) - 
    (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda))/
    (exp(lambda) - 1) ^ 2) ^ 2) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega) ^ 2))
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
  
  funcZ <- function(eta, weight, y, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- ifelse(y == 1, y, 0)
    G0 <- (z * (1 - lambda / (exp(lambda) - 1))) / (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) - (1 - z) / (1 - omega)
    
    G1 <- ((1 - z) * (y / lambda - exp(lambda) / (exp(lambda) - 1)) + 
    (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) / 
    (exp(lambda) - 1) ^ 2)) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega))
    
    G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    G0 <- G0 *  omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- solve(weight[[x]]) # more stable
      xx %*% uMatrix[x, ]
    })
    pseudoResid <- t(pseudoResid)
    dimnames(pseudoResid) <- dimnames(eta)
    pseudoResid
  }
  
  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    z <- as.numeric(y == 1)
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        -sum(weight * (-log(1 + exp(eta[, 2])) + z * log(exp(eta[, 2]) + lambda / (exp(lambda) - 1)) +
        (1 - z) * (y * log(lambda) - log(exp(lambda) - 1) - log(factorial(y)))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G0 <- (z * (1 - lambda / (exp(lambda) - 1))) / (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) - (1 - z) / (1 - omega)
        
        G1 <- ((1 - z) * (y / lambda - exp(lambda) / (exp(lambda) - 1)) + 
        (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) / 
        (exp(lambda) - 1) ^ 2)) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega))
        
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
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)

        res <- matrix(nrow = length(beta), ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        # omega^2 derivative
        domega <- (z * (1 - lambda / (exp(lambda) - 1))) / 
        (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) - 
        (1 - z) / (1 - omega)
        G00 <- -(z * (1 - lambda / (exp(lambda) - 1)) ^ 2) / 
        (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) ^ 2 - 
        (1 - z) / (1 - omega) ^ 2
        
        G00 <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)] * 
        (G00 * (omegaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2) + 
        domega * omegaLink(eta[, 2], inverse = TRUE, deriv = 2)) * weight)) %*% 
        as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        # mixed derivative
        
        G01 <- ((z * ((lambda - 1) * exp(lambda) + 1)) / 
        (omega * exp(lambda) + (1 - omega) * lambda - omega) ^ 2)
        
        G01 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * 
        G01 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
        omegaLink(eta[, 2], inverse = TRUE, deriv = 1) * weight) %*% 
        as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        # Beta^2 derivative
        G11 <- ((1 - z) * (exp(2 * lambda) / (exp(lambda) - 1) ^ 2-
        exp(lambda) / (exp(lambda) - 1) - y / lambda ^ 2) + 
        (z * ((2 * (1 - omega) * lambda * exp(2 * lambda)) / (exp(lambda) - 1) ^ 3-
        ((1 - omega) * lambda * exp(lambda)) / (exp(lambda) - 1) ^ 2 - 
        (2 * (1 - omega) * exp(lambda)) / (exp(lambda) - 1) ^ 2)) / 
        (((1 - omega) * lambda) / (exp(lambda) - 1) + omega) - 
        (z * ((1 - omega) /(exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) /
        (exp(lambda) - 1) ^ 2) ^ 2) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega) ^ 2)
        
        dlambda <- ((1 - z) * (y / lambda - exp(lambda) / (exp(lambda) - 1)) + 
        (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) / 
        (exp(lambda) - 1) ^ 2)) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega))
        
        G11 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * 
        (G11 * (lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2) + 
        dlambda * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)) * weight)) %*% 
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
    mu <- mu.eta(eta = eta)
    #idealOmega <- ifelse(y == 1, 1, 0)
    idealLambda <- ifelse(y > 1, lamW::lambertW0(-y * exp(-y)) + y, 0)
    diff <- ifelse(
      y == 1,
      -(log(omega + (1 - omega) * lambda / (exp(lambda) - 1))),
      y * (log(idealLambda) - log(lambda)) + log((exp(lambda) - 1) / (exp(idealLambda) - 1)) - log(1 - omega)
    )
    sign(y - mu) * sqrt(2 * wt * diff)
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    N <- pw / (1 - exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    bigTheta1 <- rep(0, nrow(eta)) # w.r to omega
    bigTheta2 <- pw * (exp(lambda) / (exp(lambda) - 1) ^ 2) # w.r to lambda
    bigTheta2 <- bigTheta2 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)# w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <- t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * exp(-lambda) / (1 - exp(-lambda)) ^ 2)
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    ifelse(x == 1, omega + (1 - omega) * lambda / (exp(lambda) - 1),
    (1 - omega) * (lambda ^ x) / (factorial(x) * (exp(lambda) - 1)))
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
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
    start <- stats::glm.fit(
      x = variables[wch$reg, 1:attr(Xvlm, "hwm")[1]],
      y = observed[wch$reg],
      family = stats::poisson(),
      weights = priorWeights[wch$reg],
      ...
    )$coefficients,
    if (attr(family$links, "linkNames")[1] == "neglog") start <- -start,
    if (is.null(controlMethod$omegaStart)) {
      if (controlModel$omegaFormula == ~ 1) {
        omg <- (length(observed[wch$reg]) - sum(observed == 1)) / (sum(observed[wch$reg]) - length(observed[wch$reg]))
        start <- c(start, family$links[[2]](omg))
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
      family    = "ztoipoisson",
      etaNames  = c("lambda", "omega"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

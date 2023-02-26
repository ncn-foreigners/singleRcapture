#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @export
ztoipoisson <- function(...) {
  # Fist for lambda second for omega
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2]/ (1 - x[,2]))), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),1/(exp(-x[,2]) + 1)), ncol = 2, dimnames = dimnames(x))}
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
    "nontrunc" = omega * (1 - exp(-lambda)) + lambda * (1 - omega),
    "trunc" = omega + (1 - omega) * lambda * exp(lambda) / (exp(lambda) - 1)
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
    "nontrunc" = omega * (1 - exp(-lambda)) + (1 - omega) * (lambda + lambda ^ 2),
    "trunc" = omega + (1 - omega) * (lambda ^ 2 + lambda) / (1 - exp(-lambda))
    ) - mu.eta(type = type, eta = eta) ^ 2
  }
  
  Wfun <- function(prior, y, eta, ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- omega + (1 - omega) * lambda / (exp(lambda) - 1) # expected for I's
    #z <- ifelse(y == 1, y, 0)
    G00 <- -omega * (1 - omega) + z * (exp(eta[, 2]) / (exp(eta[, 2]) + lambda / (exp(lambda) - 1)) - (exp(eta[, 2]) / (exp(eta[, 2]) + lambda / (exp(lambda) - 1))) ^ 2) * prior
    
    # mixed derivative
    G01 <- -z * exp(eta[, 2]) / ((exp(eta[, 2]) + lambda / (exp(lambda) - 1)) ^ 2) * (lambda / (exp(lambda) - 1) - (lambda / (exp(lambda) - 1) ^ 2) * exp(lambda)) * prior
    
    wwd <- 1 / (exp(lambda) - 1)
    thetaLambdaTerm <- 1 / (exp(eta[, 2]) + lambda * wwd)
    G11 <- -(1 - z) * (1 / (1 - exp(-lambda)) - lambda * exp(-lambda) / ((1 - exp(-lambda)) ^ 2)) + 
    z * (-lambda * ((wwd - exp(lambda) * lambda * wwd * wwd) ^ 2) / ((exp(eta[, 2]) + lambda * wwd) ^ 2) +
    (wwd - wwd * wwd * lambda * exp(lambda)) / (exp(eta[, 2]) + lambda * wwd) + lambda *
    (-wwd * wwd * lambda * exp(lambda) + 2 * lambda * exp(2 * lambda) * wwd * wwd * wwd - 
    2 * wwd * wwd * exp(lambda)) / (exp(eta[, 2]) + lambda * wwd)) * lambda * prior
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
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- ifelse(y == 1, y, 0)
    mm <- (exp(lambda) - 1) ^ -1
    G1 <- z * lambda * (mm - mm * mm * lambda * exp(lambda)) / (exp(eta[, 2]) + lambda * mm) + lambda * (1 - z) * (y / lambda - mm * exp(lambda))
    G0 <- -omega + z * (exp(eta[, 2]) / (exp(eta[, 2]) + lambda / (exp(lambda) - 1)))
    
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
    deriv <- deriv + 1 # to make it comfort to how swith in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        lambda <- invlink(eta)
        omega <- lambda[, 2]
        lambda <- lambda[, 1]
        -sum(weight * (-log(1 + exp(eta[, 2])) + z * log(exp(eta[, 2]) + lambda / (exp(lambda) - 1)) +
        (1 - z) * (y * log(lambda) - log(exp(lambda) - 1) - log(factorial(y)))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        lambda <- invlink(eta)
        omega <- lambda[, 2]
        lambda <- lambda[, 1]
        theta <- exp(eta[, 2])
        mm <- (exp(lambda) - 1) ^ -1
        G1 <- weight * lambda * z * (mm - mm * mm * lambda * exp(lambda)) / (exp(eta[, 2]) + lambda * mm) + weight * lambda * (1 - z) * (y / lambda - mm * exp(lambda))
        G0 <- -omega * weight + z * weight * (exp(eta[, 2]) / (exp(eta[, 2]) + lambda / (exp(lambda) - 1))) # omega derivative
        
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
        lambda <- invlink(eta)
        omega <- lambda[, 2]
        lambda <- lambda[, 1]
        theta <- exp(eta[, 2])
        res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
        
        # omega^2 derivative
        term <- lambda / (exp(lambda) - 1)
        G00 <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)] * (-omega * (1 - omega) + z * (theta / (theta + term) - (theta / (theta + term)) ^ 2)) * weight)) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        # mixed derivative
        G01 <- -t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * (term - (term ^ 2) * exp(lambda)) * z * theta / ((theta + term) ^ 2) * weight) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        # Beta^2 derivative
        wwd <- 1 / (exp(lambda) - 1)
        G11 <- (-(1 - z) * (1 / (1 - exp(-lambda)) - lambda * exp(-lambda) / ((1 - exp(-lambda)) ^ 2)) + 
        z * (-lambda * ((wwd - exp(lambda) * lambda * wwd * wwd) ^ 2) / ((theta + lambda * wwd) ^ 2) +
        (wwd - wwd * wwd * lambda * exp(lambda)) / (theta + lambda * wwd) + lambda / (theta + lambda * wwd) *
        (-wwd * wwd * lambda * exp(lambda) + 2 * lambda * exp(2 * lambda) * wwd * wwd * wwd - 2 * wwd * wwd * exp(lambda)))) * lambda * weight
        G11 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * G11)) %*% X[1:(nrow(X) / 2), 1:lambdaPredNumber]
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
    omega <- invlink(eta)
    lambda <- omega[, 1]
    omega <- omega[, 2]
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
    lambda <- invlink(eta)
    lambda <- lambda[, 1]
    N <- pw / (1 - exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    theta <- exp(eta[, 2])
    ml <- (1 - exp(-lambda)) ^ 2
    
    bigTheta1 <- rep(0, nrow(eta)) # w.r to omega
    bigTheta2 <- pw * (exp(log(lambda) - lambda) / ml) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * exp(-lambda) / ml)
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
    ifelse(x == 1, omega + (1 - omega) * lambda / (exp(lambda) - 1),
           (1 - omega) * (lambda ^ x) / (factorial(x) * (exp(lambda) - 1)))
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- invlink(eta)
    omega <- lambda[, 2]
    lambda <- lambda[, 1]
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
      linkfun = link,
      linkinv = invlink,
      mu.eta = mu.eta,
      link = c("log", "logit"),
      valideta = function (eta) {TRUE},
      variance = variance,
      Wfun = Wfun,
      funcZ = funcZ,
      devResids = devResids,
      validmu = validmu,
      pointEst = pointEst,
      popVar = popVar,
      family = "ztoipoisson",
      parNum = 2,
      etaNames = c("lambda", "omega"),
      densityFunction = dFun,
      simulate = simulate,
      getStart = getStart
    ),
    class = "family"
  )
}

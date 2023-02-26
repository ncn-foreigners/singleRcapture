#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @export
Hurdleztpoisson <- function(...) {
  # Fist for lambda second for PI
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2]/ (1 - x[,2]))), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),1/(exp(-x[,2]) + 1)), ncol = 2, dimnames = dimnames(x))}
  
  mu.eta <- function(eta, type = "trunc", ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    switch (type,
    "nontrunc" = PI + (1 - PI) * exp(-lambda) * (lambda * exp(lambda) - lambda) / (1 - lambda * exp(-lambda)),
    "trunc" = PI * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) + (1 - PI) * lambda * (exp(lambda) - 1) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    )
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    switch(type,
    "nontrunc" = PI + (1 - PI) * exp(-lambda) * lambda * (exp(lambda) * (1 + lambda) - 1) / (1 - lambda * exp(-lambda)),
    "trunc" = PI * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) + (1 - PI) * exp(-lambda) * lambda * (exp(lambda) * (1 + lambda) - 1) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    ) - (mu.eta(eta = eta, type = type) ^ 2)
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- (1 - PI) * exp(-lambda) / (1 - lambda * exp(-lambda))
    
    G00 <- -1 + PI * (1 - PI) * exp(-2 * lambda) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) ^ 2) - (1 - 2 * PI) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    G00 <- G00 * PI * (1 - PI) * prior
    
    # mixed
    
    G01 <- -PI * (1 - PI) * (-exp(-lambda) * (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) - exp(-2 * lambda) * (lambda - PI)) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) ^ 2)
    G01 <- G01 * lambda * prior
    
    # Beta^2 derivative
    G11 <- -(PI * (exp(lambda) * (lambda + 1) - 1) + exp(lambda) * (-lambda ^ 2 - 2) + exp(2 * lambda) + 1) / ((PI + exp(lambda) - lambda - 1) ^ 2) - (exp(lambda) * (lambda ^ 2 - lambda - exp(lambda) + 1) * z) / ((lambda - exp(lambda)) ^ 2)
    G11 <- G11 * lambda * prior
    
    matrix(
      -c(G11, # lambda
         G01, # mixed
         G01, # mixed
         G00  # pi
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "pi")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    z <- ifelse(y == 1, y, 0)
    G1 <- z * lambda * (lambda - 1) / (exp(lambda) - lambda)
    G1 <- G1 + (1 - z) * (y - lambda)
    G1 <- G1 - lambda * exp(-lambda) * (lambda - PI) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    G0 <- (z - PI- PI * (1 - PI) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))) # PI derivative
    
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
        PI <- lambda[, 2]
        lambda <- lambda[, 1]
        logistic <- -sum(weight * (z * log(PI) + (1 - z) * log(1 - PI)))
        rest <- -sum(z * log(1 - lambda * exp(-lambda)) +
                       (1 - z) * (y * log(lambda) - lambda - log(factorial(y))) -
                       log(1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)))
        rest + logistic
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        lambda <- invlink(eta)
        PI <- lambda[, 2]
        lambda <- lambda[, 1]
        G1 <- z * lambda * (lambda - 1) / (exp(lambda) - lambda)
        G1 <- G1 + (1 - z) * (y - lambda)
        G1 <- G1 - lambda * exp(-lambda) * (lambda - PI) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
        G1 <- G1 * weight
        G0 <- (z - PI- PI * (1 - PI) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))) # PI derivative
        G0 <- G0 * weight
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
        PI <- lambda[, 2]
        lambda <- lambda[, 1]
        Xlambda <- X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        XPI <- X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
        res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
        
        # PI^2 derivative
        G00 <- -1 + PI * (1 - PI) * exp(-2 * lambda) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) ^ 2) - (1 - 2 * PI) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
        G00 <- G00 * PI * (1 - PI)
        G00 <- t(as.data.frame(XPI * G00 * weight)) %*% as.matrix(XPI)
        
        # mixed
        
        G01 <- -PI * (1 - PI) * (-exp(-lambda) * (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) - exp(-2 * lambda) * (lambda - PI)) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) ^ 2)
        G01 <- G01 * lambda
        G01 <- t(as.data.frame(Xlambda * G01 * weight)) %*% as.matrix(XPI)
        
        # Beta^2 derivative
        G11 <- -(PI * (exp(lambda) * (lambda + 1) - 1) + exp(lambda) * (-lambda ^ 2 - 2) + exp(2 * lambda) + 1) / ((PI + exp(lambda) - lambda - 1) ^ 2) - (exp(lambda) * (lambda ^ 2 - lambda - exp(lambda) + 1) * z) / ((lambda - exp(lambda)) ^ 2)
        G11 <- G11 * lambda
        G11 <- t(as.data.frame(Xlambda * G11 * weight)) %*% Xlambda
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- G01
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- t(G01)
        
        res
      }
    )
  }
  
  validmu <- function(mu) {
    all(0 < mu & is.finite(mu))
  }
  
  devResids <- function(y, eta, wt, ...) {
    PI <- invlink(eta)
    lambda <- PI[, 1]
    PI <- PI[, 2]
    
    # when pi = 0 distribution collapses to zotpoisson
    inverseFunction <- function(y) {stats::uniroot(
      f = function(x) {
        lambda <- exp(x)
        (lambda - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda)) - y
      }, 
      lower = -log(y), upper = y * 10, 
      tol = .Machine$double.eps
    )$root}
    
    yUnq <- unique(y)
    etaSat <- sapply(yUnq, FUN = function(x) ifelse(x %in% c(1, 2), -Inf, inverseFunction(x)))
    etaSat <- sapply(y, FUN = function(x) etaSat[yUnq == x])
    idealLambda <- exp(etaSat)
    diff <- ifelse(
      y == 1,
      -(log(PI) + log(1 - lambda * exp(-lambda)) - log(1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))),
      ifelse(y == 2, log(2),
      y * etaSat - idealLambda - log(1 - exp(-idealLambda) - idealLambda * exp(-idealLambda))) - (log(1 - PI) + log(1 - lambda * exp(-lambda)) -log(1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)))
    )
    sign(y - mu.eta(eta = eta)) * sqrt(2 * wt * diff)
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    N <- pw * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    prob <- 1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)
    
    bigTheta1 <- -pw * (1 - lambda * exp(-lambda)) * exp(-lambda) / (prob ^ 2) # w.r to PI
    bigTheta2 <- pw * (prob * (lambda * exp(-lambda) - exp(-lambda)) - ((1 - PI) * exp(-lambda) + lambda * exp(-lambda) - exp(-lambda))) / (prob ^ 2) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    #f2 <- sum(pw * (1 - PI) * exp(-lambda) * (1 - lambda * exp(-lambda)) / (prob ^ 2))
    f2 <- sum(pw * (1 - lambda * exp(-lambda)) * (1 - PI) * exp(-lambda) / (prob ^ 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    ifelse(x == 1, PI * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)),
           (1 - PI) * (lambda ^ x) * exp(-lambda) / ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) * factorial(x)))
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- invlink(eta)
    PI <- lambda[, 2]
    lambda <- lambda[, 1]
    CDF <- function(x) {
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, (1 - PI) * exp(-lambda) / (1 - lambda * exp(-lambda)), 
      PI + (1 - PI) * (stats::ppois(x, lambda) - lambda * exp(-lambda)
      ) / (1 - lambda * exp(-lambda)))))
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
    if (is.null(controlMethod$piStart)) {
      cc <- colnames(Xvlm)
      cc <- cc[grepl(x = cc, pattern = "pi$")]
      cc <- unlist(strsplit(x = cc, ":pi"))
      cc <- sapply(cc, FUN = function(x) {
        ifelse(x %in% names(start), start[x], 0) # TODO: gosh this is terrible pick a better method
      })
      start <- c(start, cc)
    } else {
      start <- c(start, controlMethod$piStart)
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
      popVar= popVar,
      family = "Hurdleztpoisson",
      parNum = 2,
      etaNames = c("lambda", "pi"),
      densityFunction = dFun,
      simulate = simulate,
      getStart = getStart
    ),
    class = "family"
  )
}
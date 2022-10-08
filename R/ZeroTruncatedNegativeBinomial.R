#' Zero truncated Negative Binomial model
#'
#' @param nSim TODO
#' @param epsSim TODO
#' @param ... TODO
#' @return A object of class "family" containing objects \cr
#' makeMinusLogLike(y,X) - for creating negative likelihood function \cr
#' makeGradient(y,X) - for creating gradient function \cr
#' makeHessian(X) - for creating hessian \cr
#' linkfun - a link function to connect between linear predictor and model parameter in regression and a name of link function\cr
#' linkinv - an inverse function of link \cr
#' Dlink - a 1st derivative of link function \cr
#' mu.eta,Variance - Expected Value and Variance \cr
#' valedmu, valideta - for checking if regression arguments and valid\cr
#' family - family name\cr
#' Where: \cr
#' y is a vector of observed values \cr
#' X is a matrix / data frame of covariates
#' @importFrom stats uniroot
#' @importFrom stats dnbinom
#' @export
ztnegbin <- function(nSim = 1000, epsSim = 1e-8, ...) {
  # Fist for lambda second for alpha
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2])), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),exp(x[,2])), ncol = 2, dimnames = dimnames(x))}

  mu.eta <- function(eta, type = "trunc", ...) {
    A <- invlink(eta)
    lambda <- A[, 1]
    A <- A[, 2]
    switch (type,
      nontrunc = lambda,
      trunc = lambda / (1 - (1 + A * lambda) ** (-1 / A))
    )
  }

  variance <- function(eta, type = "nontrunc", ...) {
    A <- invlink(eta)
    lambda <- A[, 1]
    A <- A[, 2]
    switch (type,
      nontrunc = lambda * (1 + A * lambda),
      trunc = (lambda + A * (lambda ** 2) - A * (lambda ** 2) * ((1 + A * lambda) ** (-1 / A))) / ((1 - (1 + A * lambda) ** (-1 / A)) ** 2)
    )
  }
  
  compdigamma <- function(y, alpha) {
    #temp <- 0:(y-1)
    #sum(-(alpha ** (-2)) / (temp + 1 / alpha))
    (-digamma(y + 1 / alpha) + digamma(1 / alpha)) / (alpha ** 2)
  }
  
  
  comptrigamma <- function(y, alpha) {
    #temp <- 0:(y-1)
    #sum(-(temp ** 2) / (((1 + temp * alpha)) ** 2))
    (alpha ** 2 * (1 - y) + 2 * alpha * digamma(y + 1 / alpha) + trigamma(y + 1 / alpha) - 2 * alpha * digamma(1 + 1 / alpha) - trigamma(1 + 1 / alpha)) / (alpha ** 4)
  }
  
  # Computing the expected value of di/trigamma functions on (y + 1/alpha)
  
  compExpect <- function(eta) {
    alpha <- invlink(matrix(eta, ncol = 2))
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    res <- res1 <- 0
    k <- 0
    repeat{
      k <- k + 1 # 1 is the first possible y value for 0 truncated distribution
      prob <- stats::dnbinom(x = k, size = 1 / alpha, mu = lambda) / (1 - stats::dnbinom(x = 0, size = 1 / alpha, mu = lambda))
      toAdd <- compdigamma(y = k, alpha = alpha) * prob
      toAdd1 <- comptrigamma(y = k, alpha = alpha) * prob
      res <- res + toAdd
      res1 <- res1 + toAdd1
      if ((k == nSim) | ((abs(toAdd) < epsSim) & (abs(toAdd1) < epsSim))) {
        break
      }
    }
    c(res, res1)
  }
  
  Wfun <- function(prior, eta, ...) {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    z <- 1 / alpha
    M <- ((1 + lambda / z) ** z) - 1
    S <- 1 / (1 + lambda / z)
    G <- 1 / (1 - (1 + lambda / z) ** (-z))
    Ey <- mu.eta(eta = eta)
    cp4 <- z ** 2
    cp9 <- (1 / S)
    cp1 <- log(cp9) * cp4
    cp2 <- lambda * S
    cp3 <- S ** (-z)
    cp5 <- log(S) / S
    cp6 <- (S ** 2)
    cp7 <- log(S)
    cp8 <- (lambda / z + cp5)
    cp10 <- cp9 ** (-1 - z)
    cp11 <- (S / M)
    T1 <- lambda * (Ey - lambda) * cp6
    T2 <- cp10 * (lambda * (1 + z) * S + cp7 * cp4)
    T2 <- G * T2
    T3 <- cp9 ** (-z)
    T3 <- T3 * (-cp7 * cp4 - z * cp2)
    T3 <- (G ** 2) * (-T3) * cp10
    C1 <- as.numeric(((cp9 ** z) * (lambda - 1) + 1) *
                       cp6 / ((cp9 ** z - 1) ** 2))
    C2 <- (1 + Ey / z) * cp6
    Edig <- apply(X = eta, MARGIN = 1, FUN = function(x) {compExpect(x)})
    Etrig <- Edig[2,]
    Edig <- Edig[1,]
    matrix(
      c(-lambda * (C1 - C2) * prior,           # lambda predictor derivative without X matrix,
        -(-T1 + lambda * (T2 + T3)) * prior * alpha, # mixed derivative without X matrix
        -(-T1 + lambda * (T2 + T3)) * prior * alpha, # mixed derivative without X matrix
        -((2 * cp2 * cp4 + 2 * cp7 * (z ** 3) + Etrig + # alpha predictor derivative without X matrix
        (Ey + z) * (lambda ** 2) * cp6 + (z ** 3) * 2 * cp8 * cp11 +
        cp4 * (S ** (1 - z)) * (z * cp2 + cp7 * cp4) *
        cp8 / (M ** 2) + cp4 * lambda * log(cp9) * cp11 +
        cp4 * lambda * cp6 * cp8 / M) * (alpha ** 2)+
        (cp1 + Ey * z - (Ey + z) * cp2 + Edig +
        G * (1 / cp3) * (cp1 - z * cp2)) * alpha) * prior
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "alpha")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    z <- 1 / alpha
    S <- 1 / (1 + lambda / z)
    G <- 1 / (1 - (1 + lambda / z) ** (-z))
    cp1 <- log(1 / S) * (z ** 2)
    cp2 <- lambda * S
    cp3 <- S ** (-z)
    dig <- compdigamma(y = y, alpha = alpha)

    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    uMatrix <- matrix(
      c((y + (lambda - y) * cp3) * S / (1 - cp3),
        (cp1 + y * z - (y + z) * cp2 + dig +
        G * (1 / cp3) * (cp1 - z * cp2)) * alpha),
      ncol = 2)
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      #xx <- solve(weight[[x]])
      xx %*% uMatrix[x, ]
    })
    pseudoResid <- t(pseudoResid)
    dimnames(pseudoResid) <- dimnames(eta)
    pseudoResid
  }

  minusLogLike <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(beta) {
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      alpha <- invlink(eta)
      lambda <- alpha[, 1]
      alpha <- alpha[, 2]
      z <- 1 / alpha
      M <- 1 + lambda / z

      -sum(weight * (lgamma(y + z) - lgamma(z) -
      log(factorial(y)) - (y + z) * log(M) +
      y * log(lambda / z) - log(1 - (M ** (-z)))))
    }
  }


  gradient <- function(y, X, weight = 1, NbyK = FALSE, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(beta) {
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      alpha <- invlink(eta)
      lambda <- alpha[, 1]
      alpha <- alpha[, 2]
      # z is inverse of alpha in documentation
      z <- 1 / alpha
      S <- 1 / (1 + lambda / z)
      G <- 1 / (1 - (1 + lambda / z) ** (-z))
      cp1 <- log(1 / S) * (z ** 2)
      cp2 <- lambda * S
      cp3 <- S ** (-z)

      # log(alpha) derivative
      dig <- compdigamma(y = y, alpha = alpha)
      G0 <- t((cp1 + y * z - (y + z) * cp2 + dig +
              G * (1 / cp3) * (cp1 - z * cp2)) * alpha)

      # Beta derivative
      G1 <- t(((y + (lambda - y) * cp3) * S / (1 - cp3))  * weight)
      
      if (NbyK) {
        XX <- sapply(as.data.frame(X[1:nrow(eta), ]), FUN = function(x) {all(x == 0)})
        return(cbind(as.data.frame(X[1:nrow(eta), !(XX)]) * G1, as.data.frame(X[-(1:nrow(eta)), XX]) * G0))
      }

      as.numeric(c(G1, G0) %*% X)
    }
  }

  hessian <- function(y, X, weight = 1, lambdaPredNumber, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(beta) {
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      alpha <- invlink(eta)
      lambda <- alpha[, 1]
      alpha <- alpha[, 2]
      Xlambda <- X[1:(nrow(X) / 2), 1:lambdaPredNumber]
      Xalpha <- X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
      # z is inverse of alpha in documentation
      z <- 1 / alpha
      M <- ((1 + lambda / z) ** z) - 1
      S <- 1 / (1 + lambda / z)
      G <- 1 / (1 - (1 + lambda / z) ** (-z))
      res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
      cp4 <- z ** 2
      cp9 <- (1 / S)
      cp1 <- log(cp9) * cp4
      cp2 <- lambda * S
      cp3 <- S ** (-z)
      cp5 <- log(S) / S
      cp6 <- (S ** 2)
      cp7 <- log(S)
      cp8 <- (lambda / z + cp5)
      cp10 <- cp9 ** (-1 - z)
      cp11 <- (S / M)
      
      trig <- comptrigamma(y = y, alpha = alpha)
      dig <- compdigamma(y = y, alpha = alpha)
      
      # 2nd log(alpha) derivative
       G00 <- t(as.data.frame(Xalpha) *
              ((2 * cp2 * cp4 + 2 * cp7 * (z ** 3) + trig +
              (y + z) * (lambda ** 2) * cp6 + (z ** 3) * 2 * cp8 * cp11 +
              cp4 * (S ** (1 - z)) * (z * cp2 + cp7 * cp4) *
              cp8 / (M ** 2) + cp4 * lambda * log(cp9) * cp11 +
              cp4 * lambda * cp6 * cp8 / M) * (alpha ** 2) +
              (cp1 + y * z - (y + z) * cp2 + dig +
              G * (1 / cp3) * (cp1 - z * cp2)) * alpha) * weight) %*% Xalpha
      # mixed derivative
      T1 <- lambda * (y - lambda) * cp6

      T2 <- cp10 * (lambda * (1 + z) * S + cp7 * cp4)
      T2 <- G * T2

      T3 <- cp9 ** (-z)
      T3 <- T3 * (-cp7 * cp4 - z * cp2)
      T3 <- (G ** 2) * (-T3) * cp10

      G01 <- t(as.data.frame(Xlambda) * as.numeric(-T1 + lambda * (T2 + T3)) * alpha * weight) %*% as.matrix(Xalpha)

      # second beta derivative
      C1 <- as.numeric(((cp9 ** z) * (lambda - 1) + 1) *
                       cp6 / ((cp9 ** z - 1) ** 2))
      C2 <- (1 + y / z) * cp6

      G11 <- t(as.data.frame(Xlambda) * lambda * (C1 - C2) * weight) %*% Xlambda

      res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
      res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
      res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(G01)
      res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- G01

      res
    }
  }

  validmu <- function(mu) {
    all(is.finite(mu)) && all(0 < mu)
  }

  dev.resids <- function (y, eta, wt, ...) {
    disp1 <- invlink(eta)
    mu <- disp1[, 1]
    disp1 <- disp1[, 2]
    mu1 <- mu.eta(eta = eta)
    hm1y <- y
    hm1y[y == 1] <- -16
    a <- function(n) {stats::uniroot(f = function(x) {mu.eta(matrix(c(x, eta[n, 2]), ncol = 2)) - y[n]}, lower = -log(y[n]), upper = y[n] * 10, tol = .Machine$double.eps)$root}
    hm1y[y > 1] <- sapply(which(y > 1), FUN = a)
    loghm1ytdisp <- log(disp1 * exp(hm1y))
    logprobhm1y <- log(1 - ((1 + disp1 * exp(hm1y)) ** (-1 / disp1)))
    sign(y - mu1) * sqrt(-2 * wt * (-(y + 1 / disp1) * log(1 + mu * disp1) + y * log(mu * disp1) - log(1 - ((1 + mu * disp1) ** (-1/disp1))) +
                                    (y + 1 / disp1) * log(1 + exp(hm1y) * disp1) - y * loghm1ytdisp + logprobhm1y))
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    disp <- invlink(eta)
    lambda <- disp[, 1]
    disp <- disp[, 2]
    pr <- 1 - (1 + disp * lambda) ** (- 1 / disp)
    N <- pw / pr
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    z <- invlink(eta)
    lambda <- z[, 1]
    z <- z[, 2]
    pr <- 1 - (1 + z * lambda) ** (- 1 / z)
    S <- 1 / (1 + z * lambda)

    cp1 <- (1 / S)
    cp2 <- (S ** (1 - 1 / z))
    cp3 <- ((1 - cp1 ** (1 / z)) ** 2)

    bigTheta1 <- pw * z * (cp2 * (cp1 * log(cp1) - z * lambda) / ((z ** 2) * cp3)) # w.r to alpha
    bigTheta2 <- -(pw * as.numeric(lambda * cp2 / cp3)) # w.r to lambda

    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)

    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * (1 - pr) / (pr ** 2))

    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / (1 - stats::dnbinom(x = 0, mu = lambda, size = 1  / alpha))
  }

  simulate <- function(n, eta, lower = 0, upper = Inf) {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    lb <- stats::pnbinom(lower, mu = lambda, size = 1 / alpha)
    ub <- stats::pnbinom(upper, mu = lambda, size = 1 / alpha)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qnbinom(p_u, mu = lambda, size = 1 / alpha)
    sims
  }
  
  structure(
    list(
      makeMinusLogLike = minusLogLike,
      makeGradient = gradient,
      makeHessian = hessian,
      linkfun = link,
      linkinv = invlink,
      mu.eta = mu.eta,
      link = c("log", "log"),
      valideta = function (eta) {TRUE},
      variance = variance,
      Wfun = Wfun,
      funcZ = funcZ,
      dev.resids = dev.resids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      simulate = simulate,
      family = "ztnegbin",
      parNum = 2,
      etaNames = c("lambda", "alpha"),
      densityFunction = dFun
    ),
    class = "family"
  )
}

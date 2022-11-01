#' Zero One Truncated Negative Binomial model
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
#' @export
zotnegbin <- function(nSim = 1000, epsSim = 1e-8, ...) {
  # Fist for lambda second for alpha
  link <- function (x) {matrix(c(log(x[,1]),log(x[,2])), ncol = 2, dimnames = dimnames(x))}
  invlink <- function (x) {matrix(c(exp(x[,1]),exp(x[,2])), ncol = 2, dimnames = dimnames(x))}

  mu.eta <- function(eta, type = "trunc", ...) {
    A <- invlink(eta)
    lambda <- A[, 1]
    A <- A[, 2]
    switch (type,
            "nontrunc" = lambda,
            "trunc" = (lambda - lambda * ((1 + A * lambda) ** (-1 - 1 / A))) / (1 - (1 + A * lambda) ** (-1 / A) - lambda * ((1 + A * lambda) ** (-1 - 1 / A)))
    )
  }

  variance <- function(eta, type = "nontrunc", ...) {
    A <- invlink(eta)
    lambda <- A[, 1]
    A <- A[, 2]
    switch (type,
            "nontrunc" = lambda * (1 + A * lambda),
            "trunc" = {
              EX2 <- lambda * (1 + A * lambda) + lambda ** 2;
              prob <- (1 - (1 + A * lambda) ** (-1 / A) - lambda * ((1 + A * lambda) ** (-1 - 1 / A)));
              term <- lambda * ((1 + A * lambda) ** (-1 - 1 / A));
              EX <- (lambda - lambda * ((1 + A * lambda) ** (-1 - 1 / A))) / (1 - (1 + A * lambda) ** (-1 / A) - lambda * ((1 + A * lambda) ** (-1 - 1 / A)));
              (EX2 - term) / prob - EX ** 2
            }
    )
  }
  
  compExpect <- function (eta) {
    alpha <- invlink(matrix(eta, ncol = 2))
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    z <- 1 / alpha
    res <- res1 <- 0
    k <- 1
    repeat{
      k <- k + 1 # 2 is the first possible y value for 0 truncated distribution
      prob <- stats::dnbinom(x = k, size = 1 / alpha, mu = lambda) / (1 - stats::dnbinom(x = 0, size = 1 / alpha, mu = lambda) - stats::dnbinom(x = 1, size = 1 / alpha, mu = lambda))
      toAdd <- (digamma(k + z) - digamma(z)) * z * prob
      toAdd1 <- (trigamma(k + z) - trigamma(z)) * (z ** 2) * prob
      res <- res + toAdd
      res1 <- res1 + toAdd1
      if ((k == nSim) | ((abs(toAdd) < 1e-8) & (abs(toAdd1) < 1e-8))) {
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
    M <- 1 + lambda * alpha
    S <- 1 / M
    Ey <- mu.eta(eta = eta)
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    # The following vectors are meant to decrease compitation time
    cp1 <- M ** z
    cp2 <- alpha * lambda
    cp3 <- (2 * alpha + 1)
    #cp4 <- (y + lambda) No longer needed
    cp5 <- log(M)
    #cp6 <- (y - lambda)
    cp7 <- (1 + alpha)
    cp8 <- M ** (1 + z)
    cp9 <- lambda * alpha
    cp10 <- (M ** (2 + z))
    cp11 <- (alpha ** 2)
    #cp12 <- (1 + lambda)
    #cp13 <- y * alpha
    cp14 <- (lambda ** 2)
    cp15 <- lambda * S
    cp16 <- M ** 2
    cp17 <- (-1 - z)
    cp18 <- (cp8 - cp2 - lambda - 1)
    cp19 <- (cp1 - 1)
    cp20 <- cp19 ** 2
    cp21 <- lambda * (Ey + lambda) * cp11
    cp22 <- (Ey - lambda) * alpha * cp8
    cp23 <- (1 + lambda) * Ey * alpha
    cp24 <- (cp5 * cp10 - cp21 + cp22 - cp23)
    cp25 <- (1 - z * cp7)
    cp26 <- z * cp3
    cp27 <- 1 / cp10
    Edig <- apply(X = eta, MARGIN = 1, FUN = function(x) {compExpect(x)})
    Etrig <- Edig[2,]
    Edig <- Edig[1,]
    
    G11 <- lambda * ((-(alpha ** 3) * Ey * cp14 * cp20 + cp11 * lambda * (cp14 * cp1 - lambda * cp19 *
    (1 - 2 * Ey + cp1) - 2 * Ey * cp20) + cp14 * cp1 + alpha * 
    ((lambda ** 3) * cp1 + cp14 * (1 - Ey + cp1) - 2 * lambda * cp19 *
    (cp1 - Ey) - Ey * cp20) - cp20) / (cp16 * ((cp1 + lambda * (alpha * cp19 - 1) - 1) ** 2))) * prior
    
    G01 <- lambda * (-alpha * (Ey + z) * S + lambda * cp11 * (Ey + z) * (1 / cp16) +
    (cp17 * cp2 * cp27) / prob + (lambda * cp27) / prob + S +
    cp2 * cp17 * cp27 * (cp2 * (-2 - z) * S + cp5 * z) / prob -
    (alpha * cp17 * lambda * cp27 * ((1 / cp1) * (z * log(S) + cp15) - 
    lambda * (1 / cp8) * (z * cp5 + S * cp9 * cp17))) / (prob ** 2)) * prior
    
    G00 <- (Etrig + Edig + z * (((M ** cp26) * cp5 * (cp5 * (2 - cp26) + cp15 * cp3)) +
    cp22 * (cp5 * cp25 + cp15 * cp7) - 2 * cp21 + cp9 * (M ** (cp26 - 1)) +
    cp22 - cp23) / (M * cp18) - z * cp24 * (cp8 * (cp25 * cp5 + cp15 * cp7) - cp9) /
    (M * (cp18 ** 2)) - z * (cp5 * cp10 - cp21 + cp22 - cp23) /
    (M * cp18) - lambda * cp24 / (cp16 * cp18)) * prior
    
    matrix(
      -c(G11, # lambda predictor derivative without X matrix,
        G01,  # mixed derivative without X matrix
        G01,  # mixed derivative without X matrix
        G00   # alpha predictor derivative without X matrix
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
    M <- 1 + lambda * alpha
    S <- 1 / M
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    cp1 <- z * log(M)
    cp2 <- lambda * S
    cp3 <- alpha * (-1 - z)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    uMatrix <- matrix(
      c((y - alpha * (y + z) * cp2 + lambda * cp3 * lambda * (S ** (2 + z)) / prob),
        (-(digamma(y + z) - digamma(z)) * (z ** 2) - 
        z * ((S ** z) * (z * log(S) + cp2) - lambda *
        (S ** (1 + z)) * (cp1 + cp3 * cp2)) / prob +
        cp1 * z - cp2 * (y + z) + y * z) * alpha),
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
      S <- 1 / (1 + lambda / z)
      prob <- 1 - S ** z - lambda * (S ** (1 + z))

      -sum(weight * (lgamma(y + z) - lgamma(z) -
      log(factorial(y)) - (y + z) * log(1 + lambda / z) +
      y * log(lambda / z) - log(prob)))
    }
  }


  gradient <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, ...) {
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
      M <- 1 + lambda * alpha
      S <- 1 / M
      prob <- 1 - S ** z - lambda * (S ** (1 + z))
      cp1 <- z * log(M)
      cp2 <- lambda * S
      cp3 <- alpha * (-1 - z)

      # log(alpha) derivative
      G0 <- t((-(digamma(y + z) - digamma(z)) * (z ** 2)  - 
            z * ((S ** z) * (z * log(S) + cp2) - lambda *
            (S ** (1 + z)) * (cp1 + cp3 * cp2)) / prob +
            cp1 * z - cp2 * (y + z) + y * z) * alpha * weight)
      # Beta derivative
      G1 <- (weight * (y - alpha * (y + z) * cp2 + lambda * cp3 * lambda * (S ** (2 + z)) / prob))
      
      if (NbyK) {
        XX <- sapply(as.data.frame(X[1:nrow(eta), ]), FUN = function(x) {all(x == 0)})
        return(cbind(as.data.frame(X[1:nrow(eta), !(XX)]) * G1, as.data.frame(X[-(1:nrow(eta)), XX]) * G0))
      }
      if (vectorDer) {
        return(cbind(G1, G0))
      }

      as.numeric(c(G1, G0) %*% X)
    }
  }

  hessian <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)

    function(beta) {
      lambdaPredNumber <- attr(X, "hwm")[1]
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      alpha <- invlink(eta)
      lambda <- alpha[, 1]
      alpha <- alpha[, 2]
      Xlambda <- X[1:(nrow(X) / 2), 1:lambdaPredNumber]
      Xalpha <- X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]
      # z is inverse of alpha in documentation
      z <- 1 / alpha
      M <- 1 + lambda * alpha
      S <- 1 / M
      prob <- 1 - S ** z - lambda * (S ** (1 + z))
      res <- matrix(nrow = length(beta), ncol = length(beta), dimnames = list(names(beta), names(beta)))
      # The following vectors are meant to decrease compitation time
      cp1 <- M ** z
      cp2 <- alpha * lambda
      cp3 <- (2 * alpha + 1)
      #cp4 <- (y + lambda) No longer needed
      cp5 <- log(M)
      #cp6 <- (y - lambda)
      cp7 <- (1 + alpha)
      cp8 <- M ** (1 + z)
      cp9 <- lambda * alpha
      cp10 <- (M ** (2 + z))
      cp11 <- (alpha ** 2)
      #cp12 <- (1 + lambda)
      #cp13 <- y * alpha
      cp14 <- (lambda ** 2)
      cp15 <- lambda * S
      cp16 <- M ** 2
      cp17 <- (-1 - z)
      cp18 <- (cp8 - cp2 - lambda - 1)
      cp19 <- (cp1 - 1)
      cp20 <- cp19 ** 2
      cp21 <- lambda * (y + lambda) * cp11
      cp22 <- (y - lambda) * alpha * cp8
      cp23 <- (1 + lambda) * y * alpha
      cp24 <- (cp5 * cp10 - cp21 + cp22 - cp23)
      cp25 <- (1 - z * cp7)
      cp26 <- z * cp3
      cp27 <- 1 / cp10

      # 2nd log(alpha) derivative

      G00 <- t(as.data.frame(Xalpha) * weight *
             ((trigamma(y + z) - trigamma(z)) * (z ** 2) + 
             (digamma(y + z) - digamma(z)) * z + z * 
             (((M ** cp26) * cp5 * (cp5 * (2 - cp26) + cp15 * cp3)) +
             cp22 * (cp5 * cp25 + cp15 * cp7) - 2 * cp21 + cp9 * (M ** (cp26 - 1)) +
             cp22 - cp23) / (M * cp18) - z * cp24 * (cp8 * (cp25 * cp5 + cp15 * cp7) - cp9) /
             (M * (cp18 ** 2)) - z * (cp5 * cp10 - cp21 + cp22 - cp23) /
             (M * cp18) - lambda * cp24 / (cp16 * cp18))) %*% Xalpha

      # mixed derivative
      term <- (-alpha * (y + z) * S + lambda * cp11 * (y + z) * (1 / cp16) +
               (cp17 * cp2 * cp27) / prob + (lambda * cp27) / prob + S +
               cp2 * cp17 * cp27 * (cp2 * (-2 - z) * S + cp5 * z) / prob -
               (alpha * cp17 * lambda * cp27 *
               ((1 / cp1) * (z * log(S) + cp15) - lambda * (1 / cp8) *
               (z * cp5 + S * cp9 * cp17))) / (prob ** 2))

      G01 <- t(as.data.frame(Xlambda) * as.numeric(term * lambda * weight)) %*% as.matrix(Xalpha)
      # second beta derivative
      term <- ((-(alpha ** 3) * y * cp14 * cp20 +
              cp11 * lambda * (cp14 * cp1 - lambda * cp19 *
              (1 - 2 * y + cp1) - 2 * y * cp20) +
              cp14 * cp1 + alpha * ((lambda ** 3) * cp1 +
              cp14 * (1 - y + cp1) - 2 * lambda * cp19 *
              (cp1 - y) - y * cp20) - cp20) /
              (cp16 * ((cp1 + lambda * (alpha * cp19 - 1) - 1) ** 2)))

      G11 <- t(as.data.frame(Xlambda) * lambda * weight * as.numeric(term)) %*% Xlambda


      res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
      res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
      res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(G01)
      res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- G01
      
      res
    }
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  dev.resids <- function (y, eta, wt, ...) {
    # disp1 <- exp(eta[, 2])
    # mu <- invlink(eta[, 1])
    # mu1 <- mu.eta(eta = eta)
    # a <- function(n) {stats::uniroot(f = function(x) {mu.eta(matrix(c(x, eta[n, 2]), ncol = 2)) - y[n]}, lower = -log(y[n]), upper = y[n] * 10, tol = .Machine$double.eps)$root}
    # loghm1y <- y
    # loghm1y[y == 2] <- -11
    # loghm1y[y > 2] <- sapply(y[y > 2], FUN = a)
    # h1my <- exp(loghm1y)
    # logh1mydivdisp1 <- ifelse(y > 2, log(h1my / disp1), 0)
    # logprobdey <- ifelse(y > 2, log(1 - (1 + disp1 * h1my) ** (-1 / disp1) - h1my * ((1 + disp1 * h1my) ** (- 1 - 1 / disp1))), 0)
    # sign(y - mu1) * sqrt(-2 * wt * (-(y + 1 / disp1) * log(1 + mu * disp1) + y * log(mu / disp1) - log(1 - (1 + disp1 * mu) ** (-1 / disp1) - mu * ((1 + disp1 * mu) ** (- 1 - 1 / disp1))) +
    #                                  (y + 1 / disp1) * log(1 + h1my * disp1) - y * logh1mydivdisp1 + logprobdey))
    NULL
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    z <- invlink(eta)
    lambda <- z[, 1]
    z <- 1 / z[, 2]
    S <- 1 / (1 + lambda / z)
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    N <- (pw * (1 - lambda * (S ** (1 + z))) / prob)
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    z <- 1 / alpha
    M <- 1 + lambda / z
    S <- 1 / M
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    
    bigTheta1 <- (pw * alpha *  as.numeric(
                  (prob * lambda * (S ** (2 + z)) * 
                  (alpha * (alpha + 1) * lambda -M * log(M)) -
                  (1 - lambda * (S ** (1 + z))) *(-lambda * (S ** (1 + z)) *
                  (log(M) * (z ** 2) - (1 + z) * lambda * S * z) -
                  (S ** z) * (log(M) * (z ** 2) - lambda * z * S))) /
                  (prob ** 2)))
    
    bigTheta2 <- (pw * as.numeric(lambda *
                  (prob * (lambda - 1) * (S ** (2 + z)) -
                  (1 + alpha) * lambda * (S ** (2 + z)) *
                  (1 - lambda * (S ** (1 + z)))) / (prob ** 2)))
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * ((1 - lambda * (S ** (1 + z))) ** 2) *
              (1 - prob) / (prob ** 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    alpha <- invlink(eta)
    lambda <- alpha[, 1]
    alpha <- alpha[, 2]
    stats::dnbinom(x = x, mu = lambda, size = 1 / alpha) / (1 - stats::dnbinom(x = 0, mu = lambda, size = 1 / alpha) - stats::dnbinom(x = 1, mu = lambda, size = 1 / alpha))
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
      popVar = popVar,
      simulate = simulate,
      family = "zotnegbin",
      parNum = 2,
      etaNames = c("lambda", "alpha"),
      densityFunction = dFun
    ),
    class = "family"
  )
}

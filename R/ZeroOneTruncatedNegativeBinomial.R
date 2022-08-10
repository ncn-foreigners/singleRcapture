#' Zero One Truncated Negative Binomial model
#'
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
zotnegbin <- function() {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }

  mu.eta <- function(eta, type = "trunc", ...) {
    A <- exp(eta[, 2])
    lambda <- invlink(eta[, 1])
    switch (type,
            "nontrunc" = lambda,
            "trunc" = (lambda - lambda * ((1 + A * lambda) ** (-1 - 1 / A))) / (1 - ((1 + A * lambda) ** (-1 / A) + lambda * ((1 + A * lambda) ** (-1 - 1 / A))))
    )
  }

  variance <- function(mu, disp, type = "nontrunc", ...) {
    A <- exp(disp)
    switch (type,
            "nontrunc" = mu * (1 + A * mu),
            "trunc" = (mu - mu * exp(-mu)) / (1 - exp(-mu) - mu * exp(-mu))
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
  
  compExpect <- function (eta1, eta2) {
    lambda <- exp(eta1)
    alpha <- exp(eta2)
    res <- res1 <- 0
    k <- 0
    repeat{
      k <- k + 1 # 1 is the first possible y value for 0 truncated distribution
      prob <- stats::dnbinom(size = 1 / alpha, mu = lambda, x = k) / (1 - stats::dnbinom(x = 0, size = 1 / alpha, mu = lambda))
      toAdd <- compdigamma(y = k, alpha = alpha) * prob
      toAdd1 <- comptrigamma(y = k, alpha = alpha) * prob
      res <- res + toAdd
      res1 <- res1 + toAdd1
      if ((k == nSim) | ((abs(toAdd) < 1e-8) & (abs(toAdd1) < 1e-8))) {
        break
      }
    }
    c(res, res1)
  }
  
  Wfun <- function(prior, eta, disp, ...) {
    alpha <- exp(disp)
    z <- 1 / alpha
    lambda <- exp(eta)
    yexp <- mu.eta(eta = eta, disp = disp)
    M <- 1 + lambda * alpha
    S <- 1 / M
    cp1 <- M ** z
    cp11 <- (alpha ** 2)
    cp14 <- (lambda ** 2)
    cp16 <- M ** 2
    cp19 <- (cp1 - 1)
    cp20 <- cp19 ** 2
    term <- ((-(alpha ** 3) * yexp * cp14 * cp20 +cp11 * lambda * (cp14 * cp1 - lambda * cp19 * (1 - 2 * yexp + cp1) - 2 * yexp * cp20) +
                cp14 * cp1 + alpha * ((lambda ** 3) * cp1 + cp14 * (1 - yexp + cp1) - 2 * lambda * cp19 * (cp1 - yexp) - yexp * cp20) - cp20) /
               (cp16 * ((cp1 + lambda * (alpha * cp19 - 1) - 1) ** 2)))
    -lambda * term
  }
  
  funcZ <- function(eta, weight, disp, y, mu, ...) {
    alpha <- exp(disp)
    z <- 1 / alpha
    lambda <- exp(eta)
    M <- 1 + lambda * alpha
    S <- 1 / M
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    cp1 <- z * log(M)
    cp2 <- lambda * S
    cp3 <- alpha * (-1 - z)
    dig <- compdigamma(y = y, alpha = alpha)
    c((y - alpha * (y + z) * cp2 + lambda * cp3 * lambda * (S ** (2 + z)) / prob),
      (dig + cp1 - alpha * (y + z) * cp2 -
      ((S ** z) * (z * log(S) + cp2) - lambda *
      (S ** (1 + z)) * (cp1 + cp3 * cp2)) / prob +
      y)) / as.numeric(weight[,1:2])
  }

  minusLogLike <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(beta) {
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      lambda <- exp(eta[,1])
      alpha <- exp(eta[,2])
      z <- 1 / alpha
      S <- 1 / (1 + lambda / z)
      prob <- 1 - S ** z - lambda * (S ** (1 + z))

      -sum(weight * (lgamma(y + z) - lgamma(z) -
      log(factorial(y)) - (y + z) * log(1 + lambda / z) +
      y * log(lambda / z) - log(prob)))
    }
  }


  gradient <- function(y, X, weight = 1, ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(beta) {
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      lambda <- exp(eta[,1])
      alpha <- exp(eta[,2])
      # z is inverse of alpha in documentation
      z <- 1 / alpha
      M <- 1 + lambda * alpha
      S <- 1 / M
      prob <- 1 - S ** z - lambda * (S ** (1 + z))
      cp1 <- z * log(M)
      cp2 <- lambda * S
      cp3 <- alpha * (-1 - z)
      dig <- compdigamma(y = y, alpha = alpha)

      # log(alpha) derivative
      G0 <- t((dig - (-lambda * (S ** (1 + z)) * (log(M) * (z ** 2) - (1 + z) * lambda * S) -
            (S ** z) * (log(M) * (z ** 2) - cp2 * z)) / prob +
            log(M) * (z ** 2) - S * lambda * (y + z) + y * z) * alpha * weight)
      # Beta derivative
      G1 <- (weight * (y - alpha * (y + z) * cp2 + lambda * cp3 * lambda * (S ** (2 + z)) / prob))

      as.numeric(c(G1, G0) %*% X)
    }
  }

  hessian <- function(y, X, weight = 1, lambdaPredNumber) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(beta) {
      eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
      lambda <- exp(eta[,1])
      alpha <- exp(eta[,2])
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
      #cp18 <- (cp8 - cp2 - lambda - 1)
      cp19 <- (cp1 - 1)
      cp20 <- cp19 ** 2
      #cp21 <- lambda * (y + lambda) * cp11
      #cp22 <- (y - lambda) * alpha * cp8
      #cp23 <- (1 + lambda) * y * alpha
      #cp24 <- (cp5 * cp10 - cp21 + cp22 - cp23)
      #cp25 <- (1 - z * cp7)
      #cp26 <- z * cp3
      cp27 <- 1 / cp10
      dig <- compdigamma(y = y, alpha = alpha)
      trig <- comptrigamma(y = y, alpha = alpha)

      # 2nd log(alpha) derivative
      # TODO this is bad repair

      G00 <- t(as.data.frame(Xalpha) *
             (trig * (alpha ** 2) + dig * alpha + (cp9 ** 2) * (y + z) * (S ** 2) +
             (1 / cp1) * ((cp15 ** 2) * cp11 + cp15 - z * cp5) - lambda / cp10 *
             (cp11 * (1 + z) * cp14 / cp16 - cp7 * cp15 + 2 * cp15 - z * cp5) + 
             (-((-cp15 + z * cp5) ** 2)) / cp1 - lambda * ((-cp15 * cp7 + z * cp5) ** 2) / 
             (cp8 * prob) - alpha * (y + z) * cp15 + 2 * cp15 + 
             (((1 / cp1) * (cp15 - z * cp5) - (lambda / cp8) * (-cp15 * cp7 - z * cp5)) ** 2) / 
             (prob ** 2) - z * cp5)) %*% Xalpha

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

  pointEst <- function (disp, pw, lambda, contr = FALSE, ...) {
    z <- exp(-disp)
    S <- 1 / (1 + lambda / z)
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    N <- (pw * (1 - lambda * (S ** (1 + z))) / prob)
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, lambda, disp, cov, X, ...) {
    alpha <- exp(disp)
    z <- 1 / alpha
    M <- 1 + lambda / z
    S <- 1 / M
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    
    bigTheta1 <- sum(pw * alpha *  as.numeric(
                    (prob * lambda * (S ** (2 + z)) * 
                    (alpha * (alpha + 1) * lambda -M * log(M)) -
                    (1 - lambda * (S ** (1 + z))) *(-lambda * (S ** (1 + z)) *
                    (log(M) * (z ** 2) - (1 + z) * lambda * S * z) -
                    (S ** z) * (log(M) * (z ** 2) - lambda * z * S))) /
                    (prob ** 2)))
    
    bigTheta2 <- t(as.matrix(X)) %*% (pw * as.numeric(lambda *
                  (prob * (lambda - 1) * (S ** (2 + z)) -
                  (1 + alpha) * lambda * (S ** (2 + z)) *
                  (1 - lambda * (S ** (1 + z)))) / (prob ** 2)))
    
    bigTheta <- matrix(c(bigTheta1, bigTheta2), ncol = 1)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum(pw * ((1 - lambda * (S ** (1 + z))) ** 2) *
              (1 - prob) / (prob ** 2))
    
    f1 + f2
  }

  structure(
    list(
      makeMinusLogLike = minusLogLike,
      makeGradient = gradient,
      makeHessian = hessian,
      linkfun = link,
      linkinv = invlink,
      dlink = dlink,
      mu.eta = mu.eta,
      link = "log",
      valideta = function (eta) {TRUE},
      variance = variance,
      Wfun = Wfun,
      funcZ = funcZ,
      dev.resids = dev.resids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      family = "zotnegbin",
      parNum = 2
    ),
    class = "family"
  )
}

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

  mu.eta <- function(eta, disp, type = "trunc") {
    A <- exp(disp)
    lambda <- invlink(eta)
    switch (type,
            "nontrunc" = lambda,
            "trunc" = (lambda - lambda * ((1 + A * lambda) ** (-1 - 1 / A))) / (1 - ((1 + A * lambda) ** (-1 / A) + lambda * ((1 + A * lambda) ** (-1 - 1 / A))))
    )
  }

  variance <- function(mu, disp, type = "nontrunc") {
    A <- exp(disp)
    switch (type,
            "nontrunc" = mu * (1 + A * mu),
            "trunc" = (mu - mu * exp(-mu)) / (1 - exp(-mu) - mu * exp(-mu))
    )
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
    eta + ((y - alpha * (y + z) * cp2 + lambda * cp3 * lambda * (S ** (2 + z)) / prob)) / weight
  }

  minusLogLike <- function(y, X, weight = 1) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(arg) {
      alpha <- exp(arg[1])
      z <- 1 / alpha
      beta <- arg[-1]
      eta <- X %*% beta
      lambda <- exp(eta)
      S <- 1 / (1 + lambda / z)
      prob <- 1 - S ** z - lambda * (S ** (1 + z))

      -sum(weight * (lgamma(y + z) - lgamma(z) -
      log(factorial(y)) - (y + z) * log(1 + lambda / z)
      + y * log(lambda / z) - log(prob)))
    }
  }


  gradient <- function(y, X, weight = 1) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(arg) {
      alpha <- exp(arg[1])
      # z is inverse of alpha in documentation
      z <- 1 / alpha
      beta <- arg[-1]
      eta <- X %*% beta
      lambda <- exp(eta)
      M <- 1 + lambda * alpha
      S <- 1 / M
      prob <- 1 - S ** z - lambda * (S ** (1 + z))
      cp1 <- z * log(M)
      cp2 <- lambda * S
      cp3 <- alpha * (-1 - z)

      # log(alpha) derivative
      G0 <- sum((cp1 - alpha * (y + z) * cp2 -
                ((S ** z) * (z * log(S) + cp2) - lambda *
                (S ** (1 + z)) * (cp1 + cp3 * cp2)) / prob +
                y - (digamma(y + z) - digamma(z)) * z) * weight)
      # Beta derivative
      G1 <- t(X) %*% (weight * (y - alpha * (y + z) * cp2 +
                      lambda * cp3 * lambda *
                      (S ** (2 + z)) / prob))

      c(G0, G1)
    }
  }

  hessian <- function(y, X, weight = 1) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(arg) {
      alpha <- exp(arg[1])
      # z is inverse of alpha in documentation
      z <- 1 / alpha
      beta <- arg[-1]
      eta <- X %*% beta
      lambda <- exp(eta)
      M <- 1 + lambda * alpha
      S <- 1 / M
      prob <- 1 - S ** z - lambda * (S ** (1 + z))
      res <- matrix(1, nrow = length(arg), ncol = length(arg))
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

      G00 <- sum(weight * ((trigamma(y + z) - trigamma(z)) * (z ** 2) + 
            (digamma(y + z) - digamma(z)) * z +
             z * (((M ** cp26) * cp5 * (cp5 * (2 - cp26) + cp15 * cp3)) +
             cp22 * (cp5 * cp25 + cp15 * cp7) -
             2 * cp21 + cp9 * (M ** (cp26 - 1)) +
             cp22 - cp23) / (M * cp18) -
             z * cp24 * (cp8 * (cp25 * cp5 + cp15 * cp7) - cp9) /
             (M * (cp18 ** 2)) - z * (cp5 * cp10 - cp21 + cp22 - cp23) /
             (M * cp18) - lambda * cp24 / (cp16 * cp18)))

      # mixed derivative
      term <- (-alpha * (y + z) * S + lambda * cp11 * (y + z) * (1 / cp16) +
               (cp17 * cp2 * cp27) / prob + (lambda * cp27) / prob + S +
               cp2 * cp17 * cp27 * (cp2 * (-2 - z) * S + cp5 * z) / prob -
               (alpha * cp17 * lambda * cp27 *
               ((1 / cp1) * (z * log(S) + cp15) - lambda * (1 / cp8) *
               (z * cp5 + S * cp9 * cp17))) / (prob ** 2))

      G01 <- t(X) %*% as.numeric(term * lambda * weight)
      # second beta derivative
      term <- ((-(alpha ** 3) * y * cp14 * cp20 +
              cp11 * lambda * (cp14 * cp1 - lambda * cp19 *
              (1 - 2 * y + cp1) - 2 * y * cp20) +
              cp14 * cp1 + alpha * ((lambda ** 3) * cp1 +
              cp14 * (1 - y + cp1) - 2 * lambda * cp19 *
              (cp1 - y) - y * cp20) - cp20) /
              (cp16 * ((cp1 + lambda * (alpha * cp19 - 1) - 1) ** 2)))

      G11 <- t(as.data.frame(X) * lambda * weight * as.numeric(term)) %*% X


      res[1, 1] <- G00
      res[-1, -1] <- G11
      res[1, -1] <- G01
      res[-1, 1] <- G01

      res
    }
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  dev.resids <- function (y, mu, wt, disp = NULL) {
    eta <- log(mu)
    disp1 <- exp(disp)
    mu1 <- mu.eta(eta = eta, disp = disp)
    a <- function(y) {stats::uniroot(f = function(x) {mu.eta(x, disp = disp) - y}, lower = -log(y), upper = y * 10, tol = .Machine$double.eps)$root}
    loghm1y <- y
    loghm1y[y == 2] <- -11
    loghm1y[y > 2] <- sapply(y[y > 2], FUN = a)
    h1my <- exp(loghm1y)
    logh1mydivdisp1 <- ifelse(y > 2, log(h1my / disp1), 0)
    logprobdey <- ifelse(y > 2, log(1 - (1 + disp1 * h1my) ** (-1 / disp1) - h1my * ((1 + disp1 * h1my) ** (- 1 - 1 / disp1))), 0)
    sign(y - mu1) * sqrt(-2 * wt * (-(y + 1 / disp1) * log(1 + mu * disp1) + y * log(mu / disp1) - log(1 - (1 + disp1 * mu) ** (-1 / disp1) - mu * ((1 + disp1 * mu) ** (- 1 - 1 / disp1))) +
                                     (y + 1 / disp1) * log(1 + h1my * disp1) - y * logh1mydivdisp1 + logprobdey))
  }

  pointEst <- function (disp, pw, lambda, contr = FALSE) {
    z <- exp(-disp)
    S <- 1 / (1 + lambda / z)
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    N <- (pw * (1 - lambda * (S ** (1 + z))) / prob)
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (beta, pw, lambda, disp, cov, X) {
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

  simulate <- function(n, lambda, theta, lower=1, upper=Inf) {
    lb <- stats::pnbinom(lower, mu=lambda, size = theta)
    ub <- stats::pnbinom(upper, mu=lambda, size = theta)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qnbinom(p_u, mu=lambda, size = theta)
    sims
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
      simulate=simulate,
      family = "zotnegbin"
    ),
    class = "family"
  )
}

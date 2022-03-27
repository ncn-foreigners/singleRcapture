#' Zero One Truncated Negative Binomial model
#'
#' @return A object of class "family" containing objects \cr
#' make_minusloglike(y,X) - for creating negative likelihood function \cr
#' make_gradient(y,X) - for creating gradient function \cr
#' make_hessian(X) - for creating hessian \cr
#' linkfun - a link function to connect between linear predictor and model parameter in regression and a name of link function\cr
#' linkinv - an inverse function of link \cr
#' Dlink - a 1st derivative of link function \cr
#' mu.eta,Variance - Expected Value and Variance \cr
#' aic - for aic computation\cr
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

  mu.eta <- function(eta, disp) {
    A <- exp(disp)
    lambda <- invlink(eta)
    Pr <- (1 + A * lambda) ** (-1 / A)
    Pr <- Pr + lambda * ((1 + A * lambda) ** (-1 - 1 / A))
    G <- (lambda - lambda * ((1 + A * lambda) ** (-1 - 1 / A))) / (1 - Pr)
    G
  }

  variance <- function(mu, disp) {
    A <- exp(disp)
    mu * (1 + A * mu)
  }

  # These three functions are used only for the purposes of computation
  compgamma <- function(y, alpha) {
    temp <- 0:(y - 1)
    sum(log(temp + 1 / alpha))
  }
  compdigamma <- function(y, alpha) {
    temp <- 0:(y - 1)
    sum(-(alpha ** (-2)) / (temp + 1 / alpha))
  }

  compsecond <- function(y, alpha) {
    temp <- 0:(y - 1)
    sum(temp / (((1 + temp * alpha)) ** 2))
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

      -sum(weight * (sapply(y, FUN = {function (y) compgamma(y, alpha = alpha)})
           - log(factorial(y)) - (y + z) * log(1 + lambda / z)
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
                y + alpha * sapply(y,
                FUN = {function (y) compdigamma(y, alpha = alpha)})) * weight)
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

      G00 <- sum(weight * alpha * sapply(y,
             FUN = {function (y) compsecond(y, alpha = alpha)}) +
             z * (((M ** cp26) * cp5 * (cp5 * (2 - cp26) + cp15 * cp3)) +
             cp22 * (cp5 * cp25 + cp15 * cp7) -
             2 * cp21 + cp9 * (M ** (cp26 - 1)) +
             cp22 - cp23) / (M * cp18) -
             z * cp24 * (cp8 * (cp25 * cp5 + cp15 * cp7) - cp9) /
             (M * (cp18 ** 2)) - z * (cp5 * cp10 - cp21 + cp22 - cp23) /
             (M * cp18) - lambda * cp24 / (cp16 * cp18))

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
    NULL
  }

  aic <- function(y, mu, wt, dev) {
    S <- 1 / (1 + mu / dev)
    prob <- S ** dev
    prob <- prob + mu * (S ** (1 + dev))
    prob <- 1 - prob
    -2 * sum(wt * (sapply(y,
    FUN = {function (y) compgamma(y, alpha = 1 / dev)})
    -log(factorial(y)) - (y + dev) * log(1 + mu / dev) +
    y * log(mu / dev) - log(prob)))
  }

  pointEst <- function (disp, pw, lambda) {
    z <- exp(-disp)
    S <- 1 / (1 + lambda / z)
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    N <- sum(pw * (1 - lambda * (S ** (1 + z))) / prob)
    N
  }

  popVar <- function (beta, pw, lambda, disp, hess, X) {
    alpha <- exp(disp)
    z <- 1 / alpha
    M <- 1 + lambda / z
    S <- 1 / M
    prob <- 1 - S ** z - lambda * (S ** (1 + z))
    I <- as.matrix(-hess(beta))
    
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
    
    f1 <-  t(bigTheta) %*% solve(I) %*% bigTheta
    f2 <-  sum(pw * ((1 - lambda * (S ** (1 + z))) ** 2) *
                 (1 - prob) / (prob ** 2))
    
    variation <- f1 + f2
    variation
  }

  R <- list(make_minusloglike = minusLogLike,
            make_gradient = gradient,
            make_hessian = hessian,
            linkfun = link,
            linkinv = invlink,
            dlink = dlink,
            mu.eta = mu.eta,
            aic = aic,
            link = "log",
            valideta = function (eta) {TRUE},
            variance = variance,
            dev.resids = dev.resids,
            validmu = validmu,
            pointEst = pointEst,
            popVar= popVar,
            family = "zotnegbin")
  class(R) <- "family"
  R
}

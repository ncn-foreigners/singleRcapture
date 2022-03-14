#' Zero Truncated Negative Binomial model
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
ztnegbin <- function() {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }

  mu.eta <- function(eta, disp) {
    A <- exp(disp)
    lambda <- invlink(eta)
    pr <- (1 + A * lambda) ** (-1 / A)
    G <- lambda / (1 - pr)
    G
  }

  variance <- function(mu, disp) {
    A <- exp(disp)
    mu * (mu * A + 1)
  }

  # These three functions are used only for the purposes of computation
  compgamma <- function(y, alpha) {
    temp <- 0:(y-1)
    sum(log(temp + 1 / alpha))
  }
  compdigamma <- function(y, alpha) {
    temp <- 0:(y-1)
    sum(-(alpha ** (-2)) / (temp + 1 / alpha))
  }

  comptrigamma <- function(y, alpha) {
    temp <- 0:(y-1)
    sum(- (temp ** 2) / (((1 + temp * alpha)) ** 2))
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
      eta <- as.matrix(X) %*% beta
      lambda <- exp(eta)
      -sum(weight * (sapply(y, FUN = {function (y) compgamma(y, alpha = alpha)})
                    - log(factorial(y)) - (y + z) * log(1 + lambda / z) +
                    y * log(lambda / z) - log(1 - ((1 + lambda / z) ** (-z)))))
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
      S <- 1 / (1 + lambda / z)
      G <- 1 / (1 - (1 + lambda / z) ** (-z))

      # alpha derivative

      G0 <- sum(weight * ((z ** 2) * log(1 / S) +
                sapply(y, FUN = {function (y) compdigamma(y, alpha = alpha)}) +
                y * z - (y + z) * lambda * S) +
                G * (S ** z) *
                (log(1 / S) * (z ** 2) - lambda * z * S))
      # Correction for taking the derivative with respect to log(alpha)
      G0 <- G0 * alpha

      # Beta derivative
      G1 <- t(((y + (lambda - y) * (S ** (-z))) * S /
                 (1 - (1 / S) ** z))  * weight) %*% X

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
      S <- 1 / (1 + lambda / z)
      G <- 1 / (1 - (1 + lambda / z) ** (-z))
      res <- matrix(1, nrow = length(arg), ncol = length(arg))

      # alpha derivative
      G0 <- sum(weight * ((z ** 2) * log(1 / S) +
                            sapply(y, FUN = {function (y) compdigamma(y, alpha = alpha)}) +
                            y * z - (y + z) * lambda * S) +
                  G * (S ** z) *
                  (log(1 / S) * (z ** 2) - lambda * z * S))
      # 2nd alpha derivative
      M <- ((1 + lambda / z) ** z) - 1

      G00 <- sum(2 * S * lambda * (z ** 2) + 2 * log(S) * (z ** 3) +
                   sapply(y, FUN = {function (y) comptrigamma(y, alpha = alpha)}) +
                   (y + z) * (lambda ** 2) * (S ** 2) +
                   (z ** 3) * 2 * (lambda / z + log(S) / S) * (S / M) +
                   (z ** 2) * (S ** (1 - z)) * (lambda * z * S + log(S) * (z ** 2)) *
                   (lambda / z + log(S) / S) / (M ** 2) +
                   (z ** 2) * lambda * log(1 / S) * (S / M) +
                   (z ** 2) * lambda * (S ** 2) * (lambda / z + log(S) / S) / M)

      # Correction for taking the derivative with respect to log(alpha)
      G00 <- G00 * (alpha ** 2) + G0 * alpha
      # mixed derivative
      T1 <- lambda * (y - lambda) * (S ** 2)

      T2 <- ((1 / S) ** (-1 - z))
      T2 <- T2 * (lambda * (1 + z) * S + log(S) * (z ** 2))
      T2 <- G * T2

      T3 <- (1 / S) ** (-z)
      T3 <- T3 * (-log(S) * (z ** 2) - z * lambda * S)
      T3 <- (G ** 2) * (-T3) / ((1 / S) ** (1 + z))

      G01 <- t(X) %*% as.numeric(-T1 + lambda * (T2 + T3))

      G01 <- G01 * alpha
      # second beta derivative
      C1 <- as.numeric((((1 / S) ** z) * (lambda - 1) + 1) * (S ** 2) /
                         (((1 / S) ** z - 1) ** 2))
      C2 <- (1 + y / z) * (S ** 2)

      G11 <- t(as.data.frame(X) * lambda * (C1 - C2)) %*% X

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
    -2 * sum((log(gamma(y + 1 / dev)) - log(gamma(1 / dev)) -
            log(factorial(y)) - (y + 1 / dev) * log(1 + dev * mu) +
            y * log(dev * mu) - log(1 - (1 + dev * mu) ** (-1 / dev))) * wt)
  }

  pointEst <- function (disp, pw, lambda) {
    disp <- exp(disp)
    pr <- 1 - (1 + disp * lambda) ** (- 1 / disp)
    N <- sum(pw / pr)
    N
  }

  popVar <- function (beta, pw, lambda, disp, hess, X) {
    z <- exp(disp)
    pr <- 1 - (1 + z * lambda) ** (- 1 / z)
    N <- sum(1 / pr)
    S <- 1 / (1 + z * lambda)
    I <- as.matrix(-hess(beta))

    bigTheta1 <- sum(pw * ((S ** (1 - 1 / z)) *
                    ((1 / S) * log(1 / S) - z * lambda) /
                    ((z ** 2) * ((1 - (1 / S) ** (1 / z)) ** 2))))

    bigTheta2 <- -(pw * as.numeric(lambda * (S ** (1 - 1 / z)) /
                  ((1 - (1 / S) ** (1 / z)) ** 2))) %*% as.matrix(X)
    # Correction for taking the derivative with respect to log(alpha)
    # control = list(
    # reltol = .Machine$double.eps)
    bigTheta2 <- bigTheta2 * z

    bigTheta <- matrix(c(bigTheta1, bigTheta2), ncol = 1)

    f1 <-  t(bigTheta) %*% solve(I) %*% bigTheta
    f2 <-  sum(pw * (1 - pr) / (pr ** 2))

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
            family = "ztnegbin")
  class(R) <- "family"
  R
}

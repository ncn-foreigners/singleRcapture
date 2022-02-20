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
  Invlink <- exp
  Dlink <- function(Lambda) {
    1 / Lambda
  }

  mu.eta <- function(eta, disp) {
    A <- disp
    Lambda <- Invlink(eta)
    Pr <- (1 + A * Lambda) ** (-1 / A)
    G <- Lambda / (1 - Pr)
    G
  }

  Variance <- function(mu, disp) {
    A <- disp
    Lambda <- mu
    Pr <- (1 + A * Lambda) ** (-1 / A)
    G <- ((Lambda + Lambda * (A ** 2)) / (1 - Pr) -
    Pr * ((Lambda / (1 - Pr)) ** 2))
    G
  }

  # These three functions are used only for the purposes of computation
  Compgamma <- function(y, alpha) {
    Temp <- 0:(y-1)
    sum(log(Temp + 1 / alpha))
  }
  Compdigamma <- function(y, alpha) {
    Temp <- 0:(y-1)
    sum(-(alpha ** (-2)) / (Temp + 1 / alpha))
  }

  Comptrigamma <- function(y, alpha) {
    Temp <- 0:(y-1)
    sum(- (Temp ** 2) / (((1 + Temp * alpha)) ** 2))
  }

  MinusLogLike <- function(y, X, weight = 1) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(ARG) {
      alpha <- ARG[1]
      z <- 1 / alpha
      beta <- ARG[-1]
      Eta <- as.matrix(X) %*% beta
      Lambda <- exp(Eta)
      -sum(weight * (sapply(y, FUN = {function (y) Compgamma(y, alpha = alpha)})
        - log(factorial(y)) - (y + z) * log(1 + Lambda / z) +
        y * log(Lambda / z) - log(1 - ((1 + Lambda / z) ** (-z)))))
    }
  }


  Gradient <- function(y, X, weight = 1) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(ARG) {
      alpha <- ARG[1]
      # z is inverse of alpha in documentation
      z <- 1 / alpha
      beta <- ARG[-1]
      Eta <- X %*% beta
      Lambda <- exp(Eta)
      S <- 1 / (1 + Lambda / z)
      G <- 1 / (1 - (1 + Lambda / z) ** (-z))

      # alpha derivative

      G0 <- sum(weight * ((z ** 2) * log(1 / S) +
                sapply(y, FUN = {function (y) Compdigamma(y, alpha = alpha)}) +
                y * z - (y + z) * Lambda * S) +
                G * (S ** z) *
                (log(1 / S) * (z ** 2) - Lambda * z * S))


      # Beta derivative
      G1 <- t(((y + (Lambda - y) * (S ** (-z))) * S /
               (1 - (1 / S) ** z))  * weight) %*% X

      c(G0, G1)
    }
  }

  Hessian <- function(y, X, weight = 1) {
    if (is.null(weight)) {
      weight <- 1
    }
    y <- as.numeric(y)
    X <- as.matrix(X)

    function(ARG) {
      alpha <- ARG[1]
      # z is inverse of alpha in documentation
      z <- 1 / alpha
      beta <- ARG[-1]
      Eta <- X %*% beta
      Lambda <- exp(Eta)
      S <- 1 / (1 + Lambda / z)
      G <- 1 / (1 - (1 + Lambda / z) ** (-z))
      res <- matrix(1, nrow = length(ARG), ncol = length(ARG))

      # 2nd alpha derivative
      M <- ((1 + Lambda / z) ** z) - 1

      G00 <- sum(2 * S * Lambda * (z ** 2) + 2 * log(S) * (z ** 3) +
      sapply(y, FUN = {function (y) Comptrigamma(y, alpha = alpha)}) +
      (y + z) * (Lambda ** 2) * (S ** 2) +
      (z ** 3) * 2 * (Lambda / z + log(S) / S) * (S / M) +
      (z ** 2) * (S ** (1 - z)) * (Lambda * z * S + log(S) * (z ** 2)) *
      (Lambda / z + log(S) / S) / (M ** 2) +
      (z ** 2) * Lambda * log(1 / S) * (S / M) +
      (z ** 2) * Lambda * (S ** 2) * (Lambda / z + log(S) / S) / M)

      # mixed derivative
      T1 <- Lambda * (y - Lambda) * (S ** 2)

      T2 <- ((1 / S) ** (-1 - z))
      T2 <- T2 * (Lambda * (1 + z) * S + log(S) * (z ** 2))
      T2 <- G * T2

      T3 <- (1 / S) ** (-z)
      T3 <- T3 * (-log(S) * (z ** 2) - z * Lambda * S)
      T3 <- (G ** 2) * (-T3) / ((1 / S) ** (1 + z))

      G01 <- t(X) %*% as.numeric(-T1 + Lambda * (T2 + T3))

      # second beta derivative
      C1 <- as.numeric((((1 / S) ** z) * (Lambda - 1) + 1) * (S ** 2) /
                         (((1 / S) ** z - 1) ** 2))
      C2 <- (1 + y / z) * (S ** 2)

      G11 <- t(as.data.frame(X) * Lambda * (C1 - C2)) %*% X

      res[1, 1] <- G00
      res[-1, -1] <- G11
      res[1, -1] <- G01
      res[-1, 1] <- G01

      res
    }
  }

  validmu <- function(mu) {
    is.finite(mu) && all(mu > 0)
  }

  dev.resids <- function (y, mu, wt, disp = NULL) {
    NULL
  }

  aic <- function(y, mu, wt, dev) {
    -2 * sum((log(gamma(y + 1 / dev)) - log(gamma(1 / dev)) -
    log(factorial(y)) - (y + 1 / dev) * log(1 + dev * mu) +
    y * log(dev * mu) - log(1 - (1 + dev * mu) ** (-1 / dev))) * wt)
  }

  Point.est <- function (disp, pw, Lambda) {
    Pr <- 1 - (1 + disp * Lambda) ** (- 1 / disp)
    N <- sum(pw / Pr)
    N
  }

  Pop.var <- function (beta, pw, Lambda, disp, Hess, X) {
    z <- disp
    Pr <- 1 - (1 + z * Lambda) ** (- 1 / z)
    N <- sum(1 / Pr)
    S <- 1 / (1 + z * Lambda)
    Inform <- as.matrix(-Hess(beta))

    BigTheta1 <- sum(pw * ((S ** (1 - 1 / z)) *
                     ((1 / S) * log(1 / S) - z * Lambda) /
                     ((z ** 2) * ((1 - (1 / S) ** (1 / z)) ** 2))))

    BigTheta2 <- -(pw * as.numeric(Lambda * (S ** (1 - 1 / z)) /
                  ((1 - (1 / S) ** (1 / z)) ** 2))) %*% as.matrix(X)

    BigTheta <- matrix(c(BigTheta1, BigTheta2), ncol = 1)

    f1 <-  t(BigTheta) %*% solve(Inform) %*% BigTheta
    f2 <-  sum(pw * (1 - Pr) / (Pr ** 2))

    Variation <- f1 + f2
    Variation
  }

  R <- list(make_minusloglike = MinusLogLike,
            make_gradient = Gradient,
            make_hessian = Hessian,
            linkfun = link,
            linkinv = Invlink,
            Dlink = Dlink,
            mu.eta = mu.eta,
            aic = aic,
            link = "log",
            valideta = function (eta) {TRUE},
            variance = Variance,
            dev.resids = dev.resids,
            validmu = validmu,
            Point.est = Point.est,
            Pop.var= Pop.var,
            family = "ztnegbin")
  class(R) <- "family"
  R
}

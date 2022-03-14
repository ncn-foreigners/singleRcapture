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

      # alpha derivative
      G0 <- sum((z * log(M) - alpha * lambda * (y + z) * S -
            ((S ** z) * (z * log(S) + lambda * S) - lambda *
            (S ** (1 + z)) * (z * log(M) + lambda * alpha *
            (-1 - z) * S)) / prob + y + alpha * sapply(y,
            FUN = {function (y) compdigamma(y, alpha = alpha)})) * weight)
      # Beta derivative
      G1 <- t(X) %*% (weight * (y - alpha * S * (y + z) * lambda +
            lambda * alpha * (-1 - z) * lambda * (S ** (2 + z)) / prob))

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

      # 2nd alpha derivative

      G00 <- sum(alpha * sapply(y,
             FUN = {function (y) compsecond(y, alpha = alpha)}) +
             z * (((M ** (z * (1 + 2 * alpha))) * log(M) * (
             log(M) * (2 - z * (2 * alpha + 1)) +
             lambda * (2 * alpha + 1) * S)) + (y - lambda) * alpha *
             (M ** (1 + z)) * (log(M) * (1 - z * (1 + alpha)) +
             lambda * S * (1 + alpha)) - 2 * lambda * (y + lambda) *
             (alpha ** 2) + lambda * alpha * (M ** (z * (2 * alpha + 1) - 1)) +
             (y - lambda) * alpha * (M ** (1 + z)) - (1 + lambda) * y * alpha) /
             (M * (M ** (1 + z) - alpha * lambda - lambda - 1)) -
             z * (log(M) * (M ** (2 + z)) - lambda * (y + lambda) * (alpha ** 2) +
             (y - lambda) * alpha * (M ** (1 + z)) - y * alpha * (1 + lambda)) *
             ((M ** (1 + z)) * ((1 - z * (1 + alpha)) * log(M) + lambda *
             S * (1 + alpha)) - lambda * alpha) /
             (M * ((M ** (1 + z) - alpha * lambda - lambda - 1) ** 2)) -
             z * (log(M) * (M ** (2 + z)) - lambda * (y + lambda) *
             (alpha ** 2) + (y - lambda) * alpha * (M ** (1 + z)) -
             (1 + lambda) * y * alpha) /
             (M * (M ** (1 + z) - alpha * lambda - lambda - 1)) -
             lambda * (log(M) * (M ** (2 + z)) - lambda * (y + lambda) *
             (alpha ** 2) + (y - lambda) * alpha * (M ** (1 + z)) -
             (1 + lambda) * y * alpha) /
             ((M ** 2) * (M ** (1 + z) - alpha * lambda - lambda - 1)))

      # mixed derivative
      term <- (-alpha * (y + z) * S +
               lambda * (alpha ** 2) * (y + z) * (S ** 2) +
               ((-1 - z) * alpha * lambda * (S ** (z + 2))) / prob +
               (lambda * (S ** (2 + z))) / prob + S +
               alpha * lambda * (-1 - z) * (S ** (2 + z)) *
               (alpha * lambda * (-2 - z) * S + log(M) * z) / prob -
               (alpha * (-1 - z) * lambda * (S ** (2 + z)) *
               ((S ** z) * (z * log(S) + lambda * S) -
               lambda * (S ** (1 + z)) * (z * log(M) +
               S * lambda * alpha * (-1 - z)))) / (prob ** 2))

      G01 <- t(X) %*% as.numeric(term * lambda)
      # second beta derivative
      term <- ((-(alpha ** 3) * y * (lambda ** 2) *
              (((M ** z) - 1) ** 2) + (alpha ** 2) * lambda *
              ((lambda ** 2) * (M  ** z) - lambda * ((M ** z) - 1) *
              (1 - 2 * y + (M ** z)) - 2 * y * (((M ** z) - 1) ** 2)) +
              (lambda ** 2) * (M ** z) + alpha * ((lambda ** 3) * (M ** z) +
              (lambda ** 2) * (1 - y + (M ** z)) -2 * lambda * (M ** z - 1) *
              (M ** z - y) - y * ((M ** z - 1) ** 2)) - (M ** z - 1) ** 2) /
              ((M ** 2) * ((M ** z +
              lambda * (alpha * (M ** z - 1) - 1) - 1) ** 2)))

      G11 <- t(as.data.frame(X) * lambda * as.numeric(term)) %*% X


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

    bigTheta1 <- sum(pw * as.numeric((prob * (-lambda * (S ** (1 + z)) *
                     (z * log(M) + lambda * alpha * (-1 - z) * S)) -
                     ((1 - lambda * (S ** (1 + z))) *
                     (-lambda * (S ** (1 + z)) * (z * log(M) +
                     alpha * lambda * S * (-1 - z)) -
                     (z * log(M) - lambda * S) / (M ** z)))) / (prob ** 2)))

    bigTheta2 <- t(as.matrix(X)) %*% (pw * as.numeric(lambda *
                  (prob * ((lambda - 1) * (S ** (2 + z))) -
                  (1 - lambda * (S ** (1 + z))) * (lambda *
                  (1 + alpha) * (S ** (2 + z)))) / (prob ** 2)))

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

#' Zero truncated Negative Binomial model
#'
#' @return A object of class "family" containing objects \cr
#' make_minusloglike(y,X) - for creating negative likelihood function \cr
#' make_gradient(y,X) - for creating gradient function \cr
#' make_hessian(X) - for creating hessian \cr
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
#' @export
ztnegbin <- function() {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }

  mu.eta <- function(eta, disp, type = "trunc") {
    A <- exp(disp)
    lambda <- invlink(eta)
    switch (type,
      nontrunc = lambda,
      trunc = lambda / (1 - (1 + A * lambda) ** (-1 / A))
    )
  }

  variance <- function(mu, disp, type = "nontrunc") {
    A <- exp(disp)
    switch (type,
      nontrunc = mu * (1 + A * mu),
      trunc = (mu + A * (mu ** 2) - A * (mu ** 2) * ((1 + A * mu) ** (-1 / A))) / ((1 - (1 + A * mu) ** (-1 / A)) ** 2)
    )
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
      M <- 1 + lambda / z

      -sum(weight * (lgamma(y + z) - lgamma(z) -
      log(factorial(y)) - (y + z) * log(M) +
      y * log(lambda / z) - log(1 - (M ** (-z)))))
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
      cp1 <- log(1 / S) * (z ** 2)
      cp2 <- lambda * S
      cp3 <- S ** (-z)

      # log(alpha) derivative
      #G0 <- sum(weight * (cp1 +
      #          sapply(y, FUN = {function (y) compdigamma(y, alpha = alpha)}) +
      #          y * z - (y + z) * cp2 +
      #          G * (1 / cp3) * (cp1 - z * cp2)) * alpha)
      G0 <- sum(weight * (cp1 -
                            (digamma(y + z) - digamma(z)) * (z ** 2) +
                            y * z - (y + z) * cp2 +
                            G * (1 / cp3) * (cp1 - z * cp2)) * alpha)

      # Beta derivative
      G1 <- t(((y + (lambda - y) * cp3) * S / (1 - cp3))  * weight) %*% X

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
      M <- ((1 + lambda / z) ** z) - 1
      S <- 1 / (1 + lambda / z)
      G <- 1 / (1 - (1 + lambda / z) ** (-z))
      res <- matrix(1, nrow = length(arg), ncol = length(arg))
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

      # log(alpha) derivative
      #G0 <- sum(weight * (cp1 +
      #          sapply(y, FUN = {function (y) compdigamma(y, alpha = alpha)}) +
      #          y * z - (y + z) * cp2 +
      #          G * (1 / cp3) * (cp1 - z * cp2)) * alpha)
      G0 <- sum(weight * (cp1 -
                            (digamma(y + z) - digamma(z)) * (z ** 2) +
                            y * z - (y + z) * cp2 +
                            G * (1 / cp3) * (cp1 - z * cp2)) * alpha)

      # 2nd log(alpha) derivative
      #G00 <- sum((2 * cp2 * cp4 + 2 * cp7 * (z ** 3) +
      #            sapply(y, FUN = {function (y) comptrigamma(y, alpha = alpha)}) +
      #            (y + z) * (lambda ** 2) * cp6 +
      #            (z ** 3) * 2 * cp8 * cp11 +
      #            cp4 * (S ** (1 - z)) * (z * cp2 + cp7 * cp4) *
      #            cp8 / (M ** 2) +
      #            cp4 * lambda * log(cp9) * cp11 +
      #            cp4 * lambda * cp6 * cp8 / M) * (alpha ** 2) * weight)
      G00 <- sum((2 * cp2 * cp4 + 2 * cp7 * (z ** 3) +
                    (trigamma(y + z) - trigamma(z)) * (z ** 4) -
                    (digamma(y + z) - digamma(z)) * (z ** 2) +
                    (y + z) * (lambda ** 2) * cp6 +
                    (z ** 3) * 2 * cp8 * cp11 +
                    cp4 * (S ** (1 - z)) * (z * cp2 + cp7 * cp4) *
                    cp8 / (M ** 2) +
                    cp4 * lambda * log(cp9) * cp11 +
                    cp4 * lambda * cp6 * cp8 / M) * (alpha ** 2) * weight)

      # Correction for taking the derivative with respect to log(alpha)
      G00 <- G00 + G0

      # mixed derivative
      T1 <- lambda * (y - lambda) * cp6

      T2 <- cp10 * (lambda * (1 + z) * S + cp7 * cp4)
      T2 <- G * T2

      T3 <- cp9 ** (-z)
      T3 <- T3 * (-cp7 * cp4 - z * cp2)
      T3 <- (G ** 2) * (-T3) * cp10

      G01 <- t(X) %*% as.numeric(-T1 + lambda * (T2 + T3)) * alpha * weight

      # second beta derivative
      C1 <- as.numeric(((cp9 ** z) * (lambda - 1) + 1) *
                       cp6 / ((cp9 ** z - 1) ** 2))
      C2 <- (1 + y / z) * cp6

      G11 <- t(as.data.frame(X) * lambda * (C1 - C2) * weight) %*% X

      res[1, 1] <- G00
      res[-1, -1] <- G11
      res[1, -1] <- G01
      res[-1, 1] <- G01

      res
    }
  }

  validmu <- function(mu) {
    all(is.finite(mu)) && all(0 < mu)
  }

  dev.resids <- function (y, mu, wt, disp = NULL) {
    eta <- log(mu)
    disp1 <- exp(disp)
    mu1 <- mu.eta(eta = eta, disp = disp)
    a <- function(y) {stats::uniroot(f = function(x) {mu.eta(x, disp = disp) - y}, lower = -log(y), upper = y * 10, tol = .Machine$double.eps)$root}
    hm1y <- y
    hm1y[y == 1] <- -16 # This theoretically outght to be -Inf but that would cause an error this is basically the lowest value that wont cause any
                           # errors. It approximates mu.eta(eta = hm1y) for its values so well it probably wont cause any incorrect decisions.
                           # Manually setting the values of to-be logarithms to zero for values with y = 1 also causes errors.
                           # This is an imperfect solution but it will do for now.
    hm1y[y > 1] <- sapply(y[y > 1], FUN = a)
    loghm1ytdisp <- log(disp1 * exp(hm1y))
    logprobhm1y <- log(1 - ((1 + disp1 * exp(hm1y)) ** (-1 / disp1)))
    sign(y - mu1) * sqrt(-2 * wt * (-(y + 1 / disp1) * log(1 + mu * disp1) + y * log(mu * disp1) - log(1 - ((1 + mu * disp1) ** (-1/disp1))) +
                                    (y + 1 / disp1) * log(1 + exp(hm1y) * disp1) - y * loghm1ytdisp + logprobhm1y))
  }

  pointEst <- function (disp, pw, lambda, contr = FALSE) {
    disp <- exp(disp)
    pr <- 1 - (1 + disp * lambda) ** (- 1 / disp)
    N <- pw / pr
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (beta, pw, lambda, disp, hess, X) {
    z <- exp(disp)
    pr <- 1 - (1 + z * lambda) ** (- 1 / z)
    S <- 1 / (1 + z * lambda)
    I <- as.matrix(-hess(beta))

    cp1 <- (1 / S)
    cp2 <- (S ** (1 - 1 / z))
    cp3 <- ((1 - cp1 ** (1 / z)) ** 2)

    bigTheta1 <- sum(pw * (cp2 * (cp1 * log(cp1) - z * lambda) /
                          ((z ** 2) * cp3)))

    bigTheta2 <- -(pw * as.numeric(lambda * cp2 / cp3)) %*% as.matrix(X)
    # Correction for taking the derivative with respect to log(alpha)
    bigTheta1 <- bigTheta1 * z


    bigTheta <- matrix(c(bigTheta1, bigTheta2), ncol = 1)

    f1 <-  t(bigTheta) %*% solve(I) %*% bigTheta
    f2 <-  sum(pw * (1 - pr) / (pr ** 2))

    f1 + f2
  }

  structure(
    list(
      make_minusloglike = minusLogLike,
      make_gradient = gradient,
      make_hessian = hessian,
      linkfun = link,
      linkinv = invlink,
      dlink = dlink,
      mu.eta = mu.eta,
      link = "log",
      valideta = function (eta) {TRUE},
      variance = variance,
      dev.resids = dev.resids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      family = "ztnegbin"
    ),
    class = "family"
  )
}

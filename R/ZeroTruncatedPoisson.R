#' Zero truncated Poisson Model
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
#' @importFrom lamW lambertW0
#' @export
ztpoisson <- function() {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }

  mu.eta <- function(disp = NULL, eta, type = "trunc") {
    lambda <- invlink(eta)
    switch (type,
      "nontrunc" = lambda,
      "trunc" = lambda / (1 - exp(-lambda))
    )
  }

  variance <- function(disp = NULL, mu, type = "nontrunc") {
    switch (type,
      "nontrunc" = mu,
      "trunc" = mu.eta(eta = log(mu)) * (1 + mu - mu.eta(eta = log(mu)))
    )
  }

  minusLogLike <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    function(beta) {
      eta <- as.matrix(X) %*% beta
      lambda <- exp(eta)
      -sum(weight * (y * eta - log(exp(lambda) - 1) - log(factorial(y))))
    }
  }

  gradient <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      lambda <- exp(as.matrix(X) %*% beta)
      mu <- lambda / (1 - exp(-lambda))
      t(as.matrix(X)) %*% (weight * (y - mu))
    }
  }

  hessian <- function(y, X, weight = 1) {
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      lambda <- exp(as.matrix(X) %*% beta)
      eml <- exp(-lambda)
      coefficient <- (1 / (1 - eml) - lambda * eml / ((1 - eml) ** 2))

      dmu <- diag(weight * as.numeric(coefficient))
      dlam <- as.matrix(X * as.numeric(lambda))

      -((t(as.matrix(X)) %*% dmu) %*% dlam)
    }
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  dev.resids <- function(y, mu, wt, disp = NULL) {
    eta <- log(mu)
    mu1 <- mu.eta(eta = eta)
    #hm1y <- ifelse(y > 1, VGAM::lambertW(-y * exp(-y)) + y, 0)
    hm1y <- ifelse(y > 1, lamW::lambertW0(-y * exp(-y)) + y, 0)
    log1mexphm1y <- ifelse(y > 1, log(1 - exp(-hm1y)), 0)
    loghm1y <- ifelse(y > 1, log(hm1y), 0)
    #loghm1y <- ifelse(hm1y > )
    sign(y - mu1) * sqrt(-2 * wt * (y * eta - mu - log(1 - exp(-mu)) - y * loghm1y + hm1y + log1mexphm1y))
  }

  pointEst <- function (disp = NULL, pw, lambda, contr = FALSE) {
    N <- pw / (1 - exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (beta, pw, lambda, disp = NULL, hess, X) {
    X <- as.data.frame(X)
    ml <- (1 - exp(-lambda)) ** 2
    I <- -hess(beta)


    f1 <- colSums(-X * pw * (exp(log(lambda) - lambda) / ml))
    f1 <- t(f1) %*% solve(as.matrix(I)) %*% f1

    f2 <- sum(pw * exp(-lambda) / ml)

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
      family = "ztpoisson"
    ),
    class = "family"
  )
}

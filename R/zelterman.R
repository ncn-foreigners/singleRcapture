#' Zelterman's model for population estimate
#'
#' @return A object of class "family" containing objects \cr
#' makeMinusLogLike(y,X) - for creating negative likelihood function \cr
#' makeGradient(y,X) - for creating gradient function \cr
#' makeHessian(X) - for creating hessian \cr
#' linkfun - a link function to connect between linear predictor and model parameter in regression and a name of link function\cr
#' linkinv - an inverse function of link \cr
#' dlink - a 1st derivative of link function \cr
#' mu.eta,variance - Expected Value and variance \cr
#' valedmu, valideta - for checking if regression arguments and valid\cr
#' family - family name\cr
#' Where: \cr
#' y is a vector of observed values \cr
#' X is a matrix / data frame of covariates
#' @export
zelterman <- function() {
  link <- function(x) {log(x / 2)}
  invlink <- function (x) {2 * exp(x)}
  dlink <- function(lambda) {1 / lambda}
  
  mu.eta <- function(eta, ...) {
    1 / (1 + exp(-eta))
  }
  
  variance <- function(eta, ...) {
    mu <- mu.eta(eta)
    mu * (1 - mu)
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- invlink(eta)
    L1 <- lambda / 2
    (L1 / ((1 + L1) ** 2))
  }
  
  funcZ <- function(eta, weight, y, mu, ...) {
    lambda <- invlink(eta)
    L1 <- lambda / 2
    (L1 * (y - 1) + y) / (L1 + 1) / weight
  }

  minusLogLike <- function(y, X, weight = 1, ...) {
    y <- as.numeric(y)
    z <- y
    z[z == 1] <- 0
    z[z == 2] <- 1
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      eta <- as.matrix(X) %*% beta
      lambda <- invlink(eta)
      L1 <- lambda / 2
      par <- L1 / (1 + L1)
      -sum(weight * (z * log(par) + (1 - z) * log(1 - par)))
    }
  }

  gradient <- function(y, X, weight = 1, ...) {
    y <- as.numeric(y)
    z <- y - 1
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      eta <- as.matrix(X) %*% beta
      lambda <- invlink(eta)
      L1 <- lambda / 2
      t(X) %*% (weight * (L1 * (z - 1) + z) / (L1 + 1))
    }
  }

  hessian <- function(y, X, weight = 1, ...) {
    y <- as.numeric(y)
    z <- y
    z[z == 1] <- 0
    z[z == 2] <- 1
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      eta <- as.matrix(X) %*% beta
      lambda <- invlink(eta)
      L1 <- lambda / 2
      term <- -(L1 / ((1 + L1) ** 2))
      t(as.data.frame(X) * weight * term) %*% as.matrix(X)
    }
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(1 > mu)
  }

  dev.resids <- function(y, eta, wt, ...) {
    mu <- invlink(eta)
    z <- y - 1
    mu1 <- mu.eta(eta = eta)
    sign(z - mu1) * sqrt(-2 * wt * (z * log(mu1) + (1 - z) * log(1 - mu1)))
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- invlink(eta)
    N <- (pw * (1 / (1 - exp(-lambda))))
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- invlink(eta)
    Xvlm <- as.data.frame(Xvlm)
    prob <- 1 - exp(-lambda)

    f1 <- colSums(-Xvlm * pw * (exp(-lambda) * lambda / (prob ** 2)))
    f1 <- t(f1) %*% as.matrix(cov) %*% f1

    f2 <- sum(pw * (1 - prob) / (prob ** 2))

    f1 + f2
  }
  
  dFun <- function (x, eta, type = "trunc") {
    lambda <- invlink(eta)
    stats::dpois(x = x, lambda = lambda) / (1 - stats::dpois(x = 0, lambda = lambda))
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
      link = "2 * log",
      valideta = function (eta) {TRUE},
      variance = variance,
      Wfun = Wfun,
      funcZ = funcZ,
      dev.resids = dev.resids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      family = "zelterman",
      parNum = 1,
      etaNames = "lambda",
      densityFunction = dFun
    ),
    class = "family"
  )
}

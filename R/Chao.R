#' Chao model for population estimate
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
chao <- function() {
  link <- function(x) log(x / 2)
  Invlink <- function (x) 2 * exp(x)
  Dlink <- function(Lambda) {
    1 / Lambda
  }

  mu.eta <- function(disp = NULL, eta) {
    Lambda <- Invlink(eta)
    (Lambda / 2) / (1 + Lambda / 2)
  }

  Variance <- function(disp = NULL, mu) {
    ((mu / 2) / (1 + mu / 2)) * (1 / (1 + mu / 2))
  }

  MinusLogLike <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    z <- y
    z[z == 1] <- 0
    z[z == 2] <- 1
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      Eta <- as.matrix(X) %*% beta
      Lambda <- Invlink(Eta)
      L1 <- Lambda / 2
      -sum(weight * (z * log(L1 / (1 + L1)) + (1 - z) * log(1 / (1 + L1))))
    }
  }

  Gradient <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    z <- y
    z[z == 1] <- 0
    z[z == 2] <- 1
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      Eta <- as.matrix(X) %*% beta
      Lambda <- Invlink(Eta)
      L1 <- Lambda / 2
      t(X) %*% ((z - L1 / (1 + L1)) * weight)
    }
  }

  Hessian <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    z <- y
    z[z == 1] <- 0
    z[z == 2] <- 1
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      Eta <- as.matrix(X) %*% beta
      Lambda <- Invlink(Eta)
      L1 <- Lambda / 2
      term <- -(L1 / ((1 + L1) ** 2))
      t(as.data.frame(X) * weight * term) %*% as.matrix(X)
    }
  }

  validmu <- function(mu) {
    is.finite(mu) && all(1 > mu)
  }

  dev.resids <- function(y, mu, wt, disp = NULL) {
    NULL
  }

  aic <- function(y, mu, wt, dev) {
    z <- y
    z[z == 1] <- 0
    z[z == 2] <- 1
    L1 <- mu / 2
    -2 * -sum((z * log(L1 / (1 + L1)) + (1 - z) * log(1 / (1 + L1))) * wt)
  }

  Point.est <- function (disp = NULL, pw, Lambda) {
    N <- sum((1 + 1 / (Lambda + (Lambda ** 2) / 2)) * pw)
    N
  }

  Pop.var <- function (beta, pw, Lambda, disp = NULL, Hess, X) {
    X <- as.data.frame(X)
    Inform <- -Hess(beta)
    Prob <- exp(-Lambda) + Lambda * exp(-Lambda)

    f1 <- colSums(-X * pw * ((Lambda + Lambda ** 2) /
                       ((Lambda + (Lambda ** 2) / 2) ** 2)))
    f1 <- t(f1) %*% solve(as.matrix(Inform)) %*% f1

    f2 <- sum(pw * (1 - Prob) * ((1 + exp(-Lambda) / Prob) ** 2))

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
            family = "chao")
  class(R) <- "family"
  R
}

#' Zero Truncated Poisson Model Functions
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
Zero_Truncated_Poisson <- function() {
  link <- log
  Invlink <- exp
  Dlink <- function(Lambda) {
    1 / Lambda
  }

  mu.eta <- function(disp = NULL, eta) {
    Lambda <- Invlink(eta)
    Lambda / (1 - exp(-Lambda))
  }

  Variance <- function(disp = NULL, mu) {
    mu
  }

  MinusLogLike <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    function(beta) {
      Eta <- as.matrix(X) %*% beta
      Lambda <- exp(Eta)
      -sum(weight * (y * Eta - log(exp(Lambda) - 1) - log(factorial(y))))
    }
  }

  Gradient <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      Lambda <- exp(as.matrix(X) %*% beta)
      mu <- Lambda / (1 - exp(-Lambda))
      t(as.matrix(X)) %*% (weight * (y - mu))
    }
  }

  Hessian <- function(y, X, weight = 1) {
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      Lambda <- exp(as.matrix(X) %*% beta)
      coefficient <- 1 / (1 - exp(-Lambda)) - Lambda * exp(-Lambda) / ((1 - exp(-Lambda)) ** 2)
      Dmu <- diag(weight * as.numeric(coefficient))
      Dlam <- as.matrix(X * as.numeric(Lambda))

      -((t(as.matrix(X)) %*% Dmu) %*% Dlam)
    }
  }

  validmu <- function(mu) {
    is.finite(mu) && all(mu > 0)
  }

  dev.resids <- function(y, mu, wt, disp = NULL) {
    NULL
  }

  aic <- function(y, mu, wt, dev) {
    -2 * sum((y * log(mu) - log(exp(mu) - 1) - log(factorial(y))) * wt)
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
            family = "ZTP")
  class(R) <- "family"
  R
}

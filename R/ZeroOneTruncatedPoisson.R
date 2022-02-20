#' ZeroOneTruncatedPoisson model
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
zotpoisson <- function() {
  link <- log
  Invlink <- exp
  Dlink <- function(Lambda) {
    1 / Lambda
  }

  mu.eta <- function(disp = NULL, eta) {
    Lambda <- Invlink(eta)
    (Lambda - Lambda * exp(-Lambda)) / (1 - exp(-Lambda) - Lambda * exp(-Lambda))
  }

  Variance <- function(disp = NULL, mu) {
    EX2 <- ((mu + (mu ** 2) - mu * exp(-mu)) /
    (1 - exp(-mu) - mu * exp(-mu)))
    EX <- ((mu - mu * exp(-mu)) /
    (1 - exp(-mu) - mu * exp(-mu)))
    G <- EX2 - (EX ** 2)
    G
  }

  MinusLogLike <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      Eta <- as.matrix(X) %*% beta
      Lambda <- exp(Eta)
      -sum(weight * (y * Eta - Lambda - log(factorial(y)) -
           log(1 - exp(-Lambda) - Lambda * exp(-Lambda))))
    }
  }

  Gradient <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    function(beta) {
      Lambda <- exp(as.matrix(X) %*% beta)
      term <- Lambda / (exp(Lambda) - Lambda - 1)

      t(as.matrix(X)) %*% (weight * (y - Lambda - Lambda * term))
    }
  }

  Hessian <- function(y, X, weight = 1) {
    X <- as.matrix(X)
    y <- as.numeric(y)

    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      Lambda <- exp(as.matrix(X) %*% beta)

      term <- (((2 + Lambda ** 2) * exp(Lambda) - exp(2 * Lambda) - 1) /
                ((exp(Lambda) - Lambda - 1) ** 2))

      hes <- t(X) %*% as.matrix(t(t(as.data.frame(X) * Lambda * term)))
      hes
    }
  }

  validmu <- function(mu) {
    is.finite(mu) && all(mu > 0)
  }

  dev.resids <- function(y, mu, wt, disp = NULL) {
    NULL
  }

  aic <- function(y, mu, wt, dev) {
    -2 * sum((y * log(mu) - mu - log(1 - exp(-mu) - mu * exp(-mu)) -
                log(factorial(y))) * wt)
  }

  Point.est <- function (disp = NULL, pw, Lambda) {
    N <- sum(pw *(1 - Lambda * exp(-Lambda)) /
            (1 - exp(-Lambda) - Lambda * exp(-Lambda)))
    N
  }

  Pop.var <- function (beta, pw, Lambda, disp = NULL, Hess, X) {
    X <- as.data.frame(X)
    Inform <- -Hess(beta)
    Prob <- (1 - exp(-Lambda) - Lambda * exp(-Lambda))
    term <- (1 - Lambda * exp(-Lambda)) ** 2

    f1 <- t(X) %*% (as.numeric(pw * Lambda * (1 - exp(Lambda)) /
          ((1 + Lambda - exp(Lambda)) ** 2)))

    f1 <- t(f1) %*% solve(as.matrix(Inform)) %*% f1

    f2 <- sum(pw * term * (1 - Prob) / (Prob ** 2))

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
            family = "zotpoisson")
  class(R) <- "family"
  R
}

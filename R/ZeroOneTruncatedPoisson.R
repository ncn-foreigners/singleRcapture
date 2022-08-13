#' Zero-one truncated Poisson model
#'
#' @return A object of class "family" containing objects \cr
#' makeMinusLogLike(y,X) - for creating negative likelihood function \cr
#' makeGradient(y,X) - for creating gradient function \cr
#' makeHessian(X) - for creating hessian \cr
#' linkfun - a link function to connect between linear predictor and model parameter in regression and a name of link function\cr
#' linkinv - an inverse function of link \cr
#' Dlink - a 1st derivative of link function \cr
#' mu.eta,Variance - Expected Value and Variance \cr
#' valedmu, valideta - for checking if regression arguments and valid\cr
#' family - family name\cr
#' Where: \cr
#' y is a vector of observed values \cr
#' X is a matrix / data frame of covariates
#' @export
zotpoisson <- function() {
  link <- log
  invlink <- exp
  dlink <- function(lambda) {
    1 / lambda
  }

  mu.eta <- function(disp = NULL, eta, type = "trunc") {
    lambda <- invlink(eta)
    switch (type,
            "nontrunc" = lambda,
            "trunc" = (lambda - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda))
    )
  }

  variance <- function(disp = NULL, mu, type = "nontrunc") {
    switch (type,
            "nontrunc" = mu,
            "trunc" = (mu - mu * exp(-mu)) / (1 - exp(-mu) - mu * exp(-mu))
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- exp(eta)
    -lambda * (((2 + lambda ** 2) * exp(lambda) - exp(2 * lambda) - 1) /
              ((exp(lambda) - lambda - 1) ** 2))
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- exp(eta)
    eta + (y - lambda - lambda * lambda / (exp(lambda) - lambda - 1)) / weight
  }

  minusLogLike <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      eta <- as.matrix(X) %*% beta
      lambda <- exp(eta)
      -sum(weight * (y * eta - lambda - log(factorial(y)) -
           log(1 - exp(-lambda) - lambda * exp(-lambda))))
    }
  }

  gradient <- function(y, X, weight = 1) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    function(beta) {
      lambda <- exp(as.matrix(X) %*% beta)
      term <- lambda / (exp(lambda) - lambda - 1)

      t(as.matrix(X)) %*% (weight * (y - lambda - lambda * term))
    }
  }

  hessian <- function(y, X, weight = 1) {
    X <- as.matrix(X)
    y <- as.numeric(y)

    if (is.null(weight)) {
      weight <- 1
    }

    function(beta) {
      lambda <- exp(as.matrix(X) %*% beta)

      term <- (((2 + lambda ** 2) * exp(lambda) - exp(2 * lambda) - 1) /
                ((exp(lambda) - lambda - 1) ** 2))

      t(X) %*% as.matrix(t(t(as.data.frame(X) * lambda * term)))
    }
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  dev.resids <- function(y, mu, wt, disp = NULL) {
    eta <- log(mu)
    mu1 <- mu.eta(eta = eta, disp = disp)
    a <- function(y) {stats::uniroot(f = function(x) {mu.eta(x, disp = NULL) - y}, lower = -log(y), upper = y * 10, tol = .Machine$double.eps)$root}
    loghm1y <- y
    loghm1y[y == 2] <- -Inf
    loghm1y[y > 2] <- sapply(y[y > 2], FUN = a)
    loghm1y[y == 2] <- -10
    hm1y <- exp(loghm1y)
    log1mexphm1y <- log(1 - exp(-hm1y) - hm1y * exp(-hm1y))#ifelse(y > 2, log(1 - exp(-hm1y) - hm1y * exp(-hm1y)), 0)
    sign(y - mu1) * sqrt(-2 * wt * (y * eta - mu - log(1 - exp(-mu) - mu * exp(-mu)) - y * loghm1y + hm1y + log1mexphm1y))
  }

  pointEst <- function (disp = NULL, pw, lambda, contr = FALSE) {
    N <- (pw * (1 - lambda * exp(-lambda)) /
         (1 - exp(-lambda) - lambda * exp(-lambda)))
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (beta, pw, lambda, disp = NULL, cov, X) {
    X <- as.data.frame(X)
    prob <- (1 - exp(-lambda) - lambda * exp(-lambda))
    term <- (1 - lambda * exp(-lambda)) ** 2
    
    f1 <- t(X) %*% (as.numeric(pw * lambda * (1 - exp(lambda)) /
                                 ((1 + lambda - exp(lambda)) ** 2)))
    
    f1 <- t(f1) %*% as.matrix(cov) %*% f1
    
    f2 <- sum(pw * term * (1 - prob) / (prob ** 2))
    
    f1 + f2
  }

  ## is this a correct way to simulate data from this distribution?
  simulate <- function(n, lambda, lower=1) {
    NULL
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
      link = "log",
      valideta = function (eta) {TRUE},
      variance = variance,
      Wfun = Wfun,
      funcZ = funcZ,
      dev.resids = dev.resids,
      validmu = validmu,
      pointEst = pointEst,
      popVar= popVar,
      simulate = simulate,
      family = "zotpoisson"
    ),
    class = "family"
  )
}

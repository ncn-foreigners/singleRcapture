#' @rdname singleRmodels
#' @importFrom stats dbinom
#' @export
zelterman <- function(lambdaLink = "loghalf",
                      ...) {
  if (missing(lambdaLink)) lambdaLink <- "loghalf"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink)
  
  lambdaLink <- switch (lambdaLink,
    "loghalf" = singleRinternalloghalfLink
  )
  
  links[1] <- c(lambdaLink)
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    
    if (!deriv) {
      switch (type,
        "nontrunc" = lambda,
        "trunc" = (lambda / 2) / (1 + lambda / 2)
      )
    } else {
      switch (type,
        "nontrunc" = 1 * lambdaLink(eta, inverse = TRUE, deriv = 1),
        "trunc" = (2 / (lambda + 2) ^ 2) *
        lambdaLink(eta, inverse = TRUE, deriv = 1)
      )
    }
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = ((lambda / 2) / (1 + lambda / 2)) * (1 / (1 + lambda / 2))
    )
  }
  
  Wfun <- function(prior, eta, y, ...) {
    iddx <- y %in% 1:2
    lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
    z <- (lambda / 2) / (1 + lambda / 2)
    res <- iddx
    res[iddx] <- -prior[iddx] * (
      -2 / (lambda * (2 + lambda) ^ 2) *
      lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx] ^ 2)
    
    matrix(res, ncol = 1, 
           dimnames = list(rownames(eta), 
                           c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    iddx <- y %in% 1:2
    z <- y - 1
    lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
    
    res <- iddx
    res[iddx] <- prior[iddx] * (((z[iddx] - 1) * lambda + 2 * z[iddx]) / 
      (lambda * (2 + lambda))) * lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx] / weight[iddx, ]
    res
  }
  
  minusLogLike <- function(y, X, 
                           weight = 1, 
                           NbyK = FALSE, 
                           vectorDer = FALSE, 
                           deriv = 0,
                           offset,
                           eta, 
                           ...) {
    y <- as.numeric(y)
    z <- y - 1
    if (is.null(weight)) {
      weight <- 1
    }
    
    iddx <- y %in% 1:2
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X)))
    }
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (
      deriv,
      function(beta, eta) {
        if (missing(eta)) {
          eta <- as.matrix(X) %*% beta + offset
        }
        lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
        -sum(weight[iddx] * dbinom(size = 1, 
                                   prob = ((lambda / 2) / (1 + lambda / 2)),
                                   log  = TRUE,
                                   x    = z[iddx]))
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta + offset
        lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
        G0 <- iddx
        G0[iddx] <- ((z[iddx] - 1) * lambda + 2 * z[iddx]) / 
          (lambda * (2 + lambda)) * weight[iddx] * 
          lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx]
        
        #G0 <- (z - (lambda / 2) / (1 + lambda / 2)) * weight
        if (NbyK) {
          return(as.data.frame(X[iddx, , drop = FALSE]) * G0)
        }
        if (vectorDer) {
          return(matrix(G0, ncol = 1))
        }
        
        t(X) %*% G0
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta + offset
        lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
        
        G00 <- iddx
        G00[iddx] <- (-(4 * z[iddx] + 4 * z[iddx] * lambda + 
          (z[iddx] - 1) * lambda ^ 2) / ((2 + lambda) ^ 2 * lambda ^ 2)) *
          lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx] ^ 2 +
          ((z[iddx] - 1) * lambda + 2 * z[iddx]) / (lambda * (2 + lambda)) *
          lambdaLink(eta, inverse = TRUE, deriv = 2)[iddx]
        
        #G00 <- weight * (-lambda / 2 / (1 + lambda / 2) ^ 2)
        
        t(as.data.frame(X) * weight * G00) %*% as.matrix(X)
      }
    )
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(1 > mu)
  }
  
  devResids <- function(y, eta, wt, ...) {
    z <- y - 1
    lambda <- lambdaLink(eta, inverse = TRUE)
    
    ((-1) ^ y) * sqrt((-2 * wt * (z * log(lambda / 2 / (1 + lambda / 2)) + (1 - z) * log(1 / (1 + lambda / 2)))) * (y %in% 1:2))
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    N <- (pw * (1 / (1 - exp(-lambda))))
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    Xvlm <- as.matrix(Xvlm)
    prob <- 1 - exp(-lambda)

    f1 <- t((-lambdaLink(eta[,1], inverse = TRUE, deriv = 1) * pw * 
    as.numeric(exp(lambda) / (exp(lambda) - 1) ^ 2)) %*% Xvlm)
    
    f1 <- t(f1) %*% as.matrix(cov) %*% f1

    f2 <- sum(pw * (1 - prob) / (prob ^ 2))

    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta, inverse = TRUE)
    
    switch (type,
      "trunc" = {
        stats::dpois(x = x, lambda = lambda) / 
        (1 - stats::dpois(x = 0, lambda = lambda))
      },
      "nontrunc" = stats::dpois(x = x, lambda = lambda)
    )
  }

  simulate <- function(n, eta, lower = 0, upper = 2) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    lb <- stats::ppois(lower, lambda)
    ub <- stats::ppois(upper, lambda)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qpois(p_u, lambda)
    sims
  }
  
  getStart <- expression(
    if (method == "IRLS") {
      etaStart <- cbind(
        (priorWeights * (observed == 2) + .5) / (priorWeights + 1)
      ) + offset
      etaStart <- family$links[[1]](2 * (etaStart / (1 - etaStart)))
    } else if (method == "optim") {
      init <- mean(((priorWeights * (observed == 2) + .5) / (priorWeights + 1))[observed %in% 1:2])
      coefStart <- rep(family$links[[1]](2 * init / (1 - init)), NCOL(variables))
    }
  )
  
  structure(
    list(
      makeMinusLogLike = minusLogLike,
      densityFunction  = dFun,
      links     = links,
      mu.eta    = mu.eta,
      valideta  = function (eta) {TRUE},
      variance  = variance,
      Wfun      = Wfun,
      funcZ     = funcZ,
      devResids = devResids,
      validmu   = validmu,
      pointEst  = pointEst,
      popVar    = popVar,
      family    = "zelterman",
      etaNames  = c("lambda"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

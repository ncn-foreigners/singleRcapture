#' @rdname singleRmodels
#' @importFrom stats dbinom
#' @export
oichao <- function(lambdaLink = "log",
                 ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink)
  
  lambdaLink <- switch (lambdaLink,
                        "log" = singleRinternallogLink
  )
  links[1] <- c(lambdaLink)
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    if (!deriv) {
      switch(type,
             "nontrunc" = lambda,
             "trunc" = lambda / (3 + lambda))
    } else {
      switch(type,
             "nontrunc" = cbind(1 * lambdaLink(eta, inverse = TRUE, deriv = 1)),
             "trunc" = cbind((3 / (lambda + 3)^2) * lambdaLink(eta, inverse = TRUE, deriv = 1)))
    }
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    switch(type,
           "nontrunc" = lambda,
           "trunc" = (lambda / (3 + lambda)) * (3 / (3 + lambda)))
  }
  
  Wfun <- function(prior, eta, y, ...) {
    iddx <- y %in% 2:3
    lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
    z <- lambda / (3 + lambda)
    res <- iddx
    res[iddx] <- -prior[iddx] * (-9 / (lambda * (3 + lambda)^2) * lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx]^2)
    matrix(res, ncol = 1, dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    iddx <- y %in% 2:3
    z <- y - 2
    lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
    res <- iddx
    res[iddx] <- prior[iddx] * (((z[iddx] - 1) * lambda + 3 * z[iddx]) / (lambda * (3 + lambda))) * 
      lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx] / weight[iddx, ]
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
    z <- y - 2
    if (is.null(weight)) {
      weight <- 1
    }
    
    iddx <- y %in% 2:3
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X)))
    }
    
    if (!(deriv %in% c(0, 1, 2))) 
      stop("Only score function and derivatives up to 2 are supported.")
    
    # to make it conform to how switch in R works, i.e. indexing begins with 1
    deriv <- deriv + 1
    
    switch (deriv,
            function(beta, eta) {
              if (missing(eta)) {
                eta <- as.matrix(X) %*% beta + offset
              }
              lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
              -sum(weight[iddx] * dbinom(size = 1, 
                                         prob = (lambda / (3 + lambda)),
                                         log  = TRUE,
                                         x    = z[iddx]))
            },
            function(beta) {
              eta <- as.matrix(X) %*% beta + offset
              lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
              G0 <- iddx
              G0[iddx] <- ((z[iddx] - 1) * lambda + 3 * z[iddx]) / 
                (lambda * (3 + lambda)) * weight[iddx] * 
                lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx]
              
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
              G00[iddx] <- (-(9 * z[iddx] + 9 * z[iddx] * lambda + 
                                (z[iddx] - 1) * lambda ^ 2) / ((3 + lambda) ^ 2 * lambda ^ 2)) *
                lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx] ^ 2 +
                ((z[iddx] - 1) * lambda + 3 * z[iddx]) / (lambda * (3 + lambda)) *
                lambdaLink(eta, inverse = TRUE, deriv = 2)[iddx]
              
              t(as.data.frame(X) * weight * G00) %*% as.matrix(X)
            }
    )
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(1 > mu)
  }
  
  devResids <- function(y, eta, wt, ...) {
    z <- y - 2
    lambda <- lambdaLink(eta, inverse = TRUE)
    diff <- (z * log(lambda / (3 + lambda)) + (1 - z) * log(3 / (3 + lambda)))
    diff[diff > 0] <- -0
    ((-1)^y) * sqrt((-2 * wt * diff) * (y %in% 2:3))
  }
  
  # TODO: Check if N and observed_total are corresponding to theory
  
  pointEst <- function(pw, eta, contr = FALSE, y, ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    iddx <- y %in% 2:3
    f2 <- sum(pw[y == 2])
    f3 <- sum(pw[y == 3])
    
    # TODO: Possible version with bias-corrected estimator
    
    if (all(lambda == lambda[1])) {  # No covariate case
      observed_total <- sum(pw)
      f0_est <- if (f3 > 0) (2/9) * (f2^3 / f3^2) else 0
      N <- observed_total + f0_est
    } else {  # Covariate case
      N <- ((1 + iddx * 6 / (lambda^2 * (3 + lambda))) * pw)
      if (!contr) N <- sum(N)
    }
    N
  }
  
  # TODO: Check if variance is calculated correctly
  
  popVar <- function(pw, eta, cov, Xvlm, y, ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    iddx <- y %in% 2:3
    f2 <- sum(pw[iddx & y == 2])
    f3 <- sum(pw[iddx & y == 3])
    if (all(lambda == lambda[1])) {  # No covariates
      if (f3 == 0) return(0)
      variance <- (2/9)^2 * ( (3 * f2^2 / f3^2)^2 * f2 + (2 * f2^3 / f3^3)^2 * f3 )
    } else {  # Covariate case
      Xvlm <- as.matrix(Xvlm)
      prob <- (lambda^2 / 2) * exp(-lambda) + (lambda^3 / 6) * exp(-lambda)
      f1 <- t((-lambdaLink(eta[,1], inverse = TRUE, deriv = 1) * pw * 
                 as.numeric(6 * exp(-lambda) * (lambda^2 + 3 * lambda + 3) / 
                              (lambda^2 * (3 + lambda))^2))[iddx] %*% Xvlm[iddx, , drop = FALSE])
      f1 <- t(f1) %*% as.matrix(cov) %*% f1
      f2 <- sum((pw * (1 - prob) * ((1 + 6 * exp(-lambda) / prob / (lambda^2 * (3 + lambda))) ^ 2))[iddx])
      variance <- f1 + f2
    }
    variance
  }
  
  # TODO: Check what should the bounds be
  
  simulate <- function(n, eta, lower = 0, upper = 3) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    
    lb <- stats::ppois(lower, lambda)
    ub <- stats::ppois(upper, lambda)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qpois(p_u, lambda)
    sims
  }
  
  dFun <- function(x, eta, type = "trunc") {
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
  
  getStart <- expression(
    if (method == "IRLS") {
      etaStart <- cbind(
        (priorWeights * (observed == 3) + .5) / (priorWeights + 1)
      ) + offset
      etaStart <- family$links[[1]](etaStart / (3 * (1 - etaStart)))
    } else if (method == "optim") {
      init <- mean(((priorWeights * (observed == 3) + .5) / (priorWeights + 1))[observed %in% 2:3])
      coefStart <- rep(family$links[[1]](init / (3 * (1 - init))), NCOL(variables))
    }
  )
  
  structure(
    list(
      makeMinusLogLike = minusLogLike,
      densityFunction  = dFun,
      links     = links,
      mu.eta    = mu.eta,
      valideta  = function(eta) {TRUE},
      variance  = variance,
      Wfun      = Wfun,
      funcZ     = funcZ,
      devResids = devResids,
      validmu   = validmu,
      pointEst  = pointEst,
      popVar    = popVar,
      family    = "oichao",
      etaNames  = c("lambda"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}


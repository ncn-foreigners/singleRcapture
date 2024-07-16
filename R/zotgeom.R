#' @rdname singleRmodels
#' @export
zotgeom <- function(lambdaLink = c("log", "neglog"),
                    ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink)
  
  lambdaLink <- switch (lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  links[1] <- c(lambdaLink)
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    if (!deriv) {
      switch (type,
        "nontrunc" = lambda,
        "trunc" = 2 + lambda
      )
    } else {
      switch (type,
        "nontrunc" = cbind(lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)),
        "trunc" = cbind(lambdaLink(eta[, 1], inverse = TRUE, deriv = 1))
      )
    }
  }
  
  variance <- function(eta, disp, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    switch (type,
    "nontrunc" = lambda * (lambda - 1),
    "trunc" = lambda * (lambda + 1)
    )
  }
  
  Wfun <- function(prior, eta, y, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
    Ey <- mu.eta(eta)[iddx]
    
    G11 <- iddx
    G11[iddx] <- lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] ^ 2 * 
           (lambda ^ 2 + (4 - 2 * Ey) * lambda - Ey + 2) / 
           (lambda ^ 2 * (lambda + 1) ^ 2)
    
    matrix(-prior * G11, ncol = 1, dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
    
    res <- iddx
    res[iddx] <- -lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] *
    (lambda - y[iddx] + 2) / (lambda ^ 2 + lambda) * prior[iddx] / weight[iddx, ]
    res
  }
  
  minusLogLike <- function(y, X, 
                           weight    = 1, 
                           NbyK      = FALSE, 
                           vectorDer = FALSE, 
                           deriv     = 0,
                           offset, 
                           ...) {
    if (is.null(weight)) {
      weight <- 1
    }
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X)))
    }
    y <- as.numeric(y)
    iddx <- y > 1
    X <- as.matrix(X)
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- as.matrix(X) %*% beta + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        -sum(weight * iddx * ((y - 2) * log(lambda) - (y - 1) * log(1 + lambda)))
      },
      function(beta) {
        eta <- X %*% beta + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]

        G1 <- iddx
        G1[iddx] <- -lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] *
              (lambda - y[iddx] + 2) / (lambda ^ 2 + lambda)
        if (NbyK) {
          return(G1  * weight * as.data.frame(X))
        }
        if (vectorDer) {
          return(matrix(G1  * weight, ncol = 1))
        }
        t(G1  * weight) %*% X
      },
      function(beta) {
        eta <- X %*% beta + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
        
        G1 <- iddx
        G1[iddx] <- -lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)[iddx] *
              (lambda - y[iddx] + 2) / (lambda ^ 2 + lambda)
        
        G11 <- iddx
        G11[iddx] <- lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] ^ 2 * 
               (lambda ^ 2 + (4 - 2 * y[iddx]) * lambda - y[iddx] + 2) / 
               (lambda ^ 2 * (lambda + 1) ^ 2)
        
        t(as.data.frame(X) * (G1 + G11) * weight) %*% X
      }
    )
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  devResids <- function (y, eta, wt, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
    loghm1y <- ifelse(y[iddx] > 2, log(y[iddx] - 2), 0)
    
    diff <- iddx
    diff[iddx] <- (y[iddx]- 2) * log(lambda) - (y[iddx] - 1) * log(1 + lambda) - 
                  (y[iddx] - 2) * loghm1y + (y[iddx] - 1) * log(y[iddx] - 1)
    
    sign(y - mu.eta(eta = eta)) * sqrt(-2 * wt * diff)
  }
  
  pointEst <- function (pw, eta, contr = FALSE, y, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    N <- pw * (iddx * (lambda ^ 2 + lambda + 1) / lambda ^ 2 + (1 - iddx))
    if(!contr) {
      N <- sum(N)
    }
    
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, y, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    bigTheta <- t(Xvlm[iddx, , drop = FALSE]) %*% (pw * as.numeric(-(lambda + 2) / lambda ^ 3) *
                             lambdaLink(eta[, 1], inverse = TRUE, deriv = 1))[iddx]
    
    bigTheta <- as.vector(bigTheta)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    f2 <-  sum((pw * (lambda ^ 2 + lambda + 1) / lambda ^ 2)[iddx])
    
    f1 + f2
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    lb <- stats::pnbinom(lower, mu = lambda, size = 1)
    ub <- stats::pnbinom(upper, mu = lambda, size = 1)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qnbinom(p_u, mu = lambda, size = 1)
    sims
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        stats::dgeom(x = x, prob = 1 / (1 + lambda)) / 
        (1 - stats::dgeom(x = 0, prob = 1 / (1 + lambda)) - 
        stats::dgeom(x = 1, prob = 1 / (1 + lambda)))
      },
      "nontrunc" = stats::dgeom(x = x, prob = 1 / (1 + lambda))
    )
  }
  
  getStart <- expression(
    if (method == "IRLS") {
      etaStart <- cbind(
        pmin(family$links[[1]](observed), family$links[[1]](12))
      ) + offset
    } else if (method == "optim") {
      init <- c(
        family$links[[1]](weighted.mean(observed, priorWeights))
      )
      if (attr(terms, "intercept")) {
        coefStart <- c(init[1], rep(0, attr(Xvlm, "hwm")[1] - 1))
      } else {
        coefStart <- rep(init[1] / attr(Xvlm, "hwm")[1], attr(Xvlm, "hwm")[1])
      }
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
      family    = "zotgeom",
      etaNames  = c("lambda"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}
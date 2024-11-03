#' @rdname singleRmodels
#' @export
ztHurdlegeom <- function(lambdaLink = c("log", "neglog"), 
                         piLink = c("logit", "cloglog", "probit"), 
                         ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  if (missing(piLink))  piLink <- "logit"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink, piLink)
  
  lambdaLink <- switch(lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  piLink <- switch(piLink,
    "logit" = singleRinternallogitLink,
    "cloglog" = singleRinternalcloglogLink,
    "probit" = singleRinternalprobitLink
  )
  
  links[1:2] <- c(lambdaLink, piLink)
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    if (!deriv) {
      switch (type,
        "nontrunc" = (PI + (1 - PI) * (2 + lambda)) * lambda ^ 2 / 
          (lambda ^ 2 + lambda + 1),
        "trunc" = PI + (1 - PI) * (2 + lambda)
      )
    } else {
      switch (type,
        "nontrunc" = {
          matrix(c(
            -lambda * ((PI - 1) * lambda ^ 3 + (2 * PI - 2) * lambda ^ 2 + 
            (4 * PI - 5) * lambda + 2 * PI - 4) / (lambda ^ 2 + lambda + 1) ^ 2,
            -(1 + lambda) * lambda ^ 2 / (lambda ^ 2 + lambda + 1)
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
                piLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        },
        "trunc" = {
          matrix(c(
            1 - PI,
            1 - 2 - lambda
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
                piLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        }
      )
    }
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    switch (type,
    "nontrunc" = (1 - 1 / (1 + lambda)) * (PI + (1 - PI) * lambda + 2 * lambda ^ 2),
    "trunc" = (PI + (1 - PI) * (lambda * (1 + lambda) + lambda ^ 2 - lambda * (1 + lambda) ^ (-2)) / 
      (1 - 1 / (1 + lambda) - lambda * (1 + lambda) ^ (-2)))
    ) - mu.eta(eta = eta, type = type) ^ 2
  }
  
  Wfun <- function(prior, y, eta, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- PI
    YY <- mu.eta(eta) - z ## expected for (1-z)Y
    
    G00 <- (-1) * piLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 / ((1 - PI) * PI)
    
    G11 <- lambdaLink(eta[, 1], inverse = TRUE, deriv = 2) * 
      ((YY - (1 - z) * 2) / lambda - (YY - (1 - z)) / (lambda + 1)) -
      lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 *
      ((1 - z) / (lambda * (lambda + 1)) + 
      ((z - 1) * lambda + YY + (z - 1) * 2) / (lambda ^ 2 * (lambda + 1)) + 
      ((z - 1) * lambda + YY + (z - 1) * 2) / (lambda * (lambda + 1) ^ 2))
    
    matrix(
      -c(G11 * prior, # lambda
         rep(0, NROW(eta)) * prior, # mixed
         rep(0, NROW(eta)) * prior, # mixed
         G00 * prior  # pi
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "pi")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- as.numeric(y == 1)
    weight <- weight / prior
    
    G1 <- lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
      (1 - z) * ((y - 2) / lambda - (y - 1) / (lambda + 1))
    
    G0 <- piLink(eta[, 1], inverse = TRUE, deriv = 1) *  
      (z / PI - (1 - z) / (1 - PI))
    
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- diag(weight[[x]]) # in this case matrix is diagonal
      xx * uMatrix[x, ]
    })
    pseudoResid <- t(pseudoResid)
    dimnames(pseudoResid) <- dimnames(eta)
    pseudoResid
  }
  
  minusLogLike <- function(y, X, 
                           weight    = 1, 
                           NbyK      = FALSE, 
                           vectorDer = FALSE, 
                           deriv     = 0,
                           offset, 
                           ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X) / 2), rep(0, NROW(X) / 2))
    }
    z <- as.numeric(y == 1)
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        PI     <- piLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        -sum(weight * (z * log(PI) + (1 - z) * log(1 - PI) +
        (1 - z) * ((y - 2) * log(lambda) - (y - 1) * log(1 + lambda))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        PI     <- piLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- weight * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
              (1 - z) * ((y - 2) / lambda - (y - 1) / (lambda + 1))
        
        G0 <- weight * piLink(eta[, 1], inverse = TRUE, deriv = 1) *  
              (z / PI - (1 - z) / (1 - PI))
        
        if (NbyK) {
          XX <- 1:(attr(X, "hwm")[1])
          return(cbind(as.data.frame(X[1:nrow(eta), XX]) * G1, as.data.frame(X[-(1:nrow(eta)), -XX]) * G0))
        }
        if (vectorDer) {
          return(cbind(G1, G0))
        }
        as.numeric(c(G1, G0) %*% X)
      },
      function (beta) {
        lambdaPredNumber <- attr(X, "hwm")[1]
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        PI     <- piLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        res <- matrix(nrow = length(beta), ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        # PI^2 derivative
        G00 <- piLink(eta[, 1], inverse = TRUE, deriv = 2) *  
               (z / PI - (1 - z) / (1 - PI)) -
               piLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 *
               (1 / ((1 - PI) * PI) + (PI - z) / ((1 - PI) ^ 2 * PI) + 
               (z - PI) / ((1 - PI) * PI ^ 2))
        
        # Beta^2 derivative
        G11 <- lambdaLink(eta[, 1], inverse = TRUE, deriv = 2) * 
          (1 - z) * ((y - 2) / lambda - (y - 1) / (lambda + 1)) -
          lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 *
          ((1 - z) / (lambda * (lambda + 1)) + 
          ((z - 1) * (lambda - y + 2)) / (lambda ^ 2 * (lambda + 1)) + 
          ((z - 1) * (lambda - y + 2)) / (lambda * (lambda + 1) ^ 2))
        
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)] * G00 * weight)) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * G11 * weight)) %*% X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- 0
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- 0
        
        res
      }
    )
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  devResids <- function(y, eta, wt, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    # idealPI <- ifelse(y == 1, 1, 0) memmory allocation not needed when pi = 0 distribution collapses to zotgeom
    idealLambda <- ifelse(y > 1, y - 2, 0)
    diff <- ifelse(
      y == 1, -log(PI),
      ifelse(y == 2, 0,
      (y - 2) * log(idealLambda) - (y - 1) * log(1 + idealLambda)) - 
      (log(1 - PI) + (y - 2) * log(lambda) - (y - 1) * log(1 + lambda))
    )
    
    diff[diff < 0] <- 0
    sign(y - mu.eta(eta = eta)) * sqrt(2 * wt * diff)
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    N <- pw * (lambda ^ 2 + lambda + 1) / lambda ^ 2
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    bigTheta1 <- rep(0, nrow(eta)) # w.r to PI
    bigTheta2 <- (pw * as.numeric(-(lambda + 2) / lambda ^ 3) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)) # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * (lambda ^ 2 + lambda + 1) / lambda ^ 2)
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        as.numeric(x == 1) * PI + as.numeric(x > 0) * 
        (1 - PI) * stats::dnbinom(x = x, mu = lambda, size = 1) / 
        (1 - 1 / (1 + lambda) - lambda * ((1 + lambda) ^ -2))
      },
      "nontrunc" = {
        as.numeric(x == 0) * (1 + lambda) / (lambda ^ 2 + lambda + 1) +
        as.numeric(x == 1) * PI * lambda ^ 2 / (lambda ^ 2 + lambda + 1) +
        as.numeric(x > 1) * (1 - PI) * (lambda ^ x / (1 + lambda) ^ (x + 1)) *
        (lambda + 1) ^ 2 / (lambda ^ 2 + lambda + 1)
      }
    )
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    CDF <- function(x) {
      p <- lambda / (1 + lambda)
      const <- -p * (lambda * (p ^ x - 1) + p ^ x) / (1 + lambda)
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, (1 + lambda) / (1 + lambda + lambda ^ 2), 
      (1 + lambda) / (1 + lambda + lambda ^ 2) + PI * (lambda ^ 2
      ) / (1 + lambda + lambda ^ 2) +  (1 - PI) * ((1 + lambda) ^ 2
      ) * const / (lambda ^ 2 + lambda + 1))))
    }
    lb <- CDF(lower)
    ub <- CDF(upper)
    p_u <- stats::runif(n, lb, ub)
    sims <- rep(0, n)
    cond <- CDF(sims) <= p_u
    while (any(cond)) {
      sims[cond] <- sims[cond] + 1
      cond <- CDF(sims) <= p_u
    }
    sims
  }
  
  getStart <- expression(
    if (method == "IRLS") {
      etaStart <- cbind(
        pmin(family$links[[1]](observed), family$links[[1]](12)),
        family$links[[2]](weighted.mean(observed == 1, priorWeights) * (.5 + .5 * (observed == 1)) + .01)
      ) + offset
    } else if (method == "optim") {
      init <- c(
        family$links[[1]](weighted.mean(observed, priorWeights)),
        family$links[[2]](weighted.mean(observed == 1, priorWeights) + .01)
      )
      if (attr(terms, "intercept")) {
        coefStart <- c(init[1], rep(0, attr(Xvlm, "hwm")[1] - 1))
      } else {
        coefStart <- rep(init[1] / attr(Xvlm, "hwm")[1], attr(Xvlm, "hwm")[1])
      }
      if ("(Intercept):pi" %in% colnames(Xvlm)) {
        coefStart <- c(coefStart, init[2], rep(0, attr(Xvlm, "hwm")[2] - 1))
      } else {
        coefStart <- c(coefStart, rep(init[2] / attr(Xvlm, "hwm")[2], attr(Xvlm, "hwm")[2]))
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
      family    = "ztHurdlegeom",
      etaNames = c("lambda", "pi"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

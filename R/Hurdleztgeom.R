#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @export
Hurdleztgeom <- function(lambdaLink = c("log", "neglog"), 
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
        "nontrunc" = PI + (1 - PI) * lambda * lambda * (2 + lambda) / 
          (lambda ^ 2 + lambda + 1),
        "trunc" = PI * (lambda ^ 2 + lambda + 1) / 
          (lambda ^ 2 + PI * (lambda + 1)) + 
          (1 - PI) * (2 + lambda) * lambda ^ 2 / (lambda ^ 2 + PI * (lambda + 1))
      )
    } else {
      switch (type,
        "nontrunc" = {
          matrix(c(
            (1 - PI) * lambda * (lambda ^ 3 + 2 * lambda ^ 2 + 5 * lambda + 4) /
            (lambda ^ 2 + lambda + 1) ^ 2,
            1 - (lambda ^ 2 * (lambda + 2)) / (lambda ^ 2 + lambda + 1)
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
                piLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        },
        "trunc" = {
          matrix(c(
            (1 - PI) * lambda * 
            (lambda ^ 3 + 2 * PI * lambda ^ 2 + 4 * PI * lambda + 2 * PI) /
            (lambda ^ 2 + PI * lambda + PI) ^ 2,
            -lambda ^ 2 * (lambda + 1) * (lambda ^ 2 + lambda + 1) /
            ((lambda + 1) * PI + lambda ^ 2) ^ 2
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
    switch(type,
    "nontrunc" = PI + (1 - PI) * lambda * lambda * (2 * lambda * lambda + 5 * lambda + 4) / (lambda ^ 2 + lambda + 1),
    "trunc" = PI * (lambda ^ 2 + lambda + 1) / (lambda ^ 2 + PI * (lambda + 1)) + (1 - PI) * lambda * lambda * (2 * lambda * lambda + 5 * lambda + 4) / (lambda ^ 2 + PI * (lambda + 1))
    ) - (mu.eta(eta = eta, type = type) ^ 2)
  }
  
  Wfun <- function(prior, eta, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- PI * (lambda ^ 2 + lambda + 1) / (lambda ^ 2 + PI * (lambda + 1))
    #Ey <- mu.eta(eta)
    YY <- mu.eta(eta) - z ## expected for (1-z)Y
    
    G1 <- z * (2 * lambda + 1) / (lambda ^ 2 + lambda + 1) + 
      (1 - z - YY) / (lambda + 1) + YY / lambda - 
      (2 * lambda + PI) / (lambda ^ 2 + PI * (lambda + 1)) # lambda derivative
    
    G0 <- z / PI - (1 - z) / (1 - PI) - 
      (1 + lambda) / (lambda ^ 2 + PI * (1 + lambda)) # PI derivative
    
    # PI^2 derivative
    G00 <- -z / (PI ^ 2) - (1 - z) / ((1 - PI) ^ 2) + 
      ((lambda + 1) / (lambda ^ 2 + PI * (lambda + 1))) ^ 2
    G00 <- prior * (G0 * (PI * (1 - PI) * (1 - 2 * PI)) + G00 * ((PI * (1 - PI)) ^ 2))
    
    # mixed
    
    G01 <- lambda * (lambda + 2) / ((lambda ^ 2 + PI * (lambda + 1)) ^ 2) * 
           lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
           piLink(eta[, 2], inverse = TRUE, deriv = 1) * prior
    
    # Beta^2 derivative
    G11 <- ((YY - 1 + z) / (1 + lambda) ^ 2 - YY / lambda ^ 2) + 
      z * (2 * (lambda ^ 2 + lambda + 1) - (2 * lambda + 1) ^ 2) / 
      (lambda ^ 2 + lambda + 1) ^ 2 + 
      (2 * lambda + PI) ^ 2 / (lambda ^ 2 + PI * (lambda + 1)) ^ 2 - 
      2 / (lambda ^ 2 + PI * (lambda + 1))
    
    G11 <- (G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 + 
            G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)) * prior
    
    matrix(
      -c(G11, # lambda
         G01, # mixed
         G01, # mixed
         G00  # pi
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "pi")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- ifelse(y == 1, y, 0)
    
    G1 <- z * (2 * lambda + 1) / (lambda ^ 2 + lambda + 1) + 
      (1 - z) * (y / lambda - (y - 1) / (1 + lambda)) - 
      (2 * lambda + PI) / (lambda ^ 2 + PI * (lambda + 1))
    G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) # lambda derivative
    
    G0 <- z / PI - (1 - z) / (1 - PI) - 
      (1 + lambda) / (lambda ^ 2 + PI * (1 + lambda))
    G0 <- G0 * piLink(eta[, 2], inverse = TRUE, deriv = 1) # PI derivative
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- solve(weight[[x]]) # more stable
      xx %*% uMatrix[x, ]
    })
    pseudoResid <- t(pseudoResid)
    dimnames(pseudoResid) <- dimnames(eta)
    pseudoResid
  }
  
  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    z <- as.numeric(y == 1)
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        PI     <- piLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        -sum(weight * (z * (log(PI) + log(lambda ^ 2 + lambda + 1)) + 
        (1 - z) * (log(1 - PI) + y * log(lambda) - (y - 1) * log(1 + lambda)) - 
        log(lambda ^ 2 + PI * (lambda + 1))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        PI     <- piLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- z * (2 * lambda + 1) / (lambda ^ 2 + lambda + 1) + 
              (1 - z) * (y / lambda - (y - 1) / (1 + lambda)) - 
              (2 * lambda + PI) / (lambda ^ 2 + PI * (lambda + 1))
        G1 <- G1 * weight * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) # lambda derivative
        
        G0 <- z / PI - (1 - z) / (1 - PI) - 
              (1 + lambda) / (lambda ^ 2 + PI * (1 + lambda))
        G0 <- G0 * weight * piLink(eta[, 2], inverse = TRUE, deriv = 1) # PI derivative
        
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
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        PI     <- piLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        res <- matrix(nrow = length(beta), ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        G1 <- z * (2 * lambda + 1) / (lambda ^ 2 + lambda + 1) + 
          (1 - z) * (y / lambda - (y - 1) / (1 + lambda)) - 
          (2 * lambda + PI) / (lambda ^ 2 + PI * (lambda + 1))
        G1 <- G1 * weight * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) # lambda derivative
        
        G0 <- z / PI - (1 - z) / (1 - PI) - 
          (1 + lambda) / (lambda ^ 2 + PI * (1 + lambda))
        G0 <- G0 * weight * piLink(eta[, 2], inverse = TRUE, deriv = 1) # PI derivative
        
        # PI^2 derivative
        G00 <- -z / (PI ^ 2) - (1 - z) / ((1 - PI) ^ 2) + 
          ((lambda + 1) ^ 2) / ((lambda ^ 2 + PI * (lambda + 1)) ^ 2)
        G00 <- weight * (G0 * piLink(eta[, 2], inverse = TRUE, deriv = 2) + 
                         G00 * piLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2)  # second derivative of inverse logistic link
        
        # mixed
        
        G01 <- lambda * (lambda + 2) / ((lambda ^ 2 + PI * (lambda + 1)) ^ 2)
        G01 <- G01 * piLink(eta[, 2], inverse = TRUE, deriv = 1) *
               lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * weight
        
        # Beta^2 derivative
        G11 <- (1 - z) * ((y - 1) / (1 + lambda) ^ 2 - y / lambda ^ 2) + 
          z * (2 * (lambda ^ 2 + lambda + 1) - ((2 * lambda + 1) ^ 2)) /
          (lambda ^ 2 + lambda + 1) ^ 2 + 
          (2 * lambda + PI) ^ 2 / (lambda ^ 2 + PI * (lambda + 1)) ^ 2 - 
          2 / (lambda ^ 2 + PI * (lambda + 1))
        
        G11 <- (G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 + 
                G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)) * weight # second derivative of log link
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)] * G00)) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * G11)) %*% X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * G01)) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- t(t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * G01)) %*% as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]))
        
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
    #idealPI <- ifelse(y == 1, 1, 0) memmory allocation not needed
    # when pi = 0 distribution collapses to zotgeom
    idealLambda <- ifelse(y > 1, y - 2, 0)
    diff <- ifelse(
      y == 1, -(log(PI) + log(lambda ^ 2 + lambda + 1) - log(lambda ^ 2 + PI * (lambda + 1))),
      ifelse(y == 2, 0,
      (y - 2) * log(idealLambda) - (y - 1) * log(1 + idealLambda)) - (log(1 - PI) + y * log(lambda) - (y - 1) * log(1 + lambda) - log(lambda ^ 2 + PI * (lambda + 1)))
    )
    sign(y - mu.eta(eta = eta)) * sqrt(2 * wt * diff)
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    N <- pw * (lambda ^ 2 + lambda + 1) / (lambda ^ 2 + PI * (lambda + 1))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    bigTheta1 <- -pw *  piLink(eta[, 2], inverse = TRUE, deriv = 1) * 
      (lambda ^ 2 + lambda + 1) * (lambda + 1) / 
      (lambda ^ 2 + PI * (lambda + 1)) ^ 2 # w.r to PI
    bigTheta2 <- pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
      ((PI - 1) * lambda * (lambda + 2)) / (lambda ^ 2 + PI * lambda + PI) ^ 2 # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * (lambda ^ 2 + lambda + 1) * (1 - PI) * (1 + lambda) / ((lambda ^ 2 + PI * (lambda + 1)) ^ 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        ifelse(x == 1, PI * (lambda ^ 2 + lambda + 1), 
        (1 - PI) * (lambda ^ x) / ((1 + lambda) ^ (x - 1))) / 
        (lambda ^ 2 + PI * (lambda + 1))
      },
      "nontrunc" = {
        ifelse(x == 1, PI, (1 - PI) * 
        (lambda ^ x / (1 + lambda) ^ (x - 1)) / (lambda ^ 2 + lambda + 1))
      }
    )
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    CDF <- function(x) {
      p <- lambda / (1 + lambda)
      const <- -lambda * (p ^ x + lambda * (p ^ x - 1))
      polly <- lambda ^ 2 + lambda + 1
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, (1 - PI) * (1 + lambda) / polly, 
      (1 - PI) * (1 + lambda) / polly + PI + (1 - PI) * const / polly)))
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
    start <- stats::glm.fit(
      x = variables[wch$reg, 1:attr(Xvlm, "hwm")[1]],
      y = observed[wch$reg],
      family = stats::poisson(),
      weights = priorWeights[wch$reg],
      ...
    )$coefficients,
    if (attr(family$links, "linkNames")[1] == "neglog") start <- -start,
    if (is.null(controlMethod$piStart)) {
      cc <- colnames(Xvlm)
      cc <- cc[grepl(x = cc, pattern = "pi$")]
      cc <- unlist(strsplit(x = cc, ":pi"))
      cc <- sapply(cc, FUN = function(x) {
        ifelse(x %in% names(start), start[x], 0) # TODO: gosh this is terrible pick a better method
      })
      start <- c(start, cc)
    } else {
      start <- c(start, controlMethod$piStart)
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
      family    = "Hurdleztgeom",
      etaNames = c("lambda", "pi"),
      simulate  = simulate,
      getStart  = getStart,
      extraInfo = c(
        mean       = paste0(
          "PI + (1 - PI) * lambda * lambda * ",
          "(2 + lambda) / (lambda ^ 2 + lambda + 1)"
        ),
        variance   = paste0(
          "PI + (1 - PI) * lambda * lambda * ",
          "\n(2 * lambda * lambda + 5 * lambda + 4) / (lambda ^ 2 + lambda + 1)",
          " - mean ^ 2"
        ),
        popSizeEst = paste0(
          "(lambda ^ 2 + lambda + 1) / ",
          "(lambda ^ 2 + PI * (lambda + 1))"
        ),
        meanTr     = paste0(
          "PI * (lambda ^ 2 + lambda + 1) / (lambda ^ 2 + PI * (lambda + 1)) + ",
          "\n(1 - PI) * (2 + lambda) * lambda ^ 2 / (lambda ^ 2 + PI * (lambda + 1))"
        ),
        varianceTr = paste0(
          "PI * (lambda ^ 2 + lambda + 1) / (lambda ^ 2 + PI * (lambda + 1)) +",
          "\n(1 - PI) * lambda * lambda * (2 * lambda * lambda + 5 * lambda + 4)",
          " / (lambda ^ 2 + PI * (lambda + 1))",
          " - meanTr ^ 2"
        )
      )
    ),
    class = c("singleRfamily", "family")
  )
}

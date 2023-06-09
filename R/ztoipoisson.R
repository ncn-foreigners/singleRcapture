#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @export
ztoipoisson <- function(lambdaLink = c("log", "neglog"), 
                        omegaLink = c("logit", "cloglog", "probit"), 
                        ...) {
  if (missing(lambdaLink)) lambdaLink <- "log"
  if (missing(omegaLink))  omegaLink <- "logit"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink, omegaLink)

  lambdaLink <- switch(lambdaLink,
    "log"    = singleRinternallogLink,
    "neglog" = singleRinternalneglogLink
  )
  
  omegaLink <- switch(omegaLink,
    "logit" = singleRinternallogitLink,
    "cloglog" = singleRinternalcloglogLink,
    "probit" = singleRinternalprobitLink
  )
  
  links[1:2] <- c(lambdaLink, omegaLink)
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    if (!deriv) {
      switch (type,
        "nontrunc" = omega * (1 - exp(-lambda)) + lambda * (1 - omega),
        "trunc" = omega + (1 - omega) * lambda / (1 - exp(-lambda))
      )
    } else {
      switch (type,
        "nontrunc" = {
          matrix(c(
            1 + omega * exp(-lambda) - omega,
            1 - lambda - exp(-lambda)
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
             omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        },
        "trunc" = {
          matrix(c(
            (1 - omega) * exp(lambda) * (exp(lambda) - lambda - 1) /
            (exp(lambda) - 1) ^ 2,
            1 - lambda / (1 - exp(-lambda))
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
             omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        }
      )
    }
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
    "nontrunc" = omega * (1 - exp(-lambda)) + (1 - omega) * (lambda + lambda ^ 2),
    "trunc" = omega + (1 - omega) * (lambda ^ 2 + lambda) / (1 - exp(-lambda))
    ) - mu.eta(type = type, eta = eta) ^ 2
  }
  
  Wfun <- function(prior, y, eta, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- omega + (1 - omega) * lambda / (exp(lambda) - 1) # expected for I's
    #z <- ifelse(y == 1, y, 0)
    
    G00 <- prior * ((-(z * (1 - lambda / (exp(lambda) - 1)) ^ 2) / 
    (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) ^ 2 - 
    (1 - z) / (1 - omega) ^ 2) * (omegaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2) + 
    ((z * (1 - lambda / (exp(lambda) - 1))) / 
    (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) - 
    (1 - z) / (1 - omega)) * omegaLink(eta[, 2], inverse = TRUE, deriv = 2))
    
    # mixed derivative
    G01 <- prior * (((z * ((lambda - 1) * exp(lambda) + 1)) / 
    (omega * exp(lambda) + (1 - omega) * lambda - omega) ^ 2) * 
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
    omegaLink(eta[, 2], inverse = TRUE, deriv = 1))
    
    #expected value of (1-z)*y
    # XXXX <- (1 - omega) * lambda * (exp(lambda) - 1) / (exp(lambda) - 1)
    XXXX <- (1 - omega) * lambda
    
    G11 <- prior * (lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)  * 
    ((1 - z) * (-exp(lambda) / (exp(lambda) - 1)) + XXXX / lambda +
    (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) /
    (exp(lambda) - 1) ^ 2)) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega)) +
    (lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2) * 
    ((1 - z) * (exp(2 * lambda) / (exp(lambda) - 1) ^ 2 - exp(lambda) / (exp(lambda) - 1))-
    XXXX / lambda ^ 2 + (z * ((2 * (1 - omega) * lambda * exp(2 * lambda)) / (exp(lambda) - 1) ^ 3-
    ((1 - omega) * lambda * exp(lambda)) / (exp(lambda) - 1) ^ 2-
    (2 * (1 - omega) * exp(lambda)) / (exp(lambda) - 1) ^ 2)) /
    (((1 - omega) * lambda) / (exp(lambda) - 1) + omega) - 
    (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda))/
    (exp(lambda) - 1) ^ 2) ^ 2) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega) ^ 2))
    
    matrix(
      -c(G11, # lambda
         G01, # mixed
         G01, # mixed
         G00  # omega
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "omega")),
      ncol = 4
    )
  }
  
  funcZ <- function(eta, weight, y, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- ifelse(y == 1, y, 0)
    G0 <- (z * (1 - lambda / (exp(lambda) - 1))) / (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) - (1 - z) / (1 - omega)
    
    G1 <- ((1 - z) * (y / lambda - exp(lambda) / (exp(lambda) - 1)) + 
    (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) / 
    (exp(lambda) - 1) ^ 2)) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega))
    
    G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    G0 <- G0 *  omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
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
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        -sum(weight * (-log(1 + exp(eta[, 2])) + z * log(exp(eta[, 2]) + lambda / (exp(lambda) - 1)) +
        (1 - z) * (y * log(lambda) - log(exp(lambda) - 1) - lgamma(y + 1))))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2)
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G0 <- (z * (1 - lambda / (exp(lambda) - 1))) / (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) - (1 - z) / (1 - omega)
        
        G1 <- ((1 - z) * (y / lambda - exp(lambda) / (exp(lambda) - 1)) + 
        (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) / 
        (exp(lambda) - 1) ^ 2)) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega))
        
        G1 <- G1 * weight * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
        G0 <- G0 * weight *  omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
        
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
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)

        res <- matrix(nrow = length(beta), ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        # omega^2 derivative
        domega <- (z * (1 - lambda / (exp(lambda) - 1))) / 
        (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) - 
        (1 - z) / (1 - omega)
        G00 <- -(z * (1 - lambda / (exp(lambda) - 1)) ^ 2) / 
        (omega + (lambda * (1 - omega)) / (exp(lambda) - 1)) ^ 2 - 
        (1 - z) / (1 - omega) ^ 2
        
        G00 <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)] * 
        (G00 * (omegaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2) + 
        domega * omegaLink(eta[, 2], inverse = TRUE, deriv = 2)) * weight)) %*% 
        as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        # mixed derivative
        
        G01 <- ((z * ((lambda - 1) * exp(lambda) + 1)) / 
        (omega * exp(lambda) + (1 - omega) * lambda - omega) ^ 2)
        
        G01 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * 
        G01 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
        omegaLink(eta[, 2], inverse = TRUE, deriv = 1) * weight) %*% 
        as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        # Beta^2 derivative
        G11 <- ((1 - z) * (exp(2 * lambda) / (exp(lambda) - 1) ^ 2-
        exp(lambda) / (exp(lambda) - 1) - y / lambda ^ 2) + 
        (z * ((2 * (1 - omega) * lambda * exp(2 * lambda)) / (exp(lambda) - 1) ^ 3-
        ((1 - omega) * lambda * exp(lambda)) / (exp(lambda) - 1) ^ 2 - 
        (2 * (1 - omega) * exp(lambda)) / (exp(lambda) - 1) ^ 2)) / 
        (((1 - omega) * lambda) / (exp(lambda) - 1) + omega) - 
        (z * ((1 - omega) /(exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) /
        (exp(lambda) - 1) ^ 2) ^ 2) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega) ^ 2)
        
        dlambda <- ((1 - z) * (y / lambda - exp(lambda) / (exp(lambda) - 1)) + 
        (z * ((1 - omega) / (exp(lambda) - 1) - ((1 - omega) * lambda * exp(lambda)) / 
        (exp(lambda) - 1) ^ 2)) / (((1 - omega) * lambda) / (exp(lambda) - 1) + omega))
        
        G11 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * 
        (G11 * (lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2) + 
        dlambda * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)) * weight)) %*% 
        X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(G01)
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- G01
        
        res
      }
    )
  }
  
  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }
  
  devResids <- function(y, eta, wt, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    mu <- mu.eta(eta = eta)
    #idealOmega <- ifelse(y == 1, 1, 0)
    
    idealLambda <- tryCatch(
      expr = {
        suppressWarnings(
          ifelse(y > 1, lamW::lambertW0(-y * exp(-y)) + y, 0)
        )
      },
      error = function (e) {
        warning("Deviance residuals could not have been computed and zero vector will be returned instead.", call. = FALSE)
        NULL
      }
    )
    if (is.null(idealLambda)) {
      return(rep(0, length(y)))
    }
    
    diff <- ifelse(
      y == 1,
      -(log(omega + (1 - omega) * lambda / (exp(lambda) - 1))),
      y * (log(idealLambda) - log(lambda)) + log((exp(lambda) - 1) / (exp(idealLambda) - 1)) - log(1 - omega)
    )
    
    if (any(diff < 0)) {
      warning(paste0(
        "Some of differences between log likelihood in sautrated model",
        " and fitted model were positive which indicates either:\n",
        "(1): A very good model fitt or\n",
        "(2): Incorrect computation of saturated model",
        "\nDouble check deviance before proceeding"
      ))
    }
    
    sign(y - mu) * sqrt(2 * wt * diff)
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    N <- pw / (1 - exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    bigTheta1 <- rep(0, nrow(eta)) # w.r to omega
    bigTheta2 <- pw * (exp(lambda) / (exp(lambda) - 1) ^ 2) # w.r to lambda
    bigTheta2 <- bigTheta2 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)# w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <- t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * exp(-lambda) / (1 - exp(-lambda)) ^ 2)
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        (1 - omega) * (lambda ^ x) / (factorial(x) * (exp(lambda) - 1)) +
        omega * as.numeric(x == 1)
      },
      "nontrunc" = {
        stats::dpois(x = x, lambda = lambda) * 
        (as.numeric(x == 0) + as.numeric(x > 0) * (1 - omega)) +
        omega * (1 - exp(-lambda)) * as.numeric(x == 1)
      }
    )
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    CDF <- function(x) {
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, exp(-lambda), 
      exp(-lambda) + omega * (1 - exp(-lambda)) + 
      (1 - omega) * (stats::ppois(x, lambda) - exp(-lambda)))))
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
  
  # new
  getStart <- expression(
    if (!is.null(controlMethod$start)) {
      start <- controlMethod$start
    } else {
      init <- c(
        family$links[[1]](mean(observed)),
        family$links[[2]](mean(observed == 1) + .01)
      )
      if (attr(terms, "intercept")) {
        start <- c(init[1], rep(0, attr(Xvlm, "hwm")[1] - 1))
      } else {
        start <- rep(init[1] / attr(Xvlm, "hwm")[1], attr(Xvlm, "hwm")[1])
      }
      if ("(Intercept):omega" %in% colnames(Xvlm)) {
        start <- c(start, init[2], rep(0, attr(Xvlm, "hwm")[2] - 1))
      } else {
        start <- c(start, rep(init[2] / attr(Xvlm, "hwm")[2], attr(Xvlm, "hwm")[2]))
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
      family    = "ztoipoisson",
      etaNames  = c("lambda", "omega"),
      simulate  = simulate,
      getStart  = getStart,
      extraInfo = c(
        mean       = "omega * (1 - exp(-lambda)) + lambda * (1 - omega)",
        variance   = paste0(
          "omega * (1 - exp(-lambda))",
          " + (1 - omega) * (lambda + lambda ^ 2) - mean ^ 2"
        ),
        popSizeEst = "(1 + exp(-lambda)) ^ -1",
        meanTr     = "omega + (1 - omega) * lambda / (1 - exp(-lambda))",
        varianceTr = paste0(
          "omega + (1 - omega) * (lambda ^ 2 + lambda)",
          " / (1 - exp(-lambda)) - meanTr ^ 2"
        )
      )
    ),
    class = c("singleRfamily", "family")
  )
}

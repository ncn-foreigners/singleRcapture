#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @export
oiztpoisson <- function(lambdaLink = c("log", "neglog"), 
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
        "nontrunc" = omega + lambda * (1 - omega),
        "trunc" = exp(lambda) * (omega + lambda - omega * lambda) / (exp(lambda) - 1 + omega)
      )
    } else {
      switch (type,
        "nontrunc" = {
          matrix(c(
            1 - omega, 
            1 - lambda
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
            omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        },
        "trunc" = {
          matrix(c(
            (1 - omega) * exp(lambda) * (exp(lambda) + (omega - 1) * lambda - 1) /
            (exp(lambda) + omega - 1) ^ 2, 
            -exp(lambda) * ((lambda - 1) * exp(lambda) + 1) /
            (omega + exp(lambda) - 1) ^ 2
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
    "nontrunc" = omega + (1 - omega) * lambda * (lambda + 1),
    "trunc" = ((exp(lambda) * omega + lambda * (1 - omega)) / (exp(lambda) - 1 + omega) +
    (1 - omega) * lambda * ((1 + lambda) * exp(lambda) - 1) / (exp(lambda) - 1 + omega))
    ) - mu.eta(type = type, eta = eta) ^ 2
  }
  
  Wfun <- function(prior, y, eta, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    z <- (exp(lambda) * omega + lambda * (1 - omega)) / (exp(lambda) - 1 + omega) # expected for I's
    #z <- ifelse(y == 1, y, 0)
    Ey <- mu.eta(eta)
    XXX <- Ey - z ## z here is the prob of 1 and XXX is expected for (1-z)*y
    
    # omega^2 derivative
    G00 <- (-(z * (exp(lambda) - lambda) ^ 2) /
    (exp(lambda) * omega + lambda * (1 - omega)) ^2 +
    1 / (omega + exp(lambda) - 1) ^ 2 - (1 - z) / (1 - omega) ^ 2)
    
    # mixed derivative
    G01 <- (z * (exp(lambda) - 1)) / (omega * exp(lambda) + (1 - omega) * lambda) -
    (z * (exp(lambda) - lambda) * (omega * exp(lambda) - omega + 1)) /
    (omega * exp(lambda) + (1 - omega) * lambda) ^ 2 + 
    exp(lambda) / (exp(lambda) + omega - 1) ^ 2
    
    # Beta^2 derivative
    G11 <- (exp(2 * lambda) / (exp(lambda) + omega - 1) ^ 2 +
    (omega * z * exp(lambda)) / (omega * exp(lambda) + (1 - omega) * lambda) -
    (z * (omega * exp(lambda) - omega + 1) ^ 2) / 
    (omega * exp(lambda) + (1 - omega) * lambda) ^ 2 -
    exp(lambda) / (exp(lambda) + omega - 1) - XXX / lambda ^ 2)
    
    G00 <- prior * G00 * omegaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2
    
    G01 <- prior * G01 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
           omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    G11 <- prior * G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
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
  
  funcZ <- function(eta, weight, y, prior, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- ifelse(y == 1, y, 0)
    weight <- weight / prior
    
    G1 <- (z  * (exp(lambda) * omega + 1 - omega) / 
    (exp(lambda) * omega + lambda - omega * lambda) + 
    (1 - z) * y / lambda - exp(lambda) / (exp(lambda) - 1 + omega))
    G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    G0 <- ((z * (exp(lambda) - lambda)) / (exp(lambda) * omega + lambda * (1 - omega)) - 
            1 / (omega + exp(lambda) - 1) - (1 - z) / (1 - omega))
    
    G0 <- G0 * omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    uMatrix <- matrix(c(G1, G0), ncol = 2)
    
    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })
    
    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      #xx <- solve(weight[[x]])
      xx %*% uMatrix[x, ]
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
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X) / 2), rep(0, NROW(X) / 2))
    }
    
    if (is.null(weight)) {
      weight <- 1
    }
    z <- as.numeric(y == 1)
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        -sum(weight * (z * (log(exp(lambda) * omega + (1 - omega) * lambda)) +
        (1 - z) * (log(1 - omega) + y * log(lambda) - lgamma(y + 1)) - 
        log(exp(lambda) - 1 + omega)))
      },
      function(beta) {
        eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- (z  * (exp(lambda) * omega + 1 - omega) / 
        (exp(lambda) * omega + lambda - omega * lambda) + 
        (1 - z) * y / lambda - exp(lambda) / (exp(lambda) - 1 + omega))
        
        G1 <- G1 * weight * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
        
        G0 <- ((z * (exp(lambda) - lambda)) / (exp(lambda) * omega + lambda * (1 - omega)) - 
               1 / (omega + exp(lambda) - 1) - (1 - z) / (1 - omega))
        
        G0 <- G0 * weight * omegaLink(eta[, 2], inverse = TRUE, deriv = 1)
        
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
        omega  <-  omegaLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        res <- matrix(nrow = length(beta), ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        temp <- exp(lambda) * omega + lambda * (1 - omega)
        temp1 <- exp(lambda) * omega + 1 - omega
        temp2 <- exp(lambda) + 1 - omega
        
        G1 <- (z  * (exp(lambda) * omega + 1 - omega) / 
        (exp(lambda) * omega + lambda - omega * lambda) + 
        (1 - z) * y / lambda - exp(lambda) / (exp(lambda) - 1 + omega))
        
        G0 <- ((z * (exp(lambda) - lambda)) / (exp(lambda) * omega + lambda * (1 - omega)) - 
                1 / (omega + exp(lambda) - 1) - (1 - z) / (1 - omega))
        
        # omega^2 derivative
        G00 <- (-(z * (exp(lambda) - lambda) ^ 2) /
        (exp(lambda) * omega + lambda * (1 - omega)) ^2 +
        1 / (omega + exp(lambda) - 1) ^ 2 - (1 - z) / (1 - omega) ^ 2)
        
        # mixed derivative
        G01 <- (z * (exp(lambda) - 1)) / (omega * exp(lambda) + (1 - omega) * lambda) -
        (z * (exp(lambda) - lambda) * (omega * exp(lambda) - omega + 1)) /
        (omega * exp(lambda) + (1 - omega) * lambda) ^ 2 + 
        exp(lambda) / (exp(lambda) + omega - 1) ^ 2
        
        # Beta^2 derivative
        G11 <- (exp(2 * lambda) / (exp(lambda) + omega - 1) ^ 2 +
        (omega * z * exp(lambda)) / (omega * exp(lambda) + (1 - omega) * lambda) -
        (z * (omega * exp(lambda) - omega + 1) ^ 2) / 
        (omega * exp(lambda) + (1 - omega) * lambda) ^ 2 -
        exp(lambda) / (exp(lambda) + omega - 1) - (y * (1 - z)) / lambda ^ 2)
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- (
          t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)] * 
           (G00 * (omegaLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2) + 
            G0 * omegaLink(eta[, 2], inverse = TRUE, deriv = 2)) * weight)) %*% 
            as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        )
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- (
          t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] *
            (G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 +
             G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)) * weight)) %*% 
            X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        )
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(
          t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * G01 * 
            omegaLink(eta[, 2], inverse = TRUE, deriv = 1) * 
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * weight) %*% 
            as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        )
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- (
          t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) * G01 * 
            omegaLink(eta[, 2], inverse = TRUE, deriv = 1) * 
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * weight) %*% 
            as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        )
        
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
    mu1 <- mu.eta(eta = eta)
    idealOmega <- ifelse(y == 1, 1, 0)
    
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
      -(log(exp(lambda) * omega + lambda * (1 - omega)) - log(exp(lambda) - 1 + omega)),
      y * log(idealLambda) - log(exp(idealLambda) - 1) - log(1 - omega) - y * log(lambda) + log(exp(lambda) - 1 + omega)
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
    
    ## see comments in ztpoisson for explanation of pmax
    sign(y - mu1) * sqrt(2 * wt * pmax(diff, 0))
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    N <- pw / (1 - (1 - omega) * exp(-lambda))
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    bigTheta1 <- -pw * omegaLink(eta[, 2], inverse = TRUE,deriv = 1) * 
      exp(lambda)/(omega + exp(lambda) - 1) ^ 2# w.r to omega
    bigTheta2 <- pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
      ((omega - 1) * exp(lambda)) / (exp(lambda) + omega - 1) ^ 2 # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    f2 <- sum(pw * (1 - omega) * exp(-lambda) / ((1 - (1 - omega) * exp(-lambda)) ^ 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        (omega * as.numeric(x == 1) +
        (1 - omega) * stats::dpois(x = x, lambda = lambda)) / 
        (1 - (1 - omega) * stats::dpois(x = 0, lambda = lambda))
      },
      "nontrunc" = {
        (1 - omega) * stats::dpois(x = x, lambda = lambda) +
        omega * as.numeric(x == 1)
      }
    )
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    omega  <-  omegaLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    CDF <- function(x) {
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, (1 - omega) * exp(-lambda), 
      omega +  (1 - omega) * stats::ppois(x, lambda))))
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
        (sizeObserved * priorWeights * (observed == 1) + .5) / (sizeObserved * priorWeights * sum(observed == 1) + 1)
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
      if ("(Intercept):omega" %in% colnames(Xvlm)) {
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
      family    = "oiztpoisson",
      etaNames  = c("lambda", "omega"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

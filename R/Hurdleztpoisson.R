#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @export
Hurdleztpoisson <- function(lambdaLink = c("log", "neglog"), 
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
        "nontrunc" = PI + (1 - PI) * exp(-lambda) * (lambda * exp(lambda) - lambda) / (1 - lambda * exp(-lambda)),
        "trunc" = PI * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) + 
        (1 - PI) * lambda * (exp(lambda) - 1) * exp(-lambda) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
      )
    } else {
      switch (type,
        "nontrunc" = {
          matrix(c(
            (1 - PI) * exp(lambda) * (exp(lambda) - lambda ^ 2 + lambda - 1) /
            (exp(lambda) - lambda) ^ 2,
            (1 - lambda) * exp(lambda) / (exp(lambda) - lambda)
          ) * c(
            lambdaLink(eta[, 1], inverse = TRUE, deriv = 1),
                piLink(eta[, 2], inverse = TRUE, deriv = 1)
          ), ncol = 2)
        },
        "trunc" = {
          matrix(c(
            (exp(2 * lambda) + (-lambda ^ 2 + PI * lambda - 2) * exp(lambda) + 1) *
            (1 - PI) / (exp(lambda) - lambda + PI - 1) ^ 2,
            -(exp(lambda) - lambda) * ((lambda - 1) * exp(lambda) + 1) / 
            (PI + exp(lambda) - lambda - 1) ^ 2
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
    "nontrunc" = PI + (1 - PI) * exp(-lambda) * lambda * 
      (exp(lambda) * (1 + lambda) - 1) / (1 - lambda * exp(-lambda)),
    "trunc" = PI * (1 - lambda * exp(-lambda)) / 
      (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) + 
      (1 - PI) * exp(-lambda) * lambda * (exp(lambda) * (1 + lambda) - 1) / 
      (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    ) - (mu.eta(eta = eta, type = type) ^ 2)
  }
  
  Wfun <- function(prior, eta, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    z <- PI * (1 - lambda * exp(-lambda)) / (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    
    YY <- mu.eta(eta) - z ## expected for (1-z)Y
    
    # PI^2 derivative
    G00 <- exp(-2 * lambda) / (-exp(-lambda) * (1 - PI) - lambda * exp(-lambda) + 1) ^ 2 - 
      z / PI ^ 2 - (1 - z) / (1 - PI) ^ 2
    
    G00 <- G00 * piLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2
    
    # mixed
    G01 <- (exp(lambda) - 1) / (exp(lambda) - lambda + PI - 1) ^ 2
    
    G01 <- G01 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
      piLink(eta[, 2], inverse = TRUE, deriv = 1)
    
    # Beta^2 derivative
    G11 <- (lambda * exp(-lambda) + (1 - PI) * exp(-lambda) - exp(-lambda)) ^ 2 / 
      (-lambda * exp(-lambda) - (1 - PI) * exp(-lambda) + 1) ^ 2 - 
      (z * (lambda * exp(-lambda) - exp(-lambda)) ^ 2) / 
      (1 - lambda * exp(-lambda)) ^ 2 + 
      (lambda * exp(-lambda) + (1 - PI) * exp(-lambda) - 2 * exp(-lambda)) /
      (-lambda * exp(-lambda) - (1 - PI) * exp(-lambda) + 1) + 
      (z * (2 * exp(-lambda) - lambda * exp(-lambda))) /
      (1 - lambda * exp(-lambda)) - YY / lambda ^ 2
    
    G11 <- G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
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
  
  funcZ <- function(eta, weight, y, prior, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    weight <- weight / prior
    
    z <- ifelse(y == 1, y, 0)
    
    G1 <- -(lambda * exp(-lambda) + (1 - PI) * exp(-lambda) - exp(-lambda)) /
      (-lambda * exp(-lambda) - (1 - PI) * exp(-lambda) + 1) +
      (z * (lambda * exp(-lambda) - exp(-lambda))) /
      (1 - lambda * exp(-lambda)) + (1 - z) * (y / lambda - 1)
    G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
    
    G0 <- z / PI - (1 - z) / (1 - PI) - 
      exp(-lambda) / (-exp(-lambda) * (1-PI) - lambda * exp(-lambda) + 1)
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
        eta    <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        PI     <- piLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)

        -sum((z * (log(1 - lambda * exp(-lambda)) + log(PI)) +
        (1 - z) * (log(1 - PI) + y * log(lambda) - lambda - lgamma(y + 1)) -
        log(1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))) * weight)
      },
      function(beta) {
        eta    <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        PI     <- piLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- -(lambda * exp(-lambda) + (1 - PI) * exp(-lambda) - exp(-lambda)) /
          (-lambda * exp(-lambda) - (1 - PI) * exp(-lambda) + 1) +
          (z * (lambda * exp(-lambda) - exp(-lambda))) /
          (1 - lambda * exp(-lambda)) + (1 - z) * (y / lambda - 1)
        G1 <- G1 * weight * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
        
        G0 <- z / PI - (1 - z) / (1 - PI) - 
          exp(-lambda) / (-exp(-lambda) * (1-PI) - lambda * exp(-lambda) + 1)
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
        eta    <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
        PI     <- piLink(eta[, 2], inverse = TRUE)
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        res <- matrix(nrow = length(beta), ncol = length(beta), 
                      dimnames = list(names(beta), names(beta)))
        
        G1 <- -(lambda * exp(-lambda) + (1 - PI) * exp(-lambda) - exp(-lambda)) /
          (-lambda * exp(-lambda) - (1 - PI) * exp(-lambda) + 1) +
          (z * (lambda * exp(-lambda) - exp(-lambda))) /
          (1 - lambda * exp(-lambda)) + (1 - z) * (y / lambda - 1)
        
        G0 <- z / PI - (1 - z) / (1 - PI) - 
          exp(-lambda) / (-exp(-lambda) * (1-PI) - lambda * exp(-lambda) + 1)
        
        # PI^2 derivative
        G00 <- exp(-2 * lambda) / (-exp(-lambda) * (1 - PI) - lambda * exp(-lambda) + 1) ^ 2 - 
          z / PI ^ 2 - (1 - z) / (1 - PI) ^ 2
        
        G00 <- G00 * piLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2 +
               G0  * piLink(eta[, 2], inverse = TRUE, deriv = 2)
        
        # mixed
        
        G01 <- (exp(lambda) - 1) / (exp(lambda) - lambda + PI - 1) ^ 2
        G01 <- G01 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) *
          piLink(eta[, 2], inverse = TRUE, deriv = 1)
        
        # Beta^2 derivative
        G11 <- (lambda * exp(-lambda) + (1 - PI) * exp(-lambda) - exp(-lambda)) ^ 2 / 
          (-lambda * exp(-lambda) - (1 - PI) * exp(-lambda) + 1) ^ 2 - 
          (z * (lambda * exp(-lambda) - exp(-lambda)) ^ 2) / 
          (1 - lambda * exp(-lambda)) ^ 2 + 
          (lambda * exp(-lambda) + (1 - PI) * exp(-lambda) - 2 * exp(-lambda)) /
          (-lambda * exp(-lambda) - (1 - PI) * exp(-lambda) + 1) + 
          (z * (2 * exp(-lambda) - lambda * exp(-lambda))) /
          (1 - lambda * exp(-lambda)) - (y * (1 - z)) / lambda ^ 2
        G11 <- G11 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2 +
               G1  * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
        
        res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- 
          t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)] * G00 * weight)) %*% 
          as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        res[1:lambdaPredNumber, 1:lambdaPredNumber] <- 
          t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * G11 * weight)) %*% 
          X[1:(nrow(X) / 2), 1:lambdaPredNumber]
        
        res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- 
          t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * G01 * weight)) %*% 
          as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
        
        res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- 
          t(t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] * G01 * weight)) %*% 
              as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)]))
        
        res
      }
    )
  }
  
  validmu <- function(mu) {
    all(0 < mu & is.finite(mu))
  }
  
  devResids <- function(y, eta, wt, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    # when pi = 0 distribution collapses to zotpoisson
    inverseFunction <- function(y) {stats::uniroot(
      f = function(x) {
        lambda <- exp(x)
        (lambda - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda)) - y
      }, 
      lower = -log(y), upper = y * 10, 
      tol = .Machine$double.eps
    )$root}
    
    yUnq <- unique(y)
    idealLambda <- tryCatch(
      expr = {
        suppressWarnings(sapply(yUnq, 
          FUN = function(x) ifelse(x %in% c(1, 2), -Inf, inverseFunction(x))
        ))
      },
      error = function (e) {
        warning("Deviance residuals could not have been computed and zero vector will be returned instead.", call. = FALSE)
        NULL
      }
    )
    if (is.null(idealLambda)) {
      return(rep(0, length(y)))
    }
    
    idealLambda <- sapply(y, FUN = function(x) idealLambda[yUnq == x])
    idealLambda <- exp(idealLambda)
    
    diff <- ifelse(
      y == 1,
      -(log(PI) + log(1 - lambda * exp(-lambda)) - log(1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))),
      ifelse(y == 2, 0,
      y * log(idealLambda) - idealLambda - log(1 - exp(-idealLambda) - idealLambda * exp(-idealLambda)) - lgamma(y + 1)) - 
      (log(1 - PI) + y * log(lambda) - lambda - lgamma(y + 1) - log(1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)))
    )
    
    if (any(diff < 0)) {
      warning(paste0(
        "Some of differences between log likelihood in sautrated model",
        " and fitted model were positive which indicates either:\n",
        "(1): A very good model fitt or\n",
        "(2): Incorrect computation of saturated model",
        "\nDouble check deviance before proceeding"
      ))
      diff[diff < 0] <- 0
    }
    
    sign(y - mu.eta(eta = eta)) * sqrt(2 * wt * diff)
  }
  
  pointEst <- function (pw, eta, contr = FALSE, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    N <- pw * (1 - lambda * exp(-lambda)) / 
      (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda))
    
    if(!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function (pw, eta, cov, Xvlm, ...) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    prob <- 1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)
    
    bigTheta1 <- -pw * piLink(eta[, 2], inverse = TRUE, deriv = 1) * 
      ((exp(lambda) - lambda) / (PI + exp(lambda) - lambda - 1) ^ 2)# w.r to PI
    bigTheta2 <- pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
      (PI - 1) * (exp(lambda) - 1) / (exp(lambda) - lambda + PI - 1) ^ 2 # w.r to lambda
    
    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)
    
    f1 <-  t(bigTheta) %*% as.matrix(cov) %*% bigTheta
    
    #f2 <- sum(pw * (1 - PI) * exp(-lambda) * (1 - lambda * exp(-lambda)) / (prob ^ 2))
    f2 <- sum(pw * (1 - lambda * exp(-lambda)) * (1 - PI) * exp(-lambda) / (prob ^ 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        ifelse(x == 1, PI * (1 - lambda * exp(-lambda)) / 
          (1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)),
          (1 - PI) * (lambda ^ x) * exp(-lambda) / 
          ((1 - (1 - PI) * exp(-lambda) - lambda * exp(-lambda)) * factorial(x))
        )
      },
      "nontrunc" = {
        ifelse(x == 1, PI, (1 - PI) *
        stats::dpois(x, lambda) / (1 - lambda * exp(-lambda)))
      }
    )
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    PI     <- piLink(eta[, 2], inverse = TRUE)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    CDF <- function(x) {
      ifelse(x == Inf, 1, 
      ifelse(x < 0, 0, 
      ifelse(x < 1, (1 - PI) * exp(-lambda) / (1 - lambda * exp(-lambda)), 
      PI + (1 - PI) * (stats::ppois(x, lambda) - lambda * exp(-lambda)
      ) / (1 - lambda * exp(-lambda)))))
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
      family    = "Hurdleztpoisson",
      etaNames = c("lambda", "pi"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

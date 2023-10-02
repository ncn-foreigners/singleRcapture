#' @rdname singleRmodels
#' @export
zotpoisson <- function(lambdaLink = c("log", "neglog"),
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
    lambda <- lambdaLink(eta, inverse = TRUE)
    
    if (!deriv) {
      switch (type,
        "nontrunc" = lambda,
        "trunc" = (lambda - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda))
      )
    } else {
      switch (type,
        "nontrunc" = lambdaLink(eta, inverse = TRUE, deriv = 1),
        "trunc" = lambdaLink(eta, inverse = TRUE, deriv = 1) *
          (exp(2 * lambda) + (-lambda ^ 2 - 2) * exp(lambda) + 1) /
          (exp(lambda) - lambda - 1) ^ 2
      )
    }
  }

  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = (lambda - lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda))
    )
  }
  
  Wfun <- function(prior, eta, y, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
    Ey <- mu.eta(eta)[iddx]
    
    G11 <- iddx
    # This looks much scaries than it ought to :) 
    # ... to bad!
    G11[iddx] <- (-((lambda - Ey) * exp(lambda) + exp(lambda) + Ey - 1) /
    (lambda * (exp(lambda) - lambda - 1)) + 
    ((lambda - Ey) * exp(lambda) + (Ey - 1) * lambda + Ey) /
    (lambda ^ 2 * (exp(lambda) - lambda - 1)) +
    ((exp(lambda) - 1) * ((lambda - Ey) * exp(lambda) + (Ey - 1) * lambda + Ey)) /
    (lambda * (exp(lambda) - lambda - 1) ^ 2)) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] ^ 2
    
    matrix(-prior * G11, ncol = 1, dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
    
    res <- iddx
    res[iddx] <- (-lambda/(-lambda -1 + exp(lambda)) + y[iddx] / lambda - 1) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] * prior[iddx] / weight[iddx, ]
    res
  }

  minusLogLike <- function(y, X, 
                           weight    = 1, 
                           NbyK      = FALSE, 
                           vectorDer = FALSE, 
                           deriv     = 0,
                           offset, ...) {
    y <- as.numeric(y)
    iddx <- y > 1
    if (is.null(weight)) {
      weight <- 1
    }
    
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X)))
    }
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- as.matrix(X) %*% beta + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        -sum(weight * iddx * (y * log(lambda) - lambda - lgamma(y + 1) -
        log(1 - exp(-lambda) - lambda * exp(-lambda))))
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
        G1 <- iddx
        
        G1[iddx] <- (-lambda/(-lambda -1 + exp(lambda)) + y[iddx] / lambda - 1) * 
          weight[iddx] * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx]
        if (NbyK) {
          return(as.data.frame(X) * G1)
        }
        if (vectorDer) {
          return(matrix(G1, ncol = 1))
        }
        
        t(as.matrix(X)) %*% G1
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
        
        G1 <- iddx
        
        G1[iddx] <- (-lambda/(-lambda - 1 + exp(lambda)) + y[iddx] / lambda - 1) * weight[iddx] *
               lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)[iddx]
        
        G11 <- iddx
        G11[iddx] <- (-((lambda - y[iddx]) * exp(lambda)+exp(lambda) + y[iddx] - 1) /
        (lambda * (exp(lambda) - lambda - 1)) + 
         ((lambda - y[iddx]) * exp(lambda) + (y[iddx] - 1) * lambda + y[iddx]) /
        (lambda ^ 2 * (exp(lambda) - lambda - 1)) +
        ((exp(lambda) - 1) * ((lambda - y[iddx]) * exp(lambda) + (y[iddx] - 1) * lambda + y[iddx])) /
        (lambda * (exp(lambda) - lambda - 1) ^ 2)) * weight[iddx] *
        lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] ^ 2
              
        
        t(X) %*% as.matrix(t(t(as.data.frame(X) * (G1 + G11))))
      }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  devResids <- function(y, eta, wt, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    iddx <- y > 1
    #print(table(y))
    
    inverseFunction <- function(y) {stats::uniroot(
      f = function(x) {mu.eta(x) - y}, 
      lower = -log(y), upper = y * 10, 
      tol = .Machine$double.eps
    )$root}
    
    # only compute predictors in saturated model for unique values of y
    # this is faster because stats::uniroot in slow whereas lamW::lambertW0 is really fast
    # so this is not worth it in ztpoisson and yes this is what I have to do because R does 
    # not have dictionaries :( Also I checked it with rbenchmark::benchmark with many replications
    yUnq <- unique(y[iddx])
    lambdaSat <- sapply(yUnq, FUN = function(x) ifelse(x == 2, -Inf, inverseFunction(x)))
    
    idealLambda <- tryCatch(
      expr = {
        suppressWarnings(sapply(yUnq, 
          FUN = function(x) ifelse(x == 2, -Inf, inverseFunction(x))
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
    
    lambdaSat <- lambdaLink(sapply(y[iddx], FUN = function(x) lambdaSat[yUnq == x]), inverse = TRUE)
    
    diff <- iddx
    lambda <- lambda[iddx]
    diff[iddx] <- y[iddx] * log(lambda) - lambda - log(1 - exp(-lambda) - lambda * exp(-lambda)) -
    ifelse(y[iddx] == 2, lgamma(y[iddx] + 1), y[iddx] * log(lambdaSat) - lambdaSat - 
    log(1 - exp(-lambdaSat) - lambdaSat * exp(-lambdaSat)))
    
    if (any(diff > 0)) {
      warning(paste0(
        "Some of differences between log likelihood in sautrated model",
        " and fitted model were positive which indicates either:\n",
        "(1): A very good model fitt or\n",
        "(2): Incorrect computation of saturated model",
        "\nDouble check deviance before proceeding"
      ))
    }

    sign(y - mu.eta(eta = eta)) * sqrt(-2 * wt * diff)
  }

  pointEst <- function (pw, eta, contr = FALSE, y, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    N <- pw * (iddx * (1 - lambda * exp(-lambda)) / 
    (1 - exp(-lambda) - lambda * exp(-lambda)) + (1 - iddx) * 1)
    
    if(!contr) {
      N <- sum(N)
    }
    
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, y, ...) {
    iddx <- y > 1
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)[iddx]
    Xvlm <- as.data.frame(Xvlm)
    
    #prob <- (1 - exp(-lambda) - lambda * exp(-lambda))
    #term <- (1 - lambda * exp(-lambda)) ^ 2
    
    f1 <- -t(Xvlm[iddx,, drop = FALSE]) %*% 
      (as.numeric(pw[iddx] * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)[iddx] * 
                                     (exp(lambda) - 1) / (exp(lambda) - lambda - 1) ^ 2))
    
    f1 <- t(f1) %*% as.matrix(cov) %*% f1
    
    f2 <- sum((pw[iddx] * (1 - lambda * exp(-lambda)) ^ 2 * 
    (exp(-lambda) + lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda)) ^ 2))
    
    f1 + f2
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    switch (type,
      "trunc" = {
        stats::dpois(x = x, lambda = lambda) / 
        (1 - stats::dpois(x = 0, lambda = lambda) - 
        stats::dpois(x = 1, lambda = lambda))
      },
      "nontrunc" = stats::dpois(x = x, lambda = lambda)
    )
  }

  simulate <- function(n, eta, lower = 1, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    lb <- stats::ppois(lower, lambda)
    ub <- stats::ppois(upper, lambda)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qpois(p_u, lambda)
    sims
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
      family    = "zotpoisson",
      etaNames  = c("lambda"),
      simulate  = simulate,
      getStart  = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

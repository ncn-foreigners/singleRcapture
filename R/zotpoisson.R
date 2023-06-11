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
  
  Wfun <- function(prior, eta, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    Ey <- mu.eta(eta)
    
    G1 <- (-lambda / (-lambda - 1 + exp(lambda)) + Ey / lambda - 1) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
    
    G11 <- (-((lambda - Ey) * exp(lambda)+exp(lambda) + Ey - 1) /
    (lambda * (exp(lambda) - lambda - 1)) + 
    ((lambda - Ey) * exp(lambda) + (Ey - 1) * lambda + Ey) /
    (lambda ^ 2 * (exp(lambda) - lambda - 1)) +
    ((exp(lambda) - 1) * ((lambda - Ey) * exp(lambda) + (Ey - 1) * lambda + Ey)) /
    (lambda * (exp(lambda) - lambda - 1) ^ 2)) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
    matrix(-prior * (G11 + G1), 
    ncol = 1, dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    (-lambda/(-lambda -1 + exp(lambda)) + y / lambda - 1) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / weight
  }

  minusLogLike <- function(y, X, weight = 1, NbyK = FALSE, vectorDer = FALSE, deriv = 0, ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        -sum(weight * (y * log(lambda) - lambda - lgamma(y + 1) -
        log(1 - exp(-lambda) - lambda * exp(-lambda))))
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        G1 <- (-lambda/(-lambda -1 + exp(lambda)) + y / lambda - 1) * weight *
               lambdaLink(eta[, 1], inverse = TRUE, deriv = 1)
        if (NbyK) {
          return(as.data.frame(X) * G1)
        }
        if (vectorDer) {
          return(matrix(G1, ncol = 1))
        }
        
        t(as.matrix(X)) %*% G1
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        
        G1 <- (-lambda/(-lambda - 1 + exp(lambda)) + y / lambda - 1) * weight *
               lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
        
        G11 <- (-((lambda - y) * exp(lambda)+exp(lambda) + y - 1) /
        (lambda * (exp(lambda) - lambda - 1)) + 
         ((lambda - y) * exp(lambda) + (y - 1) * lambda + y) /
        (lambda ^ 2 * (exp(lambda) - lambda - 1)) +
        ((exp(lambda) - 1) * ((lambda - y) * exp(lambda) + (y - 1) * lambda + y)) /
        (lambda * (exp(lambda) - lambda - 1) ^ 2)) * weight *
        lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
              
        
        t(X) %*% as.matrix(t(t(as.data.frame(X) * (G1 + G11))))
      }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  devResids <- function(y, eta, wt, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    inverseFunction <- function(y) {stats::uniroot(
      f = function(x) {mu.eta(x) - y}, 
      lower = -log(y), upper = y * 10, 
      tol = .Machine$double.eps
    )$root}
    
    # only compute predictors in saturated model for unique values of y
    # this is faster because stats::uniroot in slow whereas lamW::lambertW0 is really fast
    # so this is not worth it in ztpoisson and yes this is what I have to do because R does 
    # not have dictionaries :( Also I checked it with rbenchmark::benchmark with many replications
    yUnq <- unique(y)
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
    
    lambdaSat <- lambdaLink(sapply(y, FUN = function(x) lambdaSat[yUnq == x]), inverse = TRUE)
    
    diff <- y * log(lambda) - lambda - log(1 - exp(-lambda) - lambda * exp(-lambda)) -
    ifelse(y == 2, log(2), # log(2) is the limit as lambda->0^+
    y * log(lambdaSat) - lambdaSat - 
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
    
    ## see comments in ztpoisson for explanation of pmin
    sign(y - mu.eta(eta = eta)) * sqrt(-2 * wt * pmin(0, diff))
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    N <- pw * (1 - lambda * exp(-lambda)) / 
    (1 - exp(-lambda) - lambda * exp(-lambda))
    
    if(!contr) {
      N <- sum(N)
    }
    
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    Xvlm <- as.data.frame(Xvlm)
    
    prob <- (1 - exp(-lambda) - lambda * exp(-lambda))
    term <- (1 - lambda * exp(-lambda)) ^ 2
    
    f1 <- -t(Xvlm) %*% (as.numeric(pw * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * 
                                     (exp(lambda) - 1) / (exp(lambda) - lambda - 1) ^ 2))
    
    f1 <- t(f1) %*% as.matrix(cov) %*% f1
    
    f2 <- sum(pw * (1 - lambda * exp(-lambda)) ^ 2 * 
    (exp(-lambda) + lambda * exp(-lambda)) / (1 - exp(-lambda) - lambda * exp(-lambda)) ^ 2)
    
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
    start <- stats::glm.fit(
      x = variables[wch$reg, ],
      y = observed[wch$reg],
      family = stats::poisson(),
      weights = priorWeights[wch$reg],
      ...
    )$coefficients,
    if (attr(family$links, "linkNames")[1] == "neglog") start <- -start
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

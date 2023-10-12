#' @rdname singleRmodels
#' @importFrom lamW lambertW0
#' @importFrom stats weighted.mean
#' @export
ztpoisson <- function(lambdaLink = c("log", "neglog"),
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
        "trunc" = lambda / (1 - exp(-lambda))
      )
    } else {
      switch (type,
        "nontrunc" = lambdaLink(eta, inverse = TRUE, deriv = 1),
        "trunc" = lambdaLink(eta, inverse = TRUE, deriv = 1) *
          (exp(lambda) * (exp(lambda) - lambda - 1)) / (exp(lambda) - 1) ^ 2
      )
    }
  }

  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    switch (type,
    "nontrunc" = lambda,
    "trunc" = mu.eta(eta = eta) * (1 + lambda - mu.eta(eta = eta))
    )
  }
  
  Wfun <- function(prior, eta, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    Ey <- mu.eta(eta)
    
    G11 <- (exp(2 * lambda) / (exp(lambda) - 1) ^ 2 -
    exp(lambda) / (exp(lambda) - 1) - Ey / lambda ^ 2) *
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
    
    matrix(-G11 * prior, ncol = 1, 
           dimnames = list(rownames(eta), c("lambda")))
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
    prior * ((y / lambda - 1 / (1 - exp(-lambda))) * 
    lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / weight)
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
      offset <- cbind(rep(0, NROW(X)))
    }
    
    if (!(deriv %in% c(0, 1, 2))) stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1 # to make it conform to how switch in R works, i.e. indexing begins with 1
    
    switch (deriv,
      function(beta) {
        lambda <- lambdaLink((as.matrix(X) %*% beta + offset)[, 1], 
                             inverse = TRUE)
        
        -sum(weight * (y * log(lambda) - log(exp(lambda) - 1) - lgamma(y + 1)))
      },
      function(beta) {
        eta <- as.matrix(X) %*% beta + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- (y / lambda - 1 / (1-exp(-lambda))) * weight * 
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
        eta <- as.matrix(X) %*% beta + offset
        lambda <- lambdaLink(eta[, 1], inverse = TRUE)
        
        G1 <- (y / lambda - 1 / (1 - exp(-lambda))) * 
               lambdaLink(eta[, 1], inverse = TRUE, deriv = 2)
        
        G11 <- (exp(2*lambda) / (exp(lambda) - 1) ^ 2 -
                exp(lambda) / (exp(lambda) - 1) - y / lambda ^ 2) *
                lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2
        
        t(as.matrix(X) * weight * (G1 + G11)) %*% as.matrix(X)
      }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu)
  }

  devResids <- function(y, eta, wt, ...) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
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
    
    
    logL <- y * log(lambda) - lambda - log(1 - exp(-lambda))
    logSat <- y * ifelse(y > 1, log(idealLambda), 0) - idealLambda - ifelse(y > 1, log(1 - exp(-idealLambda)), 0)
    
    diff <- (y * log(lambda) - lambda - log(1 - exp(-lambda)) - 
    y * ifelse(y > 1, log(idealLambda), 0) - 
    idealLambda - ifelse(y > 1, log(1 - exp(-idealLambda)), 0)
    )
    
    if (any(diff > 0)) {
      warning(paste0(
        "Some of differences between log likelihood in sautrated model",
        " and fitted model were positive which indicates either:\n",
        "(1): A very good model fitt or\n",
        "(2): Incorrect computation of saturated model",
        "\nDouble check deviance before proceeding"
      ))
    }
    
    ### Here we take pairwise minimum because in specific situations
    ### lambda and idealLambda are so close for some units that
    ### their respective likelihoods differ only by machine epsilon
    ### and rounding may cause warnings
    
    ### TLDR:: pmin must be here not because mathematical error, 
    ### rather because of a rounding error
    sign(y - mu.eta(eta = eta)) * sqrt(-2 * wt * pmin(0, logL - logSat))
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
    Xvlm <- as.data.frame(Xvlm)
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)

    bigTheta <- -(lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * pw * 
                  exp(lambda) / (exp(lambda) - 1) ^ 2) %*% as.matrix(Xvlm)
    bigTheta <- as.vector(bigTheta)
    
    f1 <- t(bigTheta) %*% as.matrix(cov) %*% bigTheta

    f2 <- sum(pw * exp(-lambda) / ((1 - exp(-lambda)) ^ 2))

    f1 + f2
  }
  
  simulate <- function(n, eta, lower = 0, upper = Inf) {
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    lb <- stats::ppois(lower, lambda)
    ub <- stats::ppois(upper, lambda)
    p_u <- stats::runif(n, lb, ub)
    sims <- stats::qpois(p_u, lambda)
    sims
  }
  
  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta[, 1], inverse = TRUE)
    
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
        pmin(family$links[[1]](observed), family$links[[1]](12))
      ) + offset
      #etaStart <- etaStart
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
  
  ratioFunc <- function(x) {x+1}
  
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
      family    = "ztpoisson",
      etaNames  = c("lambda"),
      simulate  = simulate,
      getStart  = getStart,
      ratioFunc = ratioFunc,
      extraInfo = c(
        mean       = "lambda",
        variance   = "lambda",
        popSizeEst = "(1 + exp(-lambda)) ^ -1",
        meanTr     = "lambda / (1 - exp(-lambda))",
        varianceTr = 
          "(lambda / (1 - exp(-lambda))) * (1 + lambda - lambda / (1 - exp(-lambda)))"
      )
    ),
    class = c("singleRfamily", "family")
  )
}

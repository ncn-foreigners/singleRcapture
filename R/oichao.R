#' @rdname singleRmodels
#' @author Cyprian Jurkowski, Piotr Chlebicki, Maciej Beręsewicz
#' @importFrom stats dbinom
#' @export
oichao <- function(lambdaLink = "logthird",
                   ...) {
  if (missing(lambdaLink)) lambdaLink <- "logthird"
  
  links <- list()
  attr(links, "linkNames") <- c(lambdaLink)
  
  lambdaLink <- switch(
    lambdaLink,
    "logthird" = singleRinternallogthirdLink,
    stop("Currently only lambdaLink = 'logthird' is supported.")
  )
  
  links[1] <- c(lambdaLink)
  
  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    
    if (!deriv) {
      switch(
        type,
        "nontrunc" = lambda,
        "trunc" = (lambda / 3) / (1 + lambda / 3)
      )
    } else {
      switch(
        type,
        "nontrunc" = cbind(lambdaLink(eta, inverse = TRUE, deriv = 1)),
        "trunc" = cbind(
          (3 / (lambda + 3) ^ 2) *
            lambdaLink(eta, inverse = TRUE, deriv = 1)
        )
      )
    }
  }
  
  variance <- function(eta, type = "nontrunc", ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    switch(
      type,
      "nontrunc" = lambda,
      "trunc" = ((lambda / 3) / (1 + lambda / 3)) * (1 / (1 + lambda / 3))
    )
  }
  
  Wfun <- function(prior, eta, y, ...) {
    iddx <- y %in% 2:3
    lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
    res <- iddx
    res[iddx] <- -prior[iddx] * (
      -3 / (lambda * (3 + lambda) ^ 2) *
        lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx] ^ 2
    )
    
    matrix(
      res,
      ncol = 1,
      dimnames = list(rownames(eta), c("lambda"))
    )
  }
  
  funcZ <- function(eta, weight, y, prior, ...) {
    iddx <- y %in% 2:3
    z <- y - 2
    lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
    
    res <- iddx
    res[iddx] <- prior[iddx] * (((z[iddx] - 1) * lambda + 3 * z[iddx]) /
      (lambda * (3 + lambda))) *
      lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx] /
      weight[iddx, ]
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
    
    if (!(deriv %in% c(0, 1, 2))) {
      stop("Only score function and derivatives up to 2 are supported.")
    }
    
    deriv <- deriv + 1
    
    switch(
      deriv,
      function(beta, eta) {
        if (missing(eta)) {
          eta <- as.matrix(X) %*% beta + offset
        }
        lambda <- lambdaLink(eta, inverse = TRUE)[iddx]
        -sum(weight[iddx] * dbinom(
          size = 1,
          prob = (lambda / 3) / (1 + lambda / 3),
          log = TRUE,
          x = z[iddx]
        ))
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
        G00[iddx] <- (
          -(9 * z[iddx] + 6 * z[iddx] * lambda +
              (z[iddx] - 1) * lambda ^ 2) /
            ((3 + lambda) ^ 2 * lambda ^ 2)
        ) * lambdaLink(eta, inverse = TRUE, deriv = 1)[iddx] ^ 2 +
          ((z[iddx] - 1) * lambda + 3 * z[iddx]) /
          (lambda * (3 + lambda)) *
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
    diff <- z * log(lambda / 3 / (1 + lambda / 3)) +
      (1 - z) * log(1 / (1 + lambda / 3))
    diff[diff > 0] <- -0
    
    ((-1) ^ (y - 1)) * sqrt((-2 * wt * diff) * (y %in% 2:3))
  }
  
  pointEst <- function(pw, eta, contr = FALSE, y, ...) {
    lambda <- lambdaLink(eta, inverse = TRUE)
    N <- ((1 + (y %in% 2:3) / ((lambda ^ 2) / 2 + (lambda ^ 3) / 6)) * pw)
    if (!contr) {
      N <- sum(N)
    }
    N
  }
  
  popVar <- function(pw, eta, cov, Xvlm, y, ...) {
    iddx <- y %in% 2:3
    lambda <- lambdaLink(eta, inverse = TRUE)
    Xvlm <- as.matrix(Xvlm)
    prob <- (lambda ^ 2) * exp(-lambda) / 2 + (lambda ^ 3) * exp(-lambda) / 6
    
    f1 <- t((-lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * pw *
      as.numeric((lambda + lambda ^ 2 / 2) /
        (((lambda ^ 2) / 2 + (lambda ^ 3) / 6) ^ 2)))[iddx] %*%
      Xvlm[iddx, , drop = FALSE])
    
    f1 <- t(f1) %*% as.matrix(cov) %*% f1
    
    f2 <- sum((pw * (1 - prob) * ((1 + exp(-lambda) / prob) ^ 2))[iddx])
    
    f1 + f2
  }
  
  simulate <- function(n, eta, lower = 1, upper = 3) {
    lambda <- as.numeric(lambdaLink(eta, inverse = TRUE))
    lambda <- rep_len(lambda, n)

    if (is.infinite(upper)) {
      return(stats::rpois(n = n, lambda = lambda))
    }

    lower <- as.integer(floor(lower)) + 1L
    upper <- as.integer(floor(upper))
    lower <- max(lower, 0L)

    if (lower > upper) {
      stop("Simulation support must contain at least one integer value.")
    }

    support <- seq.int(lower, upper)
    probs <- vapply(
      X = support,
      FUN = function(x) stats::dpois(x = x, lambda = lambda),
      FUN.VALUE = numeric(length(lambda))
    )
    probs <- probs / rowSums(probs)

    if (length(support) == 1L) {
      return(rep.int(support, n))
    }

    as.numeric(apply(
      X = probs,
      MARGIN = 1L,
      FUN = function(prob) sample(x = support, size = 1L, prob = prob)
    ))
  }
  
  dFun <- function(x, eta, type = "trunc") {
    if (missing(type)) type <- "trunc"
    lambda <- lambdaLink(eta, inverse = TRUE)
    
    switch(
      type,
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
      etaStart <- family$links[[1]](3 * (etaStart / (1 - etaStart)))
    } else if (method == "optim") {
      init <- mean(
        ((priorWeights * (observed == 3) + .5) / (priorWeights + 1))[observed %in% 2:3]
      )
      coefStart <- rep(family$links[[1]](3 * init / (1 - init)), NCOL(variables))
    }
  )
  
  structure(
    list(
      makeMinusLogLike = minusLogLike,
      densityFunction  = dFun,
      links            = links,
      mu.eta           = mu.eta,
      valideta         = function(eta) {TRUE},
      variance         = variance,
      Wfun             = Wfun,
      funcZ            = funcZ,
      devResids        = devResids,
      validmu          = validmu,
      pointEst         = pointEst,
      popVar           = popVar,
      family           = "oichao",
      etaNames         = c("lambda"),
      simulate         = simulate,
      simulateLower    = 1,
      simulateUpper    = 3,
      getStart         = getStart
    ),
    class = c("singleRfamily", "family")
  )
}

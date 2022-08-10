# These functions are only used internally in the package so there is no need for documenting them
#' @importFrom stats optim
signleRcaptureinternalIRLS <- function(dependent,
                                       family,
                                       covariates,
                                       start,
                                       disp = NULL,
                                       omegaTheta = NULL,
                                       weights = NULL,
                                       maxiter = 10000,
                                       eps = .Machine$double.eps,
                                       omgeps = .Machine$double.eps,
                                       epsdisp = 1e-5,
                                       dispGiven = FALSE,
                                       omegaThetaGiven = FALSE,
                                       silent = FALSE,
                                       trace = 0,
                                       stepsize = 1,
                                       ...) {
  converged <- FALSE
  
  mu.eta <- family$mu.eta
  validmu <- family$validmu
  variance <- family$variance
  famName <- family$family
  linkinv <- family$linkinv
  funcZ <- family$funcZ
  Wfun <- family$Wfun
  
  if (!is.null(weights)) {
    prior <- as.numeric(weights)
  } else {
    prior <- 1
  }
  
  logLike <- family$makeMinusLogLike(y = dependent, X = covariates, weight = prior)
  grad <- family$makeGradient(y = dependent, X = covariates, weight = prior)
  
  if (famName %in% c("chao", "zelterman")) {
    dependent <- dependent - 1
    # Modification to functions in IRLS to account for logit link in fitting
    linkinv <- function(p) {1 / (1 + exp(-p))}
  }
  
  iter <- 1
  beta <- start
  
  W <- prior
  temp <- c(disp, omegaTheta, beta)
  L <- -logLike(temp)
  dispPrev <- Inf
  omegaThetaPrev <- Inf
  
  while (!converged & (iter < maxiter)) {
    if (famName %in% c("ztnegbin", "zotnegbin") &&
        isFALSE(dispGiven) && (abs(disp - dispPrev) > epsdisp)) {
      dispPrev <- disp
      ll <- function(a) logLike(c(a, beta))
      gr <- function(a) -grad(c(a, beta))[1]
      disp <- stats::optim(par = disp,
                           lower = disp - 5 * abs(disp),
                           upper = disp + 5 * abs(disp),
                           fn = ll,
                           gr = gr,
                           method = "Brent",
                           control = list(reltol = epsdisp))$par
    }
    
    halfstepsizing <- FALSE
    WPrev <- W
    tempPrev <- temp
    betaPrev <- beta
    LPrev <- L
    
    eta <- covariates %*% beta
    mu <- mu.eta(eta = eta, disp = disp, theta = omegaTheta)
    if (!validmu(mu)) {
      stop("Fit error infinite values reached consider another model,
            mu is too close to zero/infinity")
    }
    
    W <- as.numeric(Wfun(mu = mu, prior = prior, disp = disp, eta = eta, theta = omegaTheta))
    Z <- eta + funcZ(mu = mu, y = dependent, eta = eta, weight = W, disp = disp, theta = omegaTheta)
    # This is equivalent to
    # A <- t(covariates) %*% W %*% covariates
    # B <- t(covariates) %*% W %*% Z
    # But much much faster and less memory heavy
    A <- t(covariates) %*% (covariates * W)
    B <- t(covariates) %*% (Z * W)
    beta <- betaPrev + stepsize * (solve(A, B) - betaPrev)
    
    #beta <- beta - solve(A, B, tol = .Machine$double.eps)
    
    temp <- c(disp, omegaTheta, beta)
    L <- -logLike(temp)
    
    if (L < LPrev) {
      halfstepsizing <- TRUE
      h <- stepsize * (betaPrev - beta)
      if (trace > 0) {
        cat(sep = "", "Iteration number ", iter, " log-likelihood = ", format(L, scientific = FALSE, digits = 12), "\n")
        if (trace > 1) {cat(sep = " ","Parameter vector =", temp, "\n")}
        cat("Taking a modified step....\n")
      }
      repeat {
        h <- h / 2
        beta <- betaPrev - h
        temp <- c(disp, omegaTheta, beta)
        L <- -logLike(temp)
        if (L > LPrev) {
          break
        }
        
        if (max(abs(h)) < .Machine$double.eps ** (1 / 8)) {
          if (L < LPrev) {
            L <- LPrev
            beta <- betaPrev
          }
          break
        }
      }
    }
    
    if (trace > 0) {cat(sep = "", "Iteration number ", iter, " log-likelihood = ", format(L, scientific = FALSE, digits = 20), "\n")}  
    if (trace > 1) {cat(sep = " ", "Parameter vector =", temp, "\n")}
    
    converged <- ((L - LPrev < eps) || (max(abs(beta - betaPrev)) < eps))
    
    if (!converged && (iter + 1 <= maxiter)) {
      iter <- iter + 1
    } else if (L < LPrev) {
      beta <- betaPrev
      L <- LPrev
      W <- WPrev
    }
    
    if (converged && halfstepsizing && !silent) {
      warning("Convergence at halfstepsize")
    }
    
  }
  
  if(iter == maxiter && !converged && !silent) {
    warning("Fitting algorithm (IRLS) has not converged")
  }
  
  list(coefficients = beta, iter = iter, weights = W, disp = disp, omegaTheta = omegaTheta)
}
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom stats quantile
#' @importFrom stats qnorm
signleRcaptureinternalpopulationEstimate <- function(y,
                                                     X,
                                                     Xvlm = NULL,
                                                     grad,
                                                     lambda,
                                                     beta,
                                                     weights = 1,
                                                     weights0 = NULL,
                                                     hessian,
                                                     family,
                                                     eta,
                                                     dispersion,
                                                     omegaTheta,
                                                     method = "analytic",
                                                     control) {
  if (method == "noEst") {return(NULL)}
  siglevel <- control$alpha
  trcount <- control$trcount
  numboot <- control$B
  sc <- qnorm(p = 1 - siglevel / 2)
  if (method == "analytic") {
    strappedStatistic <- "No bootstrap performed"
    
    N <- family$pointEst(disp = dispersion,
                         pw = if (family$family == "zelterman") weights0 else weights,
                         lambda = lambda,
                         theta = omegaTheta) + trcount
    cov <- switch(control$covType, # Change covariance here by adding more cases
      "observedInform" = solve(-hessian(beta)),
      "Fisher" = solve(crossprod(x = as.matrix(Xvlm) * as.numeric(family$Wfun(prior = if (family$family == "zelterman") weights0 else weights, eta = eta)), 
                                 y = as.matrix(Xvlm)))
    )
    # TODO : add Fisher for negbins
    
    variation <- as.numeric(family$popVar(beta = beta, 
                                          pw = if (family$family == "zelterman") weights0 else weights,
                                          lambda = lambda,
                                          disp = dispersion,
                                          theta = omegaTheta,
                                          cov = cov, 
                                          X = X,
                                          Xvlm = Xvlm))
    
    G <- exp(sc * sqrt(log(1 + variation / ((N - length(y)) ** 2))))
    confidenceInterval <- data.frame(t(data.frame(
      "Studentized" = c(lowerBound = max(N - sc * sqrt(variation), 
                                         length(y) + trcount), 
                        upperBound = N + sc * sqrt(variation)),
      "Logtransform" = c(lowerBound = max(length(y) + (N - length(y)) / G, 
                                          (length(y) + trcount)), 
                         upperBound = length(y) + (N - length(y)) * G)
    )))
  } else if (grepl("bootstrap", method, fixed = TRUE)) {
    funBoot <- switch(control$bootType,
                      "parametric" = parBoot,
                      "semiparametric" = semparBoot,
                      "nonparametric" = noparBoot)
    N <- family$pointEst(disp = dispersion,
                         theta = omegaTheta,
                         pw = if (family$family != "zelterman") {weights} else {weights0},
                         lambda = lambda) + trcount
    
    if (!is.null(dispersion)) {
      beta <- beta[-1]
    }
    
    strappedStatistic <- funBoot(family = family,
                                 y = y, 
                                 X = X,
                                 dispersion = dispersion,
                                 beta = beta,
                                 weights = list(weights, weights0),
                                 trcount = trcount,
                                 numboot = numboot,
                                 lambda = lambda,
                                 trace = control$traceBootstrapSize,
                                 method = control$fittingMethod,
                                 control.bootstrap.method = control$bootstrapFitcontrol)

    if (N < stats::quantile(strappedStatistic, .05)) {
      warning("bootstrap statistics unusually high, try higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)")
    } else if (N > stats::quantile(strappedStatistic, .95)) {
      warning("bootstrap statistics unusually low, try higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)")
    }
    if (max(strappedStatistic) > N ** 1.375) {
      warning("Outlier(s) in statistics from bootstrap sampling detected, consider higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)")
    }
    
    variation <- stats::var(strappedStatistic)
    
    if (control$confType == "percentilic") {
      confidenceInterval <- stats::quantile(strappedStatistic,
                                            c(siglevel / 2,
                                              1 - siglevel / 2))
      names(confidenceInterval) <- c("lowerBound", "upperBound")
    } else if (control$confType == "studentized") {
      G <- exp(sc * sqrt(log(1 + variation / ((N - length(y)) ** 2))))
      
      confidenceInterval <- data.frame(t(data.frame(
        "Studentized" = c(lowerBound = max(N - sc * sqrt(variation), 
                                           (length(y) + trcount)), 
                          upperBound = N + sc * sqrt(variation)),
        "Logtransform" = c(lowerBound = max(length(y) + (N - length(y)) / G, 
                                            (length(y) + trcount)), 
                           upperBound = length(y) + (N - length(y)) * G)
      )))
    } else if (control$confType == "basic") {
      confidenceInterval <- 2 * N - stats::quantile(strappedStatistic,
                                                    c(1 - siglevel / 2, 
                                                      siglevel / 2))
      names(confidenceInterval) <- c("lowerBound", "upperBound")
    }
  }
  
  list(pointEstimate = N,
       variance = variation,
       confidenceInterval = confidenceInterval,
       boot = if (isTRUE(control$keepbootStat)) strappedStatistic else NULL,
       control = control)
}
# multiparameter
signleRcaptureinternalIRLSmultipar <- function(dependent,
                                               family,
                                               covariates,
                                               Xvlm = NULL,
                                               start,
                                               disp = NULL,
                                               omegaTheta = NULL,
                                               weights = NULL,
                                               maxiter = 10000,
                                               eps = .Machine$double.eps,
                                               silent = FALSE,
                                               trace = 0,
                                               stepsize = 1,
                                               ...) {
  converged <- FALSE
  
  momentumFactor <- 0 # add to control
  mu.eta <- family$mu.eta
  validmu <- family$validmu
  variance <- family$variance
  famName <- family$family
  linkinv <- family$linkinv
  Zfun <- family$funcZ
  Wfun <- family$Wfun
  
  if (!is.null(weights)) {
    prior <- as.numeric(weights)
  } else {
    prior <- 1
  }
  
  iter <- 1
  betaPrev <- NULL
  beta <- start
  
  logLike <- family$makeMinusLogLike(y = dependent, X = Xvlm, weight = prior)
  grad <- family$makeGradient(y = dependent, X = Xvlm, weight = prior)
  
  W <- prior
  LPrev <- -Inf
  L <- -logLike(beta)
  eta <- Xvlm %*% beta
  eta <- matrix(eta, ncol = family$parNum)
  while (!converged & (iter < maxiter)) {
    halfstepsizing <- FALSE
    mu <- mu.eta(eta = eta, ...)
    if (!validmu(mu)) {
      stop("Fit error infinite values reached consider another model,
            mu is too close to zero/infinity")
    }
    
    WPrev <- W
    W <- Wfun(prior = prior, eta = eta, y = dependent)
    z <- eta + Zfun(eta = eta, weight = W, y = dependent)
    
    XbyW <- family$multiplyWeight(X = Xvlm, W = W, thick = family$parNum)

    # A <- t(Xvlm) %*% WW %*% (Xvlm)
    # B <- t(Xvlm) %*% WW %*% (as.numeric(z))
    A <- XbyW %*% Xvlm
    B <- XbyW %*% as.numeric(z)
    
    betaPrevPrev <- betaPrev
    betaPrev <- beta
    stepPrev <- step
    
    step <- solve(A,B) - betaPrev
    beta <- betaPrev + stepsize * (step + if ((is.null(betaPrevPrev) | !momentumFactor)) {0} else {if (L-LPrev < 1) momentumFactor * stepPrev else 0})

    eta <- Xvlm %*% beta
    eta <- matrix(eta, ncol = family$parNum)
    LPrev <- L
    L <- -logLike(beta)
    
    if (trace > 0) {cat(sep = "", "Iteration number ", iter, " log-likelihood = ", format(L, scientific = FALSE, digits = 20), "\n")}
    if (trace > 1) {cat(sep = " ", "Parameter vector =", beta, "\n")}
    if (trace > 2) {cat(sep = " ", "log-likelihood reduction: ", format(L - LPrev, scientific = FALSE, digits = 20), "\n")}
    
    if (isTRUE(L < LPrev)) {
      halfstepsizing <- TRUE
      h <- stepsize * (betaPrev - beta)
      if (trace > 0) {
        cat("Taking a modified step....\n")
      }
      repeat {
        h <- h / 2
        beta <- betaPrev - h
        L <- -logLike(beta)
        if (isTRUE(L > LPrev)) {
          break
        }
        
        if (max(abs(h)) < eps) {
          if (isTRUE(L < LPrev)) {
            L <- LPrev
            beta <- betaPrev
          }
          break
        }
      }
      if (trace > 0) {cat(sep = "", "Modified step:\nIteration number ", iter, " log-likelihood = ", format(L, scientific = FALSE, digits = 12), "\n")}
      if (trace > 1) {cat(sep = " ","Parameter vector =", beta, "\n")}
      if (trace > 2) {cat(sep = " ", "log-likelihood reduction: ", format(L - LPrev, scientific = FALSE, digits = 20), "\n")}
    }
    if (trace > 0) {cat(sep = "", "----\n")}
    converged <- ((L - LPrev < eps) || (max(abs(beta - betaPrev)) < eps))

    if (!converged && (iter + 1 <= maxiter)) {
      iter <- iter + 1
    } else if ((L < LPrev) & converged) {
      beta <- betaPrev
      L <- LPrev
      W <- WPrev
    }
    
    if (converged && halfstepsizing && !silent) {
      warning("Convergence at halfstepsize")
    }
    
  }
  
  if(iter == maxiter && !converged && !silent) {
    warning("Fitting algorithm (IRLS) has not converged")
  }
  
  list(coefficients = beta, iter = iter, weights = W, disp = disp, omegaTheta = omegaTheta)
}

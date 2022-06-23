# These functions are only used internally in the package so there is no need for documenting them
#' @importFrom stats optim
signleRcaptureinternalIRLS <- function(dependent,
                 family,
                 covariates,
                 start,
                 disp = NULL,
                 weights = NULL,
                 maxiter = 10000,
                 eps = .Machine$double.eps,
                 disp.given = FALSE,
                 silent = FALSE,
                 trace) {
  converged <- FALSE
  epsdisp <- 1e-5 # TODO add to controll
  
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
  
  loglike <- family$make_minusloglike(y = dependent,
                                      X = covariates,
                                      weight = prior)
  grad <- family$make_gradient(y = dependent,
                               X = covariates,
                               weight = prior)
  
  if (famName %in% c("chao", "zelterman")) {
    dependent <- dependent - 1
    # Modification to functions in IRLS to account for logit link in fitting
    linkinv <- function(p) {1 / (1 + exp(-p))}
  }
  
  iter <- 1
  beta <- start
  
  W <- prior
  temp <- c(disp, beta)
  L <- -loglike(temp)
  dispPrev <- Inf
  
  while (!converged && (iter < maxiter)) {
    if (famName %in% c("ztnegbin", "zotnegbin") &&
        isFALSE(disp.given) && (abs(disp - dispPrev) > epsdisp)) {
      dispPrev <- disp
      ll <- function(alpha) loglike(c(alpha, beta))
      gr <- function(alpha) -grad(c(alpha, beta))[1]
      disp <- stats::optim(par = disp,
                           lower = disp - 5 * abs(disp),
                           upper = disp + 5 * abs(disp),
                           fn = ll,
                           gr = gr,
                           method = "Brent",
                           control = list(reltol = .Machine$double.eps))$par
    }
    
    halfstepsizing <- FALSE
    WPrev <- W
    tempPrev <- temp
    betaPrev <- beta
    LPrev <- L
    
    eta <- covariates %*% beta
    mu <- mu.eta(eta = eta, disp)
    if (!validmu(mu)) {
      stop("Fit error infinite values reached consider another model,
            mu is too close to zero/infinity")
    }
    
    varY <- variance(mu = mu, disp)
    W <- as.numeric(Wfun(mu = mu, prior = prior, disp = disp, eta = eta))
    Z <- funcZ(mu = mu, y = dependent, eta = eta, weight = W, disp = disp)
    # This is equivalent to
    # A <- t(covariates) %*% W %*% covariates
    # B <- t(covariates) %*% W %*% Z
    # But much much faster and less memory heavy
    A <- t(covariates) %*% (covariates * W)
    B <- t(covariates) %*% (Z * W)
    beta <- solve(A, B, tol = .Machine$double.eps)
    #beta <- beta - solve(A, B, tol = .Machine$double.eps)
    
    temp <- c(disp, beta)
    L <- -loglike(temp)
    
    if (L < LPrev) {
      halfstepsizing <- TRUE
      h <- betaPrev - beta
      if (trace %in% c("logL", "beta")) {
        if (trace == "logL") cat(sep = "", "Iteration number ", iter, " log-likelihood = ", format(L, scientific = FALSE), "\n") else if (trace == "beta") {cat(sep = "", "Iteration number ", iter, " parameter vector = ", temp)}
        cat("Taking a modified step....\n")
      }
      repeat {
        h <- h / 2
        beta <- betaPrev - h
        temp <- c(disp, beta)
        L <- -loglike(temp)
        
        if (max(abs(h)) < .Machine$double.eps ** (1/4)) {
          if (L < LPrev) {
            L <- LPrev
            beta <- betaPrev
          }
          break
        }
        if (L > LPrev) {
          break
        }
      }
    }

    if (trace == "logL") {cat(sep = "", "Iteration number ", iter, " log-likelihood = ", format(L, scientific = FALSE), "\n")} else if (trace == "beta") {cat(sep = "", "Iteration number ", iter, " parameter vector = ", temp)}
    
    converged <- ((L - LPrev < LPrev * eps) || (max(abs(beta - betaPrev)) < eps))
    
    if (!converged) {
      iter <- iter + 1
    } else if (L < LPrev) {
      beta <- betaPrev
      L <- LPrev
      W <- WPrev
    }
    
    if(iter == maxiter && !converged && !silent) {
      warning("Fitting algorithm (IRLS) has not converged")
    }
    
    if (converged && halfstepsizing && !silent) {
      warning("Convergence at halfstepsize")
    }
    
  }
  
  list(coefficients = beta, iter = iter, weights = W, disp = disp)
}
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom stats quantile
#' @importFrom stats qnorm
signleRcaptureinternalpopulationEstimate <- function(y,
                               X,
                               grad,
                               lambda,
                               beta,
                               weights = 1,
                               weights0 = NULL,
                               hessian,
                               family,
                               dispersion,
                               method = "analytic",
                               control) {
  siglevel <- control$alpha
  trcount <- control$trcount
  numboot <- control$B
  sc <- qnorm(p = 1 - siglevel / 2)
  funBoot <- switch(control$bootType,
                    "parametric" = parBoot,
                    "semiparametric" = semparBoot,
                    "nonparametric" = noparBoot)
  if (method == "analytic") {
    strappedStatistic <- "No bootstrap performed"
    
    N <- family$pointEst(disp = dispersion,
                         pw = if (family$family == "zelterman") weights0 else weights,
                         lambda = lambda) + trcount
    
    variation <- as.numeric(family$popVar(beta = beta, 
                                          pw = if (family$family == "zelterman") weights0 else weights,
                                          lambda = lambda,
                                          disp = dispersion,
                                          hess = hessian, X = X))
    
    G <- exp(sc * sqrt(log(1 + variation / ((N - length(y)) ** 2))))
    confidenceInterval <- data.frame(t(data.frame(
      "Studentized" = c(lowerBound = max(N - sc * sqrt(variation), 
                                         length(y) + trcount), 
                        upperBound = N + sc * sqrt(variation)),
      "Logtransform" = c(lowerBound = length(y) + (N - length(y)) / G, 
                         upperBound = length(y) + (N - length(y)) * G)
    )))
  } else if (grepl("bootstrap", method, fixed = TRUE)) {
    N <- family$pointEst(disp = dispersion,
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
    
    if (control$confType == "percentilic") {
      variation <- stats::var(strappedStatistic)
      confidenceInterval <- stats::quantile(strappedStatistic,
                                            c(siglevel / 2,
                                              1 - siglevel / 2))
      names(confidenceInterval) <- c("lowerBound", "upperBound")
    } else {
      variation <- stats::var(strappedStatistic)
      G <- exp(sc * sqrt(log(1 + variation / ((N - length(y)) ** 2))))
      
      confidenceInterval <- data.frame(t(data.frame(
        "Studentized" = c(lowerBound = max(N - sc * sqrt(variation), 
                                           (length(y) + trcount)), 
                          upperBound = N + sc * sqrt(variation)),
        "Logtransform" = c(lowerBound = length(y) + (N - length(y)) / G, 
                           upperBound = length(y) + (N - length(y)) * G)
      )))
    }
  }
  
  list(pointEstimate = N,
       variance = variation,
       confidenceInterval = confidenceInterval,
       boot = if (isTRUE(control$keepbootStat)) strappedStatistic else NULL,
       control = control)
}

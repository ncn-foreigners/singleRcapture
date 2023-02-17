# These functions are only used internally in the package so there is no need for documenting them
# singleRcaptureinternalIRLS <- function(dependent,
#                                        family,
#                                        covariates,
#                                        start,
#                                        weights = NULL,
#                                        maxiter = 10000,
#                                        eps = .Machine$double.eps,
#                                        silent = FALSE,
#                                        trace = 0,
#                                        stepsize = 1,
#                                        momentumFactor,
#                                        momentumActivation,
#                                        check,
#                                        epsWeights,
#                                        crit,
#                                        ...) {
#   dg <- 8 # Add to control
#   converged <- FALSE
# 
#   mu.eta <- family$mu.eta
#   validmu <- family$validmu
#   variance <- family$variance
#   famName <- family$family
#   funcZ <- family$funcZ
#   Wfun <- family$Wfun
#   prior <- as.numeric(weights)
# 
#   logLike <- family$makeMinusLogLike(y = dependent, X = covariates, weight = prior)
#   grad <- family$makeMinusLogLike(y = dependent, X = covariates, weight = prior, deriv = 1)
# 
#   if (famName %in% c("chao", "zelterman")) {
#     dependent <- dependent - 1
#   }
# 
#   traceGreaterThanFourMessegeExpr <- expression(
#     if (trace > 4) {
#       cat(sep = " ", "\nAlgorithm will terminate if one of following conditions will be met:\n")
#       if ("abstol" %in% crit) {
#         cat("The increase to minus log-likelihood will be bellow chosen value of epsilon", eps, "\n")
#       }
#       if ("reltol" %in% crit) {
#         cat("The relative increase to minus log-likelihood will be bellow chosen value of epsilon", eps, "value at current step", format((L - LPrev) / LPrev, scientific = FALSE, digits = dg), "\n")
#       }
#       if ("coef" %in% crit) {
#         cat("Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.\nAt current step the highest change was:", format(max(abs(beta - betaPrev)), scientific = FALSE, digits = dg))
#       }
#     }
#   )
#   convergence <- expression(
#     any(
#       if ("abstol" %in% crit) {
#         L - LPrev < eps
#       } else FALSE,
#       if ("reltol" %in% crit) {
#         L - LPrev < eps * LPrev
#       } else FALSE,
#       if ("coef" %in% crit) {
#         max(abs(beta - betaPrev)) < eps
#       } else FALSE
#     )
#   )
# 
#   iter <- 1
#   step <- NULL
#   beta <- start
# 
#   W <- prior
#   L <- -logLike(beta)
# 
#   while (!converged & (iter < maxiter)) {
#     # if (famName %in% c("ztnegbin", "zotnegbin") &&
#     #     isFALSE(dispGiven) && (abs(disp - dispPrev) > epsdisp)) {
#     #   dispPrev <- disp
#     #   ll <- function(a) logLike(c(a, beta))
#     #   gr <- function(a) -grad(c(a, beta))[1]
#     #   disp <- stats::optim(par = disp,
#     #                        lower = disp - 5 * abs(disp),
#     #                        upper = disp + 5 * abs(disp),
#     #                        fn = ll,
#     #                        gr = gr,
#     #                        method = "Brent",
#     #                        control = list(reltol = epsdisp))$par
#     # }
# 
#     halfstepsizing <- FALSE
#     WPrev <- W
#     betaPrev <- beta
#     LPrev <- L
# 
#     eta <- covariates %*% beta
#     mu <- mu.eta(eta = eta)
#     # if (!validmu(mu)) {
#     #   stop("Fit error infinite values reached consider another model,
#     #         mu is too close to zero/infinity")
#     # }
# 
#     W <- Wfun(mu = mu, prior = prior, eta = eta)
#     if (any(!is.finite(W))) {
#       if (!silent) {
#         warning("NA's or NaN's or infinite values in weights matrixes detected IRLS may not work propperly.")
#       }
#       W[!is.finite(W)] <- epsWeights
#     }
#     if (check) {
#       W[, (1:family$parNum) ^ 2] <- ifelse(
#         W[, (1:family$parNum) ^ 2] < epsWeights,
#         epsWeights,
#         W[, (1:family$parNum) ^ 2]
#       )
#     }
#     Z <- eta + funcZ(mu = mu, y = dependent, eta = eta, weight = W)
#     if (any(is.nan(Z))) {
#       stop("Pseudo residuals could not be computed at current iteration, possibly infinite or non numeric values in weights appeared.")
#     }
#     # This is equivalent to
#     # A <- t(covariates) %*% W %*% covariates
#     # B <- t(covariates) %*% W %*% Z
#     # But much much faster and less memory heavy
#     A <- t(covariates) %*% (covariates * as.numeric(W))
#     B <- t(covariates) %*% (Z * as.numeric(W))
# 
#     stepPrev <- step
#     step <- solve(A,B) - betaPrev
# 
#     beta <- betaPrev + stepsize *
#     (step + if ((is.null(stepPrev) | !momentumFactor)) 0 else {
#     if (L-LPrev < 1) momentumFactor * stepPrev else 0})
# 
#     #beta <- beta - solve(A, B, tol = .Machine$double.eps)
# 
#     L <- -logLike(beta)
# 
#     if (trace > 0) {cat(sep = "", "Iteration number ", iter, " log-likelihood: ", format(L, scientific = FALSE, digits = dg))}
#     if (trace > 1) {cat(sep = " ", "\nParameter vector: ", format(beta, scientific = FALSE, digits = dg))}
#     if (trace > 2) {cat(sep = " ", "\nlog-likelihood reduction: ", format(L - LPrev, scientific = FALSE, digits = dg))}
#     if (trace > 3) {cat(sep = " ", "\nValue of gradient at current step:\n", format(grad(beta), scientific = FALSE, digits = dg))}
#     eval(traceGreaterThanFourMessegeExpr)
# 
#     if (isTRUE(L < LPrev) || is.infinite(L) || is.nan(L)) {
#       halfstepsizing <- TRUE
#       h <- stepsize * (betaPrev - beta)
#       if (trace > 0) {
#         cat("\nTaking a modified step....\n")
#       }
#       repeat {
#         h <- h / 2
#         beta <- betaPrev - h
#         L <- -logLike(beta)
#         if (isTRUE(L > LPrev) && is.finite(L)) {
#           break
#         }
# 
#         if (isTRUE(max(abs(h)) < .Machine$double.eps ^ (1 / 8))) {
#           if (isTRUE(L < LPrev)) {
#             if (!silent) {
#               warning("IRLS half-stepping terminated because the step is too small.")
#             }
#             halfstepsizing <- FALSE
#             L <- LPrev
#             beta <- betaPrev
#           }
#           break
#         }
#       }
#       if (trace > 0) {cat(sep = "", "Modified step:\nIteration number ", iter, " log-likelihood: ", format(L, scientific = FALSE, digits = dg))}
#       if (trace > 1) {cat(sep = " ","\nParameter vector:", format(beta, scientific = FALSE, digits = dg))}
#       if (trace > 2) {cat(sep = " ", "\nlog-likelihood reduction: ", format(L - LPrev, scientific = FALSE, digits = dg))}
#       if (trace > 3) {cat(sep = " ", "\nValue of gradient at current (modified) step:\n", format(grad(beta), scientific = FALSE, digits = dg))}
#       eval(traceGreaterThanFourMessegeExpr)
#     }
#     if (trace > 0) {cat(sep = "", "\n----\n")}
#     converged <- eval(convergence)
# 
#     if (!converged && (iter + 1 <= maxiter)) {
#       iter <- iter + 1
#     } else if (L < LPrev) {
#       beta <- betaPrev
#       L <- LPrev
#       W <- WPrev
#     }
# 
#     if (converged && halfstepsizing && !silent) {
#       warning("Convergence at halfstepsize")
#     }
# 
#   }
# 
#   if(iter == maxiter && !converged && !silent) {
#     warning("Fitting algorithm (IRLS) has not converged")
#   }
# 
#   if (!validmu(mu)) {
#     stop("Fit error infinite values reached consider another model,
#             mu is too close to zero/infinity")
#   }
# 
#   list(coefficients = beta, iter = iter, weights = W)
# }
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom stats quantile
#' @importFrom stats qnorm
singleRcaptureinternalpopulationEstimate <- function(y, X, grad,
                                                     beta, weights,
                                                     hessian, family,
                                                     eta, popVar,
                                                     control,
                                                     Xvlm, W, formulas,
                                                     sizeObserved,
                                                     modelFrame, cov) {
  #if (popVar == "noEst") {return(NULL)} moved to main function to avoid copying function parameters
  hwm <- attr(Xvlm, "hwm")
  siglevel <- control$alpha
  trcount <- control$trcount
  numboot <- control$B
  sc <- qnorm(p = 1 - siglevel / 2)
  
  if (popVar == "analytic") {
    strappedStatistic <- "No bootstrap performed"
    N <- family$pointEst(pw = weights, eta = eta) + trcount
    if (is.null(cov)) {
      cov <- switch(control$covType, # Change covariance here by adding more cases
      "observedInform" = solve(-hessian(beta)),
      "Fisher" = solve(singleRinternalMultiplyWeight(X = Xvlm, W = W) %*% Xvlm)
      )
    }
    
    variation <- as.numeric(family$popVar(
      eta = eta, 
      pw = weights,
      cov = cov,
      Xvlm = if (family$family == "zelterman") X else Xvlm
    ))
    
    sd <- sqrt(variation)
    if (control$sd == "normalMVUE") {
      sd <- sd / (sqrt(2 / (sizeObserved - 1)) * exp(lgamma(sizeObserved / 2) - lgamma((sizeObserved - 1) / 2)))
    }
    
    G <- exp(sc * sqrt(log(1 + variation / ((N - sizeObserved) ^ 2))))
    confidenceInterval <- data.frame(t(data.frame(
      "Studentized" = c(lowerBound = max(N - sc * sd, 
                                         sizeObserved), 
                        upperBound = N + sc * sd),
      "Logtransform" = c(lowerBound = max(sizeObserved + (N - sizeObserved) / G, 
                                          sizeObserved), 
                         upperBound = sizeObserved + (N - sizeObserved) * G)
    )))
  } else if (grepl("bootstrap", popVar, fixed = TRUE)) {
    funBoot <- switch(control$bootType,
                      "parametric" = parBoot,
                      "semiparametric" = semparBoot,
                      "nonparametric" = noparBoot)
    
    N <- family$pointEst(
      pw = if (family$family == "chao") weights[y %in% 1:2] else if (grepl(pattern = "^zot", x = family$family)) weights[y > 1] else weights,
      eta = eta) + trcount
    
    strappedStatistic <- funBoot(
      family = family, formulas = formulas,
      y = y, X = X, hwm = hwm,
      beta = beta, weights = weights,
      trcount = trcount, numboot = numboot,
      eta = eta, trace = control$traceBootstrapSize,
      visT = control$bootstrapVisualTrace,
      method = control$fittingMethod,
      controlBootstrapMethod = control$bootstrapFitcontrol,
      N = N, Xvlm = Xvlm, modelFrame = modelFrame
    )

    if (N < stats::quantile(strappedStatistic, .05)) {
      warning("bootstrap statistics unusually high, try higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)")
    } else if (N > stats::quantile(strappedStatistic, .95)) {
      warning("bootstrap statistics unusually low, try higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)")
    }
    if (max(strappedStatistic) > N ^ 1.5) {
      warning("Outlier(s) in statistics from bootstrap sampling detected, consider higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)")
    }
    
    variation <- stats::var(strappedStatistic)
    
    sd <- sqrt(variation)
    if (control$sd == "normalMVUE") {
      sd <- sd / (sqrt(2 / (sizeObserved - 1)) * exp(lgamma(sizeObserved / 2) - lgamma((sizeObserved- 1) / 2)))
    }
    
    if (control$confType == "percentilic") {
      confidenceInterval <- stats::quantile(strappedStatistic,
                                            c(siglevel / 2,
                                              1 - siglevel / 2))
      names(confidenceInterval) <- c("lowerBound", "upperBound")
    } else if (control$confType == "studentized") {
      G <- exp(sc * sqrt(log(1 + variation / ((N - sizeObserved) ^ 2))))
      
      confidenceInterval <- data.frame(t(data.frame(
        "Studentized" = c(lowerBound = max(N - sc * sd, 
                                           sizeObserved), 
                          upperBound = N + sc * sd),
        "Logtransform" = c(lowerBound = max(sizeObserved + (N - sizeObserved) / G, 
                                            sizeObserved), 
                           upperBound = sizeObserved + (N - sizeObserved) * G)
      )))
    } else if (control$confType == "basic") {
      confidenceInterval <- 2 * N - stats::quantile(strappedStatistic,
                                                    c(1 - siglevel / 2, 
                                                      siglevel / 2))
      names(confidenceInterval) <- c("lowerBound", "upperBound")
    }
  }
  
  structure(
    list(
      pointEstimate = N,
      variance = variation,
      confidenceInterval = confidenceInterval,
      boot = if (isTRUE(control$keepbootStat)) strappedStatistic else NULL,
      control = control
    ),
    class = "popSizeEstResults"
  )
}
# multiparameter
singleRcaptureinternalIRLSmultipar <- function(dependent,
                                               family,
                                               formulas,
                                               covariates,
                                               start,
                                               weights,
                                               maxiter = 10000,
                                               eps = .Machine$double.eps,
                                               silent = FALSE,
                                               trace = 0,
                                               stepsize = 1,
                                               momentumFactor,
                                               momentumActivation,
                                               check,
                                               epsWeights,
                                               crit,
                                               ...) {
  dg <- 8
  converged <- FALSE
  
  # Lowering stepsize to about .3 usually helps a great deal in IRLS fitting
  if (!silent && family$family == "zotnegbin" && stepsize == 1) {
    cat("Zero one truncated negative binomial distribution is prone to taking alpha parameter to infinity, consider lowering stepsize control parameter if fitting fails.")
  }
  
  mu.eta <- family$mu.eta
  validmu <- family$validmu
  variance <- family$variance
  famName <- family$family
  Zfun <- family$funcZ
  Wfun <- family$Wfun
  prior <- as.numeric(weights)
  
  iter <- 1
  step <- NULL
  betaPrev <- NULL
  beta <- start
  
  if (famName %in% c("chao", "zelterman")) {
    dependent <- dependent - 1
  }
  
  logLike <- family$makeMinusLogLike(y = dependent, X = covariates, weight = prior)
  grad <- family$makeMinusLogLike(y = dependent, X = covariates, weight = prior, deriv = 1)
  traceGreaterThanFourMessegeExpr <- expression(
    if (trace > 4) {
      cat(sep = " ", "\nAlgorithm will terminate if one of following conditions will be met:\n")
      if ("abstol" %in% crit) {
        cat("The increase to minus log-likelihood will be bellow chosen value of epsilon", eps, "\n")
      }
      if ("reltol" %in% crit) {
        cat("The relative increase to minus log-likelihood will be bellow chosen value of epsilon", eps, "value at current step", format((L - LPrev) / LPrev, scientific = FALSE, digits = dg), "\n")
      }
      if ("coef" %in% crit) {
        cat("Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.\nAt current step the highest change was:", format(max(abs(beta - betaPrev)), scientific = FALSE, digits = dg))
      }
    }
  )
  convergence <- expression(
    any(
      if ("abstol" %in% crit) {
        L - LPrev < eps
      } else FALSE,
      if ("reltol" %in% crit) {
        L - LPrev < eps * LPrev
      } else FALSE,
      if ("coef" %in% crit) {
        max(abs(beta - betaPrev)) < eps
      } else FALSE
    )
  )
  
  
  W <- prior
  LPrev <- -Inf
  L <- -logLike(beta)
  eta <- covariates %*% beta
  eta <- matrix(eta, ncol = family$parNum)
  while (!converged & (iter < maxiter)) {
    halfstepsizing <- FALSE
    mu <- mu.eta(eta = eta, ...)
    # if (!validmu(mu)) {
    #   stop("Fit error infinite values reached consider another model,
    #         mu is too close to zero/infinity")
    # }
    WPrev <- W
    W <- Wfun(prior = prior, eta = eta, y = dependent)
    if (any(!is.finite(W))) {
      if (!silent) {
        warning("NA's or NaN's or infinite values in weights matrixes detected IRLS may not work propperly.")
      }
      W[!is.finite(W)] <- epsWeights
    }
    if (check) {
      W[, (1:family$parNum) ^ 2] <- ifelse(
        W[, (1:family$parNum) ^ 2] < epsWeights, 
        epsWeights, 
        W[, (1:family$parNum) ^ 2]
      )
    }
    z <- eta + Zfun(eta = eta, weight = W, y = dependent)
    if (any(is.nan(z))) {
      stop("Pseudo residuals could not be computed at current iteration, possibly infinite or non numeric values in weights appeared.")
    }
    XbyW <- singleRinternalMultiplyWeight(X = covariates, W = W)
    # A <- t(Xvlm) %*% WW %*% (Xvlm)
    # B <- t(Xvlm) %*% WW %*% (as.numeric(z))
    A <- XbyW %*% covariates
    B <- XbyW %*% as.numeric(z)
    betaPrev <- beta
    stepPrev <- step
    step <- solve(A,B) - betaPrev
    beta <- betaPrev + stepsize * 
    (step + if ((is.null(stepPrev) | !momentumFactor)) 0 else {
    if (L-LPrev < momentumActivation) momentumFactor * stepPrev else 0
    })
    eta <- covariates %*% beta
    eta <- matrix(eta, ncol = family$parNum)
    LPrev <- L
    L <- -logLike(beta)
    
    if (trace > 0) {cat(sep = "", "Iteration number ", iter, " log-likelihood: ", format(L, scientific = FALSE, digits = dg))}
    if (trace > 1) {cat(sep = " ", "\nParameter vector: ", format(beta, scientific = FALSE, digits = dg))}
    if (trace > 2) {cat(sep = " ", "\nlog-likelihood reduction: ", format(L - LPrev, scientific = FALSE, digits = dg))}
    if (trace > 3) {cat(sep = " ", "\nValue of gradient at current step:\n", format(grad(beta), scientific = FALSE, digits = dg))}
    eval(traceGreaterThanFourMessegeExpr)

    if (isTRUE(L < LPrev) || is.infinite(L) || is.nan(L)) {
      halfstepsizing <- TRUE
      h <- step <- stepsize * (betaPrev - beta)
      if (trace > 0) {
        cat("\nTaking a modified step....\n")
      }
      repeat {
        h <- step <- h / 2
        beta <- betaPrev - h
        L <- -logLike(beta)
        if (isTRUE(L > LPrev) && is.finite(L)) {
          break
        }

        if (isTRUE(max(abs(h)) < .Machine$double.eps)) {
          if (isTRUE(L < LPrev)) {
            if (!silent) {
              warning("IRLS half-stepping terminated because the step is too small.")
            }
            halfstepsizing <- FALSE
            L <- LPrev
            beta <- betaPrev
          }
          break
        }
      }
      if (trace > 0) {cat(sep = "", "Modified step:\nIteration number ", iter, " log-likelihood: ", format(L, scientific = FALSE, digits = dg))}
      if (trace > 1) {cat(sep = " ","\nParameter vector:", format(beta, scientific = FALSE, digits = dg))}
      if (trace > 2) {cat(sep = " ", "\nlog-likelihood reduction: ", format(L - LPrev, scientific = FALSE, digits = dg))}
      if (trace > 3) {cat(sep = " ", "\nValue of gradient at current (modified) step:\n", format(grad(beta), scientific = FALSE, digits = dg))}
      eval(traceGreaterThanFourMessegeExpr)
    }
    if (trace > 0) {cat(sep = "", "\n----\n")}
    converged <- eval(convergence)

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
  
  if (!validmu(mu)) {
    stop("Fit error infinite values reached consider another model,
          mu is too close to zero/infinity")
  }
  
  list(coefficients = beta, iter = iter, weights = W)
}
# make Xvlm matrix
singleRinternalGetXvlmMatrix <- function(X, nPar, formulas, parNames) {
  formulas[[1]][[2]] <- NULL
  Xses <- list()
  parentFrame <- X
  for (k in 1:nPar) {
    Xses[[k]] <- model.matrix(formulas[[k]], data = parentFrame)
    if (k != 1) {
      colnames(Xses[[k]]) <- paste0(colnames(Xses[[k]]), ":", parNames[k])
    }
  }
  hwm <- sapply(Xses, ncol)
  # TODO: sparse matrix, maybe use Matrix
  Xvlm <- matrix(0, nrow = nPar * nrow(X), ncol = sum(hwm))
  colnames(Xvlm) <- unlist(sapply(X = Xses, FUN = colnames))
  row <- 0
  col <- 0
  for (k in Xses) {
    Xvlm[(row+1):(row + nrow(k)), (col+1):(col + ncol(k))] <- as.matrix(k)
    row <- row + nrow(k)
    col <- col + ncol(k)
  }
  attr(Xvlm, "hwm") <- hwm
  Xvlm
}
# Chosing data for estimation/regression
singleRcaptureinternalDataCleanupSpecialCases <- function (family, observed, popVar) {
  if (grepl("zot", family$family)) {
    trr <- sum(observed == 1)
    wch1 <- wch2 <- (observed > 1)
    if (popVar != "analytic") {
      # in bootstrap we need all
      wch2 <- rep(TRUE, length(observed))
    }
  } else if (family$family == "chao") {
    trr <- sum(observed > 2)
    wch1 <- wch2 <- (observed %in% c(1, 2))
    if (popVar != "analytic") {
      # in bootstrap we need all
      wch2 <- rep(TRUE, length(observed))
    }
  } else if (family$family == "zelterman") {
    # In zelterman model regression is indeed based only on 1 and 2 counts
    # but estimation is based on ALL counts
    wch1 <- (observed %in% c(1, 2))
    wch2 <- rep(TRUE, length(observed))
    trr <- 0
  } else {
    trr <- 0
    wch1 <- wch2 <- rep(TRUE, length(observed))
  }
  list(reg = wch1, # which rows for regression
       est = wch2, # which rows for estimation
       trr = trr)  # add to trcount
}
#' @importFrom stats reformulate
singleRinternalMergeFormulas <- function(ff) {
  # This code was inspired by: https://stevencarlislewalker.wordpress.com/2012/08/06/merging-combining-adding-together-two-formula-objects-in-r/
  #ff is a list with many formulas
  env <- environment(ff[[1]])
  for (k in 1:length(ff)) { # eliminate left hand side
    if (length(ff[[k]]) == 3) {
      resp <- ff[[k]][[2]]
      ff[[k]][[2]] <- NULL
    }
  }
  dotCheck <- sapply(ff, FUN = function(x) {x == ~.})
  if (any(dotCheck)) {
    out <- reformulate(".", resp)
  } else {
    for (k in 1:length(ff)) { # extract right hand side and create string vector
      ff[[k]] <- strsplit(deparse(ff[[k]][[2]]), " \\+ ")[[1]]
    }
    ff <- unlist(ff)
    out <- reformulate(ff, resp)
  }
  
  environment(out) <- env
  out
}

# This is almost certainly an overkill but is supports arbitrary number of linear predictors
singleRinternalMultiplyWeight <- function (X, W, ...) {
  hwm <- attr(X, "hwm")
  thick <- sqrt(ncol(W))
  XbyW <- matrix(0, nrow = ncol(X), ncol = nrow(X))
  wch <- c(0, cumsum(hwm))
  whichVector <- t(matrix(1:(thick ^ 2), ncol = thick))
  index1 <- 1:(nrow(X) / thick)
  for (i in 1:thick) {
    index <- 1:(nrow(X) / thick)
    for (j in 1:thick) {
      XbyW[(wch[i]+1):wch[i+1], index] <- t(X[index1, (wch[i]+1):wch[i+1]] * W[, whichVector[i, j]])
      index <- index + nrow(X) / thick
    }
    index1 <- index1 + nrow(X) / thick
  }
  XbyW
}
# TODO
cholFroW <- function(W, prior) {
  if (NROW(W) != NROW(prior)) stop("Error in estimatePopsize.fit, working weights and prior weights suggest different number of observations.")
  L <- list()
  for (k in 1:NROW(W)) {
    L[[k]] <- chol(matrix(prior[k] * W[k,], ncol = 2, nrow = 2))
  }
  L
}
# TODO
MultibyCholW <- function(Z, X, cholW, which) {
  # zmień na większą liczbe parametrów
  W <- cholW
  if (!is.null(Z)) {# 1 2 to par 3 to mieszane
    return(c(Z[1:NROW(W)]*cholW[,1]+Z[NROW(W)+1:2*NROW(W)]*cholW[,3], Z[1:NROW(W)]*cholW[,2]))
  } else if (!is.null(X)) {
    # to samo co wyżej może wykorzystaj istniejące funkcje
    # to jest explicite używając własnoći dekompozycji choleskiego
    return(cbind(
      rbind(cholW[,1]*X[1:(NROW(X)/2),1:which[1]]+cholW[,3]*X[-(1:(NROW(X)/2)),-(1:which[1])],
            cholW[,2]*X[-(1:(NROW(X)/2)),1:which[1]]),
      rbind(cholW[,1]*X[1:(NROW(X)/2),-(1:which[1])]+cholW[,3]*X[-(1:(NROW(X)/2)),-(1:which[1])],
            cholW[,2]*X[-(1:(NROW(X)/2)),-(1:which[1])])
    ))
  }
}
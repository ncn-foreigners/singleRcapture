# These functions are only used internally in the package so there is no need for documenting them
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
        "observedInform" = {
          tryCatch(
            expr = {suppressWarnings(solve(-hessian(beta)))},
            error = function (e) "e"
          )
        },
        "Fisher" = solve(singleRinternalMultiplyWeight(X = Xvlm, W = W) %*% Xvlm)
      )
      
      if (isTRUE(cov == "e")) {
        warning(
          "The analytically computed hessian of log-likelihood is numerically singular.",
          "Fisher information matrix will be used in population size estimation instead.",
          sep = " "
        )
        control$covType <- "Fisher"
        cov <- solve(singleRinternalMultiplyWeight(X = Xvlm, W = W) %*% Xvlm)
      }
    }
    
    variation <- as.numeric(family$popVar(
      eta = eta, 
      pw = weights,
      cov = cov,
      Xvlm = if (family$family == "zelterman") X else Xvlm
    ))
    
    if (!is.finite(variation))
      stop("Computed variance is infinite/NaN/NULL")
    
    sd <- sqrt(variation)
    if (control$sd == "normalMVUE") {
      sd <- sd / (sqrt(2 / (sizeObserved - 1)) * 
      exp(lgamma(sizeObserved / 2) - lgamma((sizeObserved - 1) / 2)))
    }
    
    G <- exp(sc * sqrt(log(1 + variation / ((N - sizeObserved) ^ 2))))
    confidenceInterval <- data.frame(t(data.frame(
      "normal" = c(lowerBound = max(N - sc * sd, sizeObserved), 
                   upperBound = N + sc * sd),
      "logNormal" = c(lowerBound = sizeObserved + max((N - sizeObserved) / G, 0), 
                      upperBound = sizeObserved + (N - sizeObserved) * G)
    )))
  } else if (grepl("bootstrap", popVar, fixed = TRUE)) {
    funBoot <- switch(
      control$bootType,
      "parametric" = parBoot,
      "semiparametric" = semparBoot,
      "nonparametric" = noparBoot
    )
    
    N <- family$pointEst(
      pw = if (family$family == "chao") weights[y %in% 1:2] 
      else if (grepl(pattern = "^zot", x = family$family)) weights[y > 1] 
      else weights,
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
      warning("bootstrap statistics unusually high, try higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)\n")
    } else if (N > stats::quantile(strappedStatistic, .95)) {
      warning("bootstrap statistics unusually low, try higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)\n")
    }
    if (max(strappedStatistic) > N ^ 1.5) {
      warning("Outlier(s) in statistics from bootstrap sampling detected, consider higher maxiter/lower epsilon for fitting bootstrap samples (bootstrapFitcontrol)\n")
    }
    
    variation <- stats::var(strappedStatistic)
    
    
    if (!is.finite(variation))
      stop("Computed variance is infinite/NaN/NULL")
    
    sd <- sqrt(variation)
    if (control$sd == "normalMVUE") {
      sd <- sd / (sqrt(2 / (sizeObserved - 1)) * exp(lgamma(sizeObserved / 2) - lgamma((sizeObserved- 1) / 2)))
    }
    
    if (control$confType == "percentilic") {
      confidenceInterval <- stats::quantile(strappedStatistic,
                                            c(siglevel / 2,
                                              1 - siglevel / 2))
      names(confidenceInterval) <- c("lowerBound", "upperBound")
    } else if (control$confType == "normal") {
      G <- exp(sc * sqrt(log(1 + variation / ((N - sizeObserved) ^ 2))))
      
      confidenceInterval <- data.frame(t(data.frame(
        "normal" = c(lowerBound = max(N - sc * sd, 
                                           sizeObserved), 
                          upperBound = N + sc * sd),
        "logNormal" = c(lowerBound = max(sizeObserved + (N - sizeObserved) / G, 
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
                                               printOften,
                                               saveLog,
                                               ...) {
  dg <- 8 # add to controll
  converged <- FALSE
  
  # Lowering stepsize to about .3 usually helps a great deal in IRLS fitting
  if (!silent && family$family == "zotnegbin" && stepsize == 1) {
    cat("Zero one truncated negative binomial distribution is prone to taking alpha parameter to infinity.",
        "\nConsider lowering stepsize control parameter if fitting fails.\n")
  }
  
  mu.eta   <- family$mu.eta
  validmu  <- family$validmu
  variance <- family$variance
  famName  <- family$family
  Zfun     <- family$funcZ
  Wfun     <- family$Wfun
  prior    <- as.numeric(weights)
  
  iter <- 1
  step <- NULL
  betaPrev <- NULL
  beta <- start
  
  if (famName %in% c("chao", "zelterman")) {
    dependent <- dependent - 1
  }
  
  logLike <- family$makeMinusLogLike(y = dependent, X = covariates, weight = prior)
  grad    <- family$makeMinusLogLike(y = dependent, X = covariates, weight = prior, deriv = 1)
  
  logg <- NULL
  if (isTRUE(saveLog)) {
    logg <- data.frame()
  }
  
  addToLog <- expression(
    if (isTRUE(saveLog)) {
      if (trace == 1) {logg[NROW(logg) + 1, c(1, 2, 3)] <- c(iter, halfstepsizing, L)}
      if (trace == 2) {logg[NROW(logg) + 1, c(1:3, 4:(3+length(beta)))] <- c(iter, halfstepsizing, L, beta)}
      if (trace  > 2) {logg[NROW(logg) + 1, c(1:3, 4:(3+2*length(beta)))] <- c(iter, halfstepsizing, L, beta, grad(beta))}
    }
  )
  
  traceGreaterThanFourMessegeExpr <- expression(
    if (trace > 4) {
      cat(sep = " ", "\nAlgorithm will terminate if one of following conditions will be met:\n")
      if ("abstol" %in% crit) {
        cat("The increase to minus log-likelihood will be bellow chosen value of epsilon", eps, "\n")
      }
      if ("reltol" %in% crit) {
        cat("The relative increase to minus log-likelihood will be bellow chosen value of epsilon", eps, 
            "value at current step", format((L - LPrev) / LPrev, scientific = FALSE, digits = dg), "\n")
      }
      if ("coef" %in% crit) {
        cat("Maximum change to the vector of regression parameters will be bellow the chosen value of epsilon.\nAt current step the highest change was:", 
            format(max(abs(beta - betaPrev)), scientific = FALSE, digits = dg))
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
  
  parNum <- length(family$etaNames)
  W <- prior
  LPrev <- -Inf
  L <- -logLike(beta)
  eta <- covariates %*% beta
  eta <- matrix(eta, ncol = parNum)
  while (!converged & (iter < maxiter)) {
    halfstepsizing <- FALSE
    mu <- mu.eta(eta = eta, ...)
    if (!validmu(mu)) {
      mu <- mu.eta(eta = eta, ...)
      stop(paste0(
        "Fit error infinite values reached consider another model,",
        "mu is too close to zero/infinity.\n"
      ))
    }

    if (any(!is.finite(
      sapply(1:length(family$etaNames), FUN = function(x) {
        family$links[[x]](eta[, x], inverse = TRUE)
      })
    )))
      stop(paste0(
        "Infinite values of distribution ",
        "parameters obtained for some IRLS iteration.\n"
      ))
    
    WPrev <- W
    tryCatch(
      expr = {
        W <- Wfun(prior = prior, eta = eta, y = dependent)
      },
      error = function (e) {
        stop(cat(
          "Working weight matrixes at iteration:", 
          iter, 
          "could not have been computer.", sep = " "
        ))
      }
    )
    
    if (any(!is.finite(W))) {
      if (!silent) {
        warning(paste0(
          "NA's or NaN's or infinite values in weights matrixes ",
          "detected IRLS may not work propperly.\n",
          "Replacing these values by weightsEpsilon control parameter.\n"
        ))
      }
      W[!is.finite(W)] <- epsWeights
    }
    if (check) {
      W[, (1:parNum) ^ 2] <- ifelse(
        W[, (1:parNum) ^ 2] < epsWeights, 
        epsWeights, 
        W[, (1:parNum) ^ 2]
      )
    }
    
    
    err <- tryCatch(
      expr = {
        z <- eta + Zfun(eta = eta, weight = W, y = dependent)
        FALSE
      },
      error = function (e) {
        TRUE
      }
    )
    
    if (isTRUE(err)) {
      stop(paste0(
        "Pseudo residuals of IRLS algorithm could not have been computed at iteration: ",
        iter, "\nMost likely working weight matrixes could not have been inverted."
      ), call. = FALSE)
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
    eta <- matrix(eta, ncol = parNum)
    
    LPrev <- L
    L <- -logLike(beta)
    
    if ((iter - 1) %% printOften == 0) {
      if (trace > 0) {cat(sep = "", "Iteration number ", iter, 
                          " log-likelihood: ", 
                          format(L, scientific = FALSE, 
                                 digits = dg))}
      if (trace > 1) {cat(sep = " ", "\nParameter vector: ", 
                          format(beta, scientific = FALSE, 
                                 digits = dg))}
      if (trace > 2) {cat(sep = " ", "\nlog-likelihood reduction: ", 
                          format(L - LPrev, scientific = FALSE, 
                                 digits = dg))}
      if (trace > 3) {cat(sep = " ", "\nValue of gradient at current step:\n", 
                          format(grad(beta), scientific = FALSE, 
                                 digits = dg))}
      eval(traceGreaterThanFourMessegeExpr)
      eval(addToLog)
    }

    if (isTRUE(L < LPrev) || is.infinite(L) || is.nan(L)) {
      halfstepsizing <- TRUE
      h <- step <- stepsize * (betaPrev - beta)
      if ((trace > 0) && (iter - 1) %% printOften == 0) {
        cat("\nTaking a modified step....\n")
      }
      repeat {
        h <- step <- h / 2
        beta <- betaPrev - h
        L <- -logLike(beta)
        if (isTRUE(L > LPrev) && is.finite(L)) {
          break
        }

        if (isTRUE(max(abs(h)) < .Machine$double.eps / 10^6)) {
          halfstepsizing <- FALSE
          if (isTRUE(L < LPrev)) {
            if (!silent) {
              warning("IRLS half-stepping terminated because the step is too small.")
            }
            L <- LPrev
            beta <- betaPrev
          } else {
            if (!silent) {
              warning("IRLS half-stepping terminated because change to score function was too small.")
            }
          }
          break
        }
      }
      if ((iter - 1) %% printOften == 0) {
        if (trace > 0) {cat(sep = "", "Iteration number ", iter, 
                            " log-likelihood: ", 
                            format(L, scientific = FALSE, 
                                   digits = dg))}
        if (trace > 1) {cat(sep = " ", "\nParameter vector: ", 
                            format(beta, scientific = FALSE, 
                                   digits = dg))}
        if (trace > 2) {cat(sep = " ", "\nlog-likelihood reduction: ", 
                            format(L - LPrev, scientific = FALSE, 
                                   digits = dg))}
        if (trace > 3) {cat(sep = " ", "\nValue of gradient at current step:\n", 
                            format(grad(beta), scientific = FALSE, 
                                   digits = dg))}
        eval(traceGreaterThanFourMessegeExpr)
        eval(addToLog)
      }
    }
    
    eta <- covariates %*% beta
    eta <- matrix(eta, ncol = parNum)
    
    if (trace > 0 && (iter - 1) %% printOften == 0) {cat(sep = "", "\n----\n")}
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
  
  if (trace > 4) {
    hhh <- family$makeMinusLogLike(y = dependent, X = covariates, weight = prior, deriv = 2)(beta)
    cat("Value of analytically computed hessian at fitted regression coefficients:\n")
    print(hhh)
    cat("The matrix above has the following eigen values:\n", 
        eigen(hhh, only.values = TRUE)$values, "\n", sep = " ")
    
    if (isTRUE(saveLog)) attr(logg, "hessian") <- hhh
  }
  
  if (isTRUE(saveLog)) {
    if (trace == 1) colnames(logg) <- c("iterationNumber", "halfStep", 
                                        "Log-likelihood")
    if (trace == 2) colnames(logg) <- c("iterationNumber", "halfStep", 
                                        "Log-likelihood", colnames(covariates))
    if (trace  > 2) colnames(logg) <- c("iterationNumber", "halfStep", 
                                        "Log-likelihood", colnames(covariates), 
                                        paste0("gradient -- ", colnames(covariates)))
  }
  
  mu <- mu.eta(eta = eta, ...)
  if (!validmu(mu)) {
    stop("Fit error infinite values reached consider another model,
          mu is too close to zero/infinity")
  }

  list(coefficients = beta, iter = iter, weights = W, logg = logg)
}
# make Xvlm matrix
#' @importFrom stats terms
singleRinternalGetXvlmMatrix <- function(X, formulas, parNames, contrasts = NULL) {
  if (length(formulas[[1]]) == 3) {
    formulas[[1]][[2]] <- NULL
  }
  if (attr(attr(X, "terms"), "response") != 0) {
    X <- X[, colnames(X)[-attr(attr(X, "terms"), "response")], drop = FALSE]
  }
  nPar <- length(parNames)
  Xses <- list()
  
  for (k in 1:nPar) {
    # TODO:: Add contrasts here
    if (length(attr(terms(formulas[[k]], data = X), "term.labels")) != 0) {
      if (attr(terms(formulas[[k]], data = X), "intercept") == 0) {
        Xses[[k]] <- model.matrix(
          ~ . - 1,
          data = X[, attr(terms(formulas[[k]], data = X), "term.labels"), drop = FALSE]
        )
      } else {
        Xses[[k]] <- model.matrix(
          ~ .,
          data = X[, attr(terms(formulas[[k]], data = X), "term.labels"), drop = FALSE]
        )
      }
    } else {
      Xses[[k]] <- model.matrix(
        ~ 1,
        X[, attr(terms(formulas[[k]], data = X), "term.labels"), drop = FALSE]
      )
      if (attr(terms(formulas[[k]], data = X), "intercept") == 0)
        warning(paste0(
          "One of formulas for the model has no variable ",
          "and no intercept and will be coerced to ~ 1."
        ))
    }
    if (k != 1) {
      colnames(Xses[[k]]) <- paste0(colnames(Xses[[k]]), ":", parNames[k])
    }
  }
  hwm <- sapply(Xses, ncol)
  
  # TODO:: Low priority but this could be much better 
  # (without the need to allocate memmory for each X in Xses) 
  # and for Xvlm which just stacks them
  Xvlm <- matrix(0, nrow = nPar * nrow(X), ncol = sum(hwm))
  colnames(Xvlm) <- unlist(sapply(X = Xses, FUN = colnames))
  row <- 0
  col <- 0
  for (k in Xses) {
    Xvlm[(row + 1):(row + nrow(k)), (col + 1):(col + ncol(k))] <- as.matrix(k)
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
  list(reg = as.logical(wch1), # which rows for regression
       est = as.logical(wch2), # which rows for estimation
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
# cholFroW <- function(W, prior) {
#   if (NROW(W) != NROW(prior)) 
#     stop(paste0(
#       "Error in estimatePopsize.fit, working ",
#       "weights and prior weights suggest different number of observations."
#     ))
#   L <- list()
#   for (k in 1:NROW(W)) {
#     L[[k]] <- chol(matrix(prior[k] * W[k,], ncol = 2, nrow = 2))
#   }
#   L
# }
# TODO
# MultibyCholW <- function(Z, X, cholW, which) {
#   # zmień na większą liczbe parametrów
#   W <- cholW
#   if (!is.null(Z)) {# 1 2 to par 3 to mieszane
#     return(c(Z[1:NROW(W)]*cholW[,1]+Z[NROW(W)+1:2*NROW(W)]*cholW[,3], Z[1:NROW(W)]*cholW[,2]))
#   } else if (!is.null(X)) {
#     # to samo co wyżej może wykorzystaj istniejące funkcje
#     # to jest explicite używając własnoći dekompozycji choleskiego
#     return(cbind(
#       rbind(cholW[,1]*X[1:(NROW(X)/2),1:which[1]]+cholW[,3]*X[-(1:(NROW(X)/2)),-(1:which[1])],
#             cholW[,2]*X[-(1:(NROW(X)/2)),1:which[1]]),
#       rbind(cholW[,1]*X[1:(NROW(X)/2),-(1:which[1])]+cholW[,3]*X[-(1:(NROW(X)/2)),-(1:which[1])],
#             cholW[,2]*X[-(1:(NROW(X)/2)),-(1:which[1])])
#     ))
#   }
# }
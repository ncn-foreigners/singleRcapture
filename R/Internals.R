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
                                                     modelFrame, 
                                                     cov, offset,
                                                     weightsFlag) {
  hwm <- attr(Xvlm, "hwm")
  siglevel <- control$alpha
  numboot <- control$B
  sc <- qnorm(p = 1 - siglevel / 2)
  
  if (popVar == "analytic") {
    strappedStatistic <- "No bootstrap performed"
    N <- family$pointEst(pw = weights, eta = eta, y = y)
    if (is.null(cov)) {
      
      cov <- switch(control$covType,
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
      Xvlm = if (family$family == "zelterman") X else Xvlm,
      y = y
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
    
    N <- family$pointEst(pw = weights, eta = eta, y = y)
    
    if (control$cores > 1) {
      funBoot <- switch(
        control$bootType,
        "parametric" = parBootMultiCore,
        "semiparametric" = semparBootMultiCore,
        "nonparametric" = noparBootMultiCore
      )
      
      strappedStatistic <- funBoot(
        family = family, formulas = formulas,
        y = y, X = X, hwm = hwm,
        beta = beta, weights = weights, numboot = numboot,
        eta = eta, cores = control$cores,
        method = control$fittingMethod,
        controlBootstrapMethod = control$bootstrapFitcontrol,
        N = N, Xvlm = Xvlm, modelFrame = modelFrame,
        offset = offset, weightsFlag = weightsFlag,
        visT = control$bootstrapVisualTrace
      )
    } else {
      funBoot <- switch(
        control$bootType,
        "parametric" = parBoot,
        "semiparametric" = semparBoot,
        "nonparametric" = noparBoot
      )
      
      strappedStatistic <- funBoot(
        family = family, formulas = formulas,
        y = y, X = X, hwm = hwm,
        beta = beta, weights = weights, numboot = numboot,
        eta = eta, trace = control$traceBootstrapSize,
        visT = control$bootstrapVisualTrace,
        method = control$fittingMethod,
        controlBootstrapMethod = control$bootstrapFitcontrol,
        N = N, Xvlm = Xvlm, modelFrame = modelFrame,
        offset = offset, weightsFlag = weightsFlag
      )
    }

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
# multiparameter IRLS
singleRcaptureinternalIRLSmultipar <- function(dependent,
                                               family,
                                               formulas,
                                               covariates,
                                               etaStart,
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
                                               offset,
                                               ...) {
  dg <- 8
  converged <- FALSE
  
  # Lowering stepsize to about .3 usually helps a great deal in IRLS fitting
  if (!silent && family$family == "zotnegbin" && stepsize == 1) {
    cat("Zero one truncated negative binomial distribution is prone to taking alpha parameter to infinity.",
        "\nConsider lowering stepsize control parameter if fitting fails.\n")
  }
  
  mu.eta    <- family$mu.eta
  validmu   <- family$validmu
  variance  <- family$variance
  famName   <- family$family
  Zfun      <- family$funcZ
  Wfun      <- family$Wfun
  prior     <- as.numeric(weights)
  dependent <- as.numeric(dependent)
  
  iter <- 1
  step <- NULL
  betaPrev <- NULL
  beta <- rep(0, NCOL(covariates))
  
  logLike <- family$makeMinusLogLike(
    y      = dependent, 
    X      = covariates, 
    weight = prior, 
    offset = offset
  )
  
  grad    <- family$makeMinusLogLike(
    y      = dependent, 
    X      = covariates, 
    weight = prior, 
    deriv  = 1, 
    offset = offset
  )
  
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
  L   <- -Inf
  eta <- etaStart
  
  while (!converged & (iter < maxiter)) {
    halfstepsizing <- FALSE
    mu <- mu.eta(eta = eta, ...)
    if (!validmu(mu)) {
      mu <- mu.eta(eta = eta, ...)
      stop(paste0(
        "Fit error infinite values reached consider another model, mu is too close to zero/infinity.\n"
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
        stop(
          "Working weight matrixes at iteration: ",
          iter,
          " could not have been computed."
        )
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
    
    z <- NULL
    try(
      expr = {z <- eta + Zfun(eta = eta, weight = W, y = dependent, prior = prior) - offset}
    )
    
    XbyW     <- singleRinternalMultiplyWeight(X = covariates, W = W)
    A        <- XbyW %*% covariates
    B        <- XbyW %*% as.numeric(z)
    # The code above does basically:
    # A <- t(Xvlm) %*% WW %*% (Xvlm)
    # B <- t(Xvlm) %*% WW %*% (as.numeric(z))
    # but does not allocate memmory to the full weight matrix
    
    if (any(!is.finite(A)) || any(!is.finite(B))) {
      stop("IRLS step could not have been computed at iteration ", iter, call. = FALSE)
    }
    
    betaPrev <- beta
    stepPrev <- step
    step     <- solve(A, B)
    
    if (!is.null(betaPrev)) {
      step <- step - betaPrev
      beta <- betaPrev + stepsize * 
        (step + if ((is.null(stepPrev) | !momentumFactor)) 0 else {
          if (L-LPrev < momentumActivation) momentumFactor * stepPrev else 0
        })
    } else {
      beta <- step
    }
    
    eta <- covariates %*% beta
    eta <- matrix(eta, ncol = parNum) + offset
    
    LPrev <- L
    L     <- -logLike(beta)
    
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
        h    <- step <- h / 2
        beta <- betaPrev - h
        L    <- -logLike(beta)
        if (isTRUE(L > LPrev) && is.finite(L)) {
          break
        }

        if (isTRUE(max(abs(h)) < .Machine$double.eps / 10^6)) {
          halfstepsizing <- FALSE
          if (isTRUE(L < LPrev)) {
            if (!silent) {
              warning("IRLS half-stepping terminated because the step is too small.")
            }
            L    <- LPrev
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
    eta <- matrix(eta, ncol = parNum) + offset
    
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
    hhh <- family$makeMinusLogLike(y = dependent, X = covariates, weight = prior, deriv = 2, offset = offset)(beta)
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
    stop("Fit error infinite values reached consider another model, mu is too close to zero/infinity")
  }

  list(coefficients = beta, iter = iter, weights = W, logg = logg)
}
#' @importFrom stats terms
singleRinternalGetXvlmMatrix <- function(X, formulas, parNames, contrasts = NULL) {
  nPar <- length(parNames)
  Xses <- list()
  
  for (k in 1:nPar) {
    # Contrasts could be added in this place
    if (length(attr(terms(formulas[[k]], data = X), "term.labels")) != 0) {
      Xses[[k]] <- tryCatch(
        expr = {model.matrix(
          terms(formulas[[k]], data = X),
          data = X
        )},
        error = function (e) {
          ff <- formulas[[k]]
          ff[[2]] <- NULL
          model.matrix(
            terms(ff, data = X),
            data = X
          )
        }
      )
    } else {
      Xses[[k]] <- model.matrix(
        ~ 1,
        X
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
#' @importFrom stats reformulate
singleRinternalMergeFormulas <- function(ff) {
  # This code was inspired by: https://stevencarlislewalker.wordpress.com/2012/08/06/merging-combining-adding-together-two-formula-objects-in-r/
  env <- environment(ff[[1]])
  for (k in 1:length(ff)) { 
    # eliminate left hand side
    if (length(ff[[k]]) == 3) {
      resp <- ff[[k]][[2]]
      ff[[k]][[2]] <- NULL
    }
  }
  dotCheck <- sapply(ff, FUN = function(x) {x == ~.})
  if (any(dotCheck)) {
    out <- reformulate(".", resp)
  } else {
    for (k in 1:length(ff)) { 
      # extract right hand side and create string vector
      ff[[k]] <- strsplit(deparse(ff[[k]][[2]]), " \\+ ")[[1]]
    }
    ff <- unlist(ff)
    out <- reformulate(ff, resp)
  }
  
  environment(out) <- env
  out
}
# This is almost certainly an overkill but it supports arbitrary number of linear predictors
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

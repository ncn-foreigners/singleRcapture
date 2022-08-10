#' Function to fit single-source capture-recapture models
#'
#' @param y vector of dependent variables
#' @param X model matrix
#' @param family same as model in estimate_popsize
#' @param control control parameters
#' @param method method of estimation same as in estimate_popsize
#' @param prior.weights vector of weights its the same argument as weights
#' in estimate_popsize
#' @param start start for regression fitting
#' @param dispersion start for dispersion
#' @param omegaTheta start for theta predictors for omega
#' @param ... arguments to pass to other methods
#' @return
#' list of object connected to regression
#' @export
estimate_popsize.fit <- function(y,
                                 X,
                                 Xvlm,
                                 family,
                                 control,
                                 method,
                                 prior.weights,
                                 start,
                                 dispersion,
                                 omegaTheta,
                                 ...) {
  tbgname <- colnames(X)
  X <- as.matrix(X)

  if (grepl("hurdle", family$family)) {
    FITT <- signleRcaptureinternalIRLS(dependent = y,
                                       covariates = X,
                                       eps = control$epsilon,
                                       epsdisp = control$dispEpsilon,
                                       family = singleRcapture::chao(),
                                       maxiter = control$maxiter,
                                       disp = dispersion,
                                       weights = prior.weights,
                                       start = start,
                                       silent = control$silent,
                                       dispGiven = control$dispGiven,
                                       trace = control$verbose,
                                       stepsize = control$stepsize)
    
    gamma <- FITT$coefficients
    
    FITT <- signleRcaptureinternalIRLS(dependent = y,
                                       covariates = X,
                                       eps = control$epsilon,
                                       epsdisp = control$dispEpsilon,
                                       family = chao(),
                                       maxiter = control$maxiter,
                                       disp = dispersion,
                                       weights = prior.weights,
                                       start = start,
                                       silent = control$silent,
                                       dispGiven = control$dispGiven,
                                       trace = control$verbose,
                                       stepsize = control$stepsize)
    
    dispersion <- FITT$disp
    weights <- FITT$weights
    beta <- c(dispersion, FITT$coefficients)
  } else {
    if (method == "robust") {
      
      if (family$family == "ztoipoisson") {
        stop("Robust not yet implemented for one inflated")
      }
      
      if (family$parNum == 1) {
        FittingFunction <- signleRcaptureinternalIRLS
      } else {
        FittingFunction <- signleRcaptureinternalIRLSmultipar
      }
      
      FITT <- FittingFunction(
        dependent = y,
        covariates = X,
        Xvlm = Xvlm,
        eps = control$epsilon,
        epsdisp = control$dispEpsilon,
        family = family,
        maxiter = control$maxiter,
        disp = dispersion,
        omegaTheta = omegaTheta,
        weights = prior.weights,
        start = start,
        silent = control$silent,
        dispGiven = control$dispGiven,
        omegaThetaGiven = control$thetaGiven,
        trace = control$verbose,
        stepsize = control$stepsize
      )
      
      iter <- FITT$iter
      dispersion <- FITT$disp
      omegaTheta <- FITT$omegaTheta
      weights <- FITT$weights
      beta <- FITT$coefficients
    } else if (method == "mle") {
      logLike <- family$makeMinusLogLike(y = y, X = Xvlm, weight = prior.weights)
      grad <- family$makeGradient(y = y, X = Xvlm, weight = prior.weights, lambdaPredNumber = length(start) - 1)
      
      weights <- prior.weights
      methodopt <- control$mleMethod
      
      if (!isFALSE(control$optimPass)) {
        ctrl <- control$optimPass
      } else {
        if (methodopt == "L-BFGS-B") {
          ctrl <- list(
            factr = control$epsilon,
            maxit = control$maxiter,
            trace = control$verbose,
            trace = if (is.numeric(control$trace)) control$trace else 0
          )
        } else {
          ctrl <- list(
            reltol = control$epsilon,
            maxit = control$maxiter,
            trace = control$verbose
          )
        }
      }

      FITT <- stats::optim(fn = logLike,
                           par = start,
                           gr = function(x) -grad(x),
                           method = methodopt,
                           control = ctrl)
      beta <- FITT$par
      iter <- FITT$counts
      if (FITT$convergence == 1) {
        warning("Convergence not obtained in ", control$maxiter, " iterations of mle fitting algorithm", sep = "")
      }
    } else {
      stop("Method not implemented")
    }
  }

  if (is.null(FITT)) {
    stop("fitting error try another model
          (negative binomial models are highly volitile)")
  }

  beta <- as.vector(beta)

  list(
    beta = beta,
    weights = weights,
    iter = iter
  )
}

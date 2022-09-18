#' Function to fit single-source capture-recapture models
#'
#' @param y vector of dependent variables
#' @param X model matrix
#' @param family same as model in estimate_popsize
#' @param hwm a vector containing information on how many columns in X matrix are associated with each linear predictor
#' @param control control parameters
#' @param method method of estimation same as in estimate_popsize
#' @param prior.weights vector of weights its the same argument as weights
#' in estimate_popsize
#' @param start start for regression fitting
#' @param ... arguments to pass to other methods
#' @return
#' list of object connected to regression
#' @export
estimate_popsize.fit <- function(y, X,
                                 family,
                                 hwm,
                                 control,
                                 method,
                                 prior.weights,
                                 start,
                                 ...) {
  tbgname <- colnames(X)
  X <- as.matrix(X)

  if (grepl("hurdle", family$family)) {
    FITT <- singleRcaptureinternalIRLS(dependent = y,
                                       covariates = X,
                                       eps = control$epsilon,
                                       family = singleRcapture::chao(),
                                       maxiter = control$maxiter,
                                       weights = prior.weights,
                                       start = start,
                                       silent = control$silent,
                                       trace = control$verbose,
                                       stepsize = control$stepsize)
    
    gamma <- FITT$coefficients
    
    FITT <- singleRcaptureinternalIRLS(dependent = y,
                                       covariates = X,
                                       eps = control$epsilon,
                                       family = chao(),
                                       maxiter = control$maxiter,
                                       weights = prior.weights,
                                       start = start,
                                       silent = control$silent,
                                       trace = control$verbose,
                                       stepsize = control$stepsize)
    
    weights <- FITT$weights
    beta <- FITT$coefficients
  } else {
    if (method == "robust") {
      
      if (family$parNum == 1) {
        FittingFunction <- singleRcaptureinternalIRLS
      } else {
        FittingFunction <- singleRcaptureinternalIRLSmultipar
      }
      
      FITT <- FittingFunction(
        dependent = y,
        covariates = X,
        eps = control$epsilon,
        family = family,
        maxiter = control$maxiter,
        weights = prior.weights,
        start = start,
        silent = control$silent,
        trace = control$verbose,
        stepsize = control$stepsize,
        hwm = hwm,
        momentumFactor = control$momentumFactor,
        momentumActivation = control$momentumActivation,
      )
      
      iter <- FITT$iter
      weights <- FITT$weights
      beta <- FITT$coefficients
    } else if (method == "mle") {
      logLike <- family$makeMinusLogLike(y = y, X = X, weight = prior.weights)
      grad <- family$makeGradient(y = y, X = X, weight = prior.weights, lambdaPredNumber = length(start) - 1)
      
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
                           # hessian = TRUE,
                           control = ctrl)
      beta <- FITT$par
      # print(FITT)
      # print(rootSolve::gradient(f = grad, x = beta))
      # Commented lines of code are used in verification of computed analytic hessian
      iter <- FITT$counts
      if (FITT$convergence == 1 && !control$silent) {
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

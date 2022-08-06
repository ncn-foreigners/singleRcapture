#' Function to fit single-source capture-recapture models
#'
#' @param y Vector of dependent variables.
#' @param X A model matrix.
#' @param family Same as model in \code{estimate_popsize}.
#' @param control Control parameters.
#' @param method Method of estimation same as in \code{estimate_popsize}.
#' @param prior.weights Vector of weights its the same argument as weights
#' in \code{estimate_popsize}.
#' @param start Start for regression fitting.
#' @param dispersion Start for dispersion.
#' @param ... Arguments to pass to other methods.
#' @return
#' list of object connected to regression
#' @export
estimate_popsize.fit <- function(y,
                                 X,
                                 family,
                                 control,
                                 method,
                                 prior.weights,
                                 start,
                                 dispersion,
                                 ...) {
  tbgname <- colnames(X)
  X <- as.matrix(X)

  if (method == "robust") {

    if (!is.null(dispersion)) {
      start <- start[-1]
    }

    FITT <- signleRcaptureinternalIRLS(dependent = y,
                                       covariates = X,
                                       eps = control$epsilon,
                                       family = family,
                                       maxiter = control$maxiter,
                                       disp = dispersion,
                                       weights = prior.weights,
                                       start = start,
                                       silent = control$silent,
                                       disp.given = control$disp.given,
                                       trace = control$verbose)

    iter <- FITT$iter
    dispersion <- FITT$disp
    weights <- FITT$weights
    beta <- c(dispersion, FITT$coefficients)
  } else if (method == "mle") {
    log_like <- family$make_minusloglike(y = y, X = X,
                                         weight = prior.weights)
    grad <- family$make_gradient(y = y, X = X,
                                 weight = prior.weights)
    
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

    FITT <- stats::optim(fn = log_like,
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

  if (is.null(FITT)) {
    return(NULL)
    stop("fitting error try another model
          (negative binomial models are highly volitile)")
  }

  beta <- as.vector(beta)

  if (is.null(dispersion)) {
    names(beta) <- tbgname
  } else {
    names(beta) <- c("log(dispersion)", tbgname)
  }

  list(
    beta = beta,
    weights = weights,
    iter = iter
  )
}

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
#' @param ... arguments to pass to other methods
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

  log_like <- family$make_minusloglike(y = y, X = X,
                                       weight = prior.weights)
  grad <- family$make_gradient(y = y, X = X,
                               weight = prior.weights)
  hessian <- family$make_hessian(y = y, X = X,
                                 weight = prior.weights)

  df.reduced <- length(y) - dim(X)[2]

  if (family$family %in% c("zotnegbin", "ztnegbin")) {
    df.reduced <- 2 * df.reduced
  }

  if (method == "robust") {

    if (!is.null(dispersion)) {
      start <- start[-1]
    }

    FITT <- IRLS(dependent = y,
                 covariates = X,
                 eps = control$epsilon,
                 family = family,
                 maxiter = control$maxiter,
                 disp = dispersion,
                 weights = prior.weights,
                 start = start,
                 silent = control$silent,
                 disp.given = control$disp.given,
                 trace = control$trace)

    iter <- FITT$iter
    dispersion <- FITT$disp
    weights <- FITT$weights
    beta <- c(dispersion, FITT$coefficients)
  } else if (method == "mle") {
    weights <- prior.weights
    methodopt <- control$mleMethod
    
    if (!isFALSE(control$optimPass)) {
      ctrl <- control$optimPass
    } else {
      ctrl <- list(
        factr = control$epsilon,
        maxit = control$maxiter,
        trace = if (is.numeric(control$trace)) control$trace else 0
      )
    }
    
    FITT <- stats::optim(fn = log_like,
                         par = start,
                         gr = function(x) -grad(x),
                         method = methodopt,
                         control = ctrl)
    beta <- FITT$par
    iter <- FITT$counts
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

  hess <- hessian(beta)

  list(ll = log_like,
       grad = grad,
       hessian = hessian,
       beta = beta,
       weights = weights,
       hess = hess,
       iter = iter,
       degf = df.reduced)
}

#' Regression fitting for singleRcapture
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
                 eps = .Machine$double.eps,
                 family = family,
                 disp = dispersion,
                 weights = prior.weights,
                 start = start)

    iter <- FITT$iter
    dispersion <- FITT$disp
    weights <- FITT$weights
    beta <- c(dispersion, FITT$coefficients)
  } else if (method == "mle") {
    weights <- prior.weights
    methodopt <- "L-BFGS-B"
    
    ctrl <- list(factr = .Machine$double.eps,
                 maxit = 5000)
    
    if (family$family == "ztnegbin") {
      methodopt <- "Nelder-Mead"
      ctrl <- list(reltol = .Machine$double.eps,
                   maxit = 5000)
    }
    
    if (dim(X)[2] == 1) {
      FITT <- stats::optim(par = start,
                           fn = log_like,
                           gr = function(x) -grad(x),
                           method = "Nelder-Mead",
                           control = list(reltol = .Machine$double.eps,
                                          warn.1d.NelderMead = FALSE))
      beta <- FITT$par
      iter <- FITT$counts
    } else {
      FITT <- stats::optim(fn = log_like,
                           par = start,
                           gr = function(x) -grad(x),
                           method = methodopt,
                           control = ctrl)
      beta <- FITT$par
      iter <- FITT$counts
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

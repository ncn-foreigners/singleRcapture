#' Iteratively Reweighted Least Squares
#'
#' A method for fitting regression parameters in generalized linear model
#'
#' @param dependent A vector of Dependent variables
#' @param covariates An array of Covariates
#' @param start A place to start the numerical estimation
#' @param maxiter Maximum number of iterations
#' @param disp dispersion parameter if family needs it
#' @param eps Precision level
#' @param family Family of distributions used in regression
#' @param weights Optional object of weights used in fitting the model
#' @param disp.given FALSE is default, set true if dispersion has already been
#' estimated and there is no need for further estimation
#' @return Returns a list of fitted Coefficients and number of iterations and final weights
#' @importFrom stats optim
#' @export
IRLS <- function(dependent,
                 family,
                 covariates,
                 start,
                 disp = NULL,
                 weights = NULL,
                 maxiter = 1000,
                 eps = .Machine$double.eps ** .25,
                 disp.given = FALSE) {
  converged <- FALSE
  iter <- 1
  beta <- start

  if (!is.null(weights)) {
    prior <- as.numeric(weights)
  } else {
    prior <- 1
  }
  W <- prior

  loglike <- family$make_minusloglike(y = dependent,
                                      X = covariates,
                                      weight = prior)
  grad <- family$make_gradient(y = dependent,
                               X = covariates,
                               weight = prior)
  temp <- c(disp, beta)
  L <- -loglike(temp)
  dispPrev <- Inf

  while (!converged && (iter < maxiter)) {
    if (family$family %in% c("ztnegbin", "zotnegbin") &&
        isFALSE(disp.given) && (abs(disp - dispPrev) > 1e-5)) {
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

    WPrev <- W
    tempPrev <- temp
    betaPrev <- beta
    LPrev <- L

    eta <- covariates %*% beta
    mu <- family$mu.eta(eta = eta, disp)
    if (!family$validmu(mu)) {
      stop("Fit error infinite values reached consider another model,
            mu is too close to zero/infinity")
    }

    varY <- family$variance(mu = mu, disp)
    Z <- eta + (dependent - mu) / mu
    W <- as.numeric(prior * (mu ** 2) / varY)
    # This is equivalent to
    # A <- t(covariates) %*% W %*% covariates
    # B <- t(covariates) %*% W %*% Z
    # But much much faster and less memory heavy
    A <- t(covariates) %*% (covariates * W)
    B <- t(covariates) %*% (Z * W)
    beta <- solve(A, B, tol = .Machine$double.eps)

    temp <- c(disp, beta)
    L <- -loglike(temp)

    converged <- (((L - LPrev) < eps) || (max(abs(beta - betaPrev)) < eps))

    if (!converged) {
      iter <- iter + 1
    } else if ((L - LPrev) < 0) {
      beta <- betaPrev
      L <- LPrev
      W <- WPrev
    }

  }

  list(coefficients = beta, iter = iter, weights = W, disp = disp)
}

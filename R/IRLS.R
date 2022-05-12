#' Iteratively Reweighted Least Squares
#'
#' A method for fitting regression parameters in generalized linear model. \cr
#' Used internally in singleRcapture package.
#'
#' @param dependent A vector of dependent variables
#' @param covariates An array of covariates
#' @param start Vector of starting points for the the numerical procedure
#' @param maxiter Maximum number of iterations
#' @param disp dispersion parameter if family needed
#' @param eps Tolerance level for the numerical procedure
#' @param family Family of distributions used in regression
#' @param weights Optional object of weights used in fitting the model
#' @param disp.given FALSE is default, set true if dispersion has already been
#' estimated and there is no need for further estimation
#' @param silent Boolean value indicating whether warnings should be suppressed
#' @return Returns a list of fitted Coefficients and number of iterations and final weights
#' @importFrom stats optim
#' @export
IRLS <- function(dependent,
                 family,
                 covariates,
                 start,
                 disp = NULL,
                 weights = NULL,
                 maxiter = 10000,
                 eps = .Machine$double.eps,
                 disp.given = FALSE,
                 silent = FALSE) {
  converged <- FALSE
  epsdisp <- 1e-5 # TODO add to controll
  
  mu.eta <- family$mu.eta
  validmu <- family$validmu
  variance <- family$variance
  famName <- family$family
  linkinv <- family$linkinv
  funcZ <- function(mu, y, eta) {eta + (y - mu) / mu}
  Wfun <- function(mu, prior, varY, eta) {prior * (mu ** 2) / varY}

  if (famName %in% c("chao", "zelterman")) {
    dependent <- dependent - 1
    # Modification to functions in IRLS to account for logit link in fitting
    linkinv <- function(p) {1 / (1 + exp(-p))}
    funcZ <- function(mu, y, eta) {eta + (y - mu) / variance(mu = mu, disp)}
    Wfun <- function(mu, prior, varY, eta) {prior * varY}
  }

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
    Z <- funcZ(mu = mu, y = dependent, eta = eta)
    W <- as.numeric(Wfun(mu = mu, prior = prior, 
                         varY = varY, eta = eta))
    # This is equivalent to
    # A <- t(covariates) %*% W %*% covariates
    # B <- t(covariates) %*% W %*% Z
    # But much much faster and less memory heavy
    A <- t(covariates) %*% (covariates * W)
    B <- t(covariates) %*% (Z * W)
    beta <- solve(A, B, tol = .Machine$double.eps)

    temp <- c(disp, beta)
    L <- -loglike(temp)

    #print(L - LPrev)
    converged <- (((L - LPrev) < eps) || (max(abs(beta - betaPrev)) < eps))

    if (!converged) {
      iter <- iter + 1
    } else if ((L - LPrev) < 0) {
      beta <- betaPrev
      L <- LPrev
      W <- WPrev
    }

    if(iter == maxiter && !converged && !silent) {
      warning("Fitting algorithm (IRLS) has not converged")
    }

  }

  list(coefficients = beta, iter = iter, weights = W, disp = disp)
}

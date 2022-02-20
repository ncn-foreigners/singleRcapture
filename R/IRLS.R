#' Iteratively Reweighted Least Squares
#'
#' A method for fitting regression parameters in generalized linear model
#'
#' @param Dependent A vector of Dependent variables
#' @param Covariates An array of Covariates
#' @param start A place to start the numerical estimation
#' @param Maxiter Maximum number of iterations
#' @param disp dispersion parameter if family needs it
#' @param eps Precision level
#' @param family Family of distributions used in regression
#' @param weights Optional object of weights used in fitting the model
#' @param disp.given FALSE is default, set true if dispersion has already been
#' estimated and there is no need for further estimation
#' @return Returns a list of fitted Coefficients and number of iterations and final weights
#' @export
IRLS <- function(Dependent,
                 family,
                 Covariates,
                 start,
                 disp = NULL,
                 weights = NULL,
                 Maxiter = 10000,
                 eps = 1e-10,
                 disp.given = FALSE) {
  converged <- FALSE
  Iter <- 1
  Beta <- start
  W <- NULL

  if (!is.null(weights)) {
    Prior <- as.numeric(weights)
  } else {
    Prior <- 1 / length(Dependent)
  }

  Loglike <- family$make_minusloglike(y = Dependent,
                                      X = Covariates,
                                      weight = Prior)
  Grad <- family$make_gradient(y = Dependent,
                               X = Covariates,
                               weight = Prior)

  if (family$family %in% c("ztnegbin", "zotnegbin")) {
    Temp <- c(disp, Beta)
  } else {
    Temp <- Beta
  }
  L <- -Loglike(Temp)
  dispPrev <- Inf

  while (!converged && (Iter < Maxiter)) {
    # Think about making your own optim method
    if (family$family %in% c("ztnegbin", "zotnegbin") &&
        isFALSE(disp.given) &&
        (abs(disp - dispPrev) > eps)) {
      dispPrev <- disp
      ll <- function(alpha) Loglike(c(alpha, Beta))
      gr <- function(alpha) Grad(c(alpha, Beta))[1]
      disp <- stats::optimize(f = ll,
                              interval = c(disp - disp / 2,
                                           disp + disp / 2))$minimum
    }

    PrevWeight <- W
    PrevBeta <- Beta
    PrevL <- L

    Eta <- Covariates %*% Beta
    Parameter <- family$linkinv(Eta)
    mu <- family$mu.eta(eta = Eta, disp)
    if (!family$validmu(mu)) {
      stop("Fit error infinite values reached consider another model,
            mu is too close to zero/infinity")
    }

    VarY <- family$variance(mu = Parameter, disp)
    Z <- Eta + (Dependent - mu) * family$Dlink(mu)
    W <- as.numeric(Prior / ((family$Dlink(mu) ** 2) * VarY))
    # This is equivalent to
    # A <- t(Covariates) %*% W %*% Covariates
    # B <- t(Covariates) %*% W %*% Z
    # But much much faster and less memory heavy
    A <- t(Covariates) %*% (Covariates * W)
    B <- t(Covariates) %*% (Z * W)
    Beta <- solve(A, B)

    if (family$family %in% c("ztnegbin", "zotnegbin")) {
      Temp <- c(disp, Beta)
    } else {
      Temp <- Beta
    }
    L <- -Loglike(Temp)

    converged <- (((L - PrevL) < eps) || (max(abs(Beta - PrevBeta)) < eps))

    if (!converged) {
      Iter <- Iter + 1
    } else if ((L - PrevL) < 0) {
      Beta <- PrevBeta
      L <- PrevL
      W <- PrevWeight
    }

  }

  list(Coefficients = Beta, iter = Iter, Weights = W, disp = disp)
}

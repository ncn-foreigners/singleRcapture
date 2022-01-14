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
#'
#' @return Returns a list of fitted Coefficients and number of iterations and final weights
#' @export
IRLS <- function(Dependent,
                 family,
                 Covariates,
                 start,
                 disp = NULL,
                 weights = NULL,
                 Maxiter = 10000,
                 eps = 1e-10) {
  converged <- FALSE
  Iter <- 1
  Beta <- start

  if (!is.null(weights)) {
    Prior <- as.numeric(weights)
  } else {
    Prior <- 1
  }

  Loglike <- family$make_minusloglike(y = Dependent,
                                      X = Covariates,
                                      weight = Prior)
  if (family$family == "ZTNB") {
    Temp <- c(disp, Beta)
  } else {
    Temp <- Beta
  }
  L <- -Loglike(Temp)


  while (!converged && (Iter < Maxiter)) {
    PrevBeta <- Beta
    PrevL <- L

    Eta <- Covariates %*% Beta
    Parameter <- family$linkinv(Eta)
    mu <- family$mu.eta(eta = Eta, disp)
    VarY <- family$variance(mu = mu, disp)
    Z <- Eta + (Dependent - mu) * family$Dlink(mu)
    W <- as.numeric(Prior / ((family$Dlink(mu) ** 2) * VarY))

    # This is equivalent to
    # A <- t(Covariates) %*% W %*% Covariates
    # B <- t(Covariates) %*% W %*% Z
    # But much much faster and less memory heavy
    A <- t(Covariates) %*% (Covariates * W)
    B <- t(Covariates) %*% (Z * W)
    Beta <- solve(A, B)

    if (family$family == "ZTNB") {
      Temp <- c(disp, Beta)
    } else {
      Temp <- Beta
    }
    L <- -Loglike(Temp)

    converged <- ((L-PrevL) < eps) || (max(abs(Beta - PrevBeta)) < eps)

    if (!converged) {
      Iter <- Iter + 1
    }

  }
  list("Coefficients" = Beta, "iter" = Iter, Weights = W)
}

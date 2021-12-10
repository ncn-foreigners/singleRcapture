#' Itterativly Reweighted Least Squares
#'
#' A method for fitting regression parameters in generalised linear model
#'
#' @param Dependent A vector of Dependent variables
#' @param Covariates An array of Covariates
#' @param start A place to start the numerical estimation
#' @param Maxiter Maximum number of itterations
#' @param eps Precision level
#' @param Family Family of distributions used in regression
#'
#' @return Returns a list of fitted Coefficients and number of itterations
#' @export
IRLS <- function(Dependent,Family,Covariates,start,Maxiter = 100,eps = 1e-10){
  converged <- FALSE
  Iter <- 1
  Beta <- start
  Loglike <- Family()$make_minusloglike(y = Dependent,X = Covariates)

  Eta <- Covariates %*% Beta
  Parameter <- Family()$Invlink(Eta)
  mu <- Family()$mu(Parameter)
  L <- -Loglike(Beta)

  while(!converged && (Iter < Maxiter)){
    PrevBeta <- Beta
    PrevL <- L

    Eta <- Covariates %*% Beta
    Parameter <- Family()$Invlink(Eta)
    mu <- Family()$mu(Parameter)

    Z <- Covariates %*% Beta + (Dependent-mu) / mu
    W <- diag(as.numeric(mu))

    A <- (t(Covariates) %*% W) %*% Covariates
    B <- (t(Covariates) %*% W) %*% Z
    Beta <- solve(A,B)

    L <- -Loglike(Beta)

    converged <- ((L-PrevL) < eps) || (max(abs(Beta - PrevBeta)) < eps)

    if(!converged){
      Iter <- Iter + 1
    }

  }
  return(list('Coefficients' = Beta,'Number_of_itterations' = Iter))
}

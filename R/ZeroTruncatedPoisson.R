#' Zero Truncated Poisson Model Functions
#'
#' Used internally withing a package creates gradient hessian and minus log likelihood
#' for Zero Truncated Poisson
#'
#' @return A list of 3 objects \cr
#' make_minusloglike(y,X) - for creating negative likelihood function \cr
#' make_gradient(y,X) - for creating gradient function \cr
#' make_hessian(X) - for creating hessian \cr
#' Where: \cr
#' y is a vector of observed values \cr
#' X is a matrix / data frame of covariates
#' @export
ZeroTruncatedPoisson <- function(){
  link <- log
  Invlink <- exp

  MinusLogLike <- function(y,X){
    y <- as.numeric(y)
    function(beta){
      Eta <- as.matrix(X) %*% beta
      Lambda <- exp(as.matrix(X) %*% beta)
      -sum(y * Eta - log(exp(Lambda) - 1) - log(factorial(y)))
    }
  }

  Gradient <- function(y,X){
    y <- as.numeric(y)
    function(beta){
      Lambda <- exp(as.matrix(X) %*% beta)
      mu <- Lambda / (1 - exp(-Lambda))
      t(as.matrix(X)) %*% (y - mu)
    }
  }

  Hessian <- function(X){
    function(beta){
      Lambda <- exp(as.matrix(X) %*% beta)
      coefficient <- 1/(1-exp(-Lambda)) - Lambda * exp(-Lambda) / ((1-exp(-Lambda)) ** 2)
      Dmu <- diag(as.numeric(coefficient))
      Dlam <- as.matrix(X * Lambda)

      -( (t(as.matrix(X)) %*% Dmu) %*% Dlam)
    }
  }
  mu <- function(Lambda){
    Lambda / ( 1 - exp(-Lambda))
  }

  return(list(make_minusloglike = MinusLogLike,make_gradient = Gradient,
              make_hessian = Hessian,link = link,Invlink = Invlink,mu = mu))
}

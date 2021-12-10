#' Total Population Point and Interval estimate
#'
#' Creates a point and interval estimate for total and missing population from zero truncated poisson model.
#'
#' @param y Observed values
#' @param X A covariate matrix
#' @param Grad A gradient of a model with respect to regression parameters
#' @param lambda A lambda Parameter for the model
#' @param beta Fitted regression values
#' @param Hess Hessian of a model
#'
#' @return Returns a list of size 3 with: \cr
#' (1) Point estimate \cr
#' (2) Variance \cr
#' (3) Confidence interval constructed using Point estimate plus minus 1.96 times standard deviation
#' @export
#'
PoissonPopEst <- function(y,X,Grad,lambda,beta,Hess){
  X <- as.data.frame(X)
  N <- sum(1/(1-exp(-lambda)))
  Inform <- Hess(beta)

  f1 <- colSums(-X * (exp(log(lambda)-lambda)/((1-exp(-lambda)) ** 2)))
  f1 <- t(f1) %*% solve(as.matrix(Inform)) %*% f1

  f2 <- sum(exp(-lambda)/((1-exp(-lambda)) ** 2))

  VV <- f1 + f2

  return(list(Point_estimate = N,Variation = VV,
              ConfidenceInterval = c(N-1.96*sqrt(VV),N+1.96*sqrt(VV))))
}
